classdef sampHighOrder
% Perform MRI sampling using direct summation of complex exponentials.
% Supports B0 map, time-varying spherical harmonic distributions of phase,
% and time-varying distributions of phase that are typical for concomitant
% gradients. Currently 2D only.
%
% Usage: 
%   Define sampHighOrder object using 
%       S = sampHighOrder(b0,sampTimes,phs_spha,phs_conc,phs_grid,imMask,useGPU,useSingle,useInterp,svdThresh,subFact)
%   Sample image using data = S*image. x can have more dimensions than
%       b0, but summations are only performed over first two dims (i.e., 2D
%       imaging is assumed). Fairly high memory overhead due to 
%       precomputations, but quite fast.
%   Perform adjoint of sampling using S'*data. Required for iterative
%       solvers like lsqr.
%
% Inputs: 
%   b0: B0 map in units of rad/s. Caution: make sure orientation
%       of b0 is consistent with phs_grid.x and phs_grid.y.
%   sampTimes: sampling times in units of s. Can be multi-dimensional.
%   phs_spha:  [Ncoeff_spha x size(sampTimes)] array of coefficients for
%       spherical harmonic phase. Units are rad/<spatial>, where <spatial>
%       can be unitless (DC), m (normal k-space units), m^2, etc, depending
%       on the coefficient.
%   phs_conc:  [Ncoeff_conc x size(sampTimes)] array of coefficients for
%       concomitant grad phase. Units are rad/<spatial>, where <spatial>
%       can be unitless (DC), m (normal k-space units), m^2, etc, depending
%       on the coefficient.
%   phs_grid: struct with x, y, and z positions for each voxel in b0. 
%       Should be created using ndgrid or meshgrid. If created with
%       meshgrid, phs_grid.x will have positions varying along 2nd dim
%       (with ndgrid, variation along 1st dim), and vice versa for
%       phs_grid.y. Each field of phs_grid must be the same size as b0.
%           phs_grid.x: value of x at each position in units of m. 
%           phs_grid.y: value of y at each position in units of m. 
%           phs_grid.z: value of z at each position in units of m. 
%   imMask: mask denoted expected location of signal. b0map and 2nd
%       order+ spherical harmonic phase terms are set to 0 outside the
%       mask. This improves the performance of interpolation (useInterp
%       option), because it removes regions with quickly varying phase that
%       are irrevelant to image recon (since there's no signal there).
%   useGPU: whether to use GPU. Default = true.
%   useSingle: if true, use single instead of double precision. Default =
%       false. 
%   useInterp: uses an interpolated approach. Always faster than direct
%       approach on CPU (i.e., when useGPU = 1). Not much benefit for GPU,
%       unless not enough memory for direct approach. 
%		See Wilm et al DOI: 10.1109/TMI.2012.2190991
%   svdThresh (default = 0.05): trades off accuracy with speed when
%       useInterp=1. Decrease to improve accuracy.
%   subFact (default = 5): another trade-off for accuracy and speed for 
%       useInterp=1. Decrease for accuracy. 5 or less should have
%       negligible error. Min val = 1.
%
%
%   Within the class, basis functions for each index of phs_spha and
%   phs_conc are computed using basisFuncSphHarm and basisFuncConcGrad,
%   respectively. 
%
%   (c) Corey Baron 2020

	properties
		NDim = 2;   % Number of dims to do sums over in both image domain and k-space. 
        useSingle = 0;
		useGPU = 1;
		useInterp = 0; 
        noPrecomp = 0;
        svdThresh = 0.05; % Threshold for svd used in interpolated method
        subFact = 5;      % Factor to subsample in time by for interpolated approach
	end

	properties (SetAccess = protected)
		adjoint = 0;
		b0 = []; % rad/sec
        b0mask = [];
        phs_spha = []; % 0th, 1st, 2nd and 3rd order spherical harmonic terms for phase over time. Must have dimensions 16 x size(sampTimes). 
        phs_conc = []; % conc grad terms for phase over time. Must have dimensions 4 x size(sampTimes). 
        phs_grid = []; % Struct with fields phs_grid.X, phs_grid.Y, phs_grid.Z. Each must have dimensions equivalent to b0. MUST be in magnet frame.
		sampTimes = []; % sec
		kbase = []; % spatial part of spherical harmonics that is multiplied with phs_spha or phs_conc terms
        phiDiv = []; % temporal derivative of kbase. Can be used to find global delays. See DOI: 10.1002/mrm.29460
		traj = [];		  % Precomputed values for interpolated approach
        svdSpace = [];    % Precomputed values for interpolated approach
        svdTime = [];     % Precomputed values for interpolated approach
        phsShft = [];
		kSize = [];
		imSize = [];
        trajFromRaw = [];
        noPrecompNdiv = [];
        ksphaDiv = [];
        kconcDiv = [];
	end

	methods

		function obj = sampHighOrder(b0,sampTimes,phs_spha,phs_conc,phs_grid,b0mask,useGPU,useSingle,useInterp,svdThresh,subFact,noPrecomp)
			if nargin == 0
				obj.tests;
				return;
            end
            if isstruct(b0)
                % Overloaded for computing trajectory from raw
                obj.trajFromRaw = obj.computeTraj(b0.probe_positions,b0.probe_raw,b0.gammaProbes,b0.gammaMRI,b0.dt,b0.B0,b0.fieldOffsets,b0.coilParams,b0.nonLinSphHarm,b0.fitOrder);
                return;
            end
            if nargin>5
				obj.b0mask = b0mask;
            end
            if nargin<6 || isempty(useGPU)
                % Determine options based on memory available
                % Precomputations using GPU requires about 3.3 times numel(b0)*numel(sampTimes)*8 bytes of memory 
                % After this calc, sampHighOrder object requires numel(b0)*numel(sampTimes)*8 bytes
                % If the calc is successful, there is enough working memory to run
                % multiplications with sampHighOrder object.
                % To be safe, we allow for 4 times numel(b0)*numel(sampTimes)*8
                if ~(gpuDeviceCount>0)
                    obj.useGPU = 0;
                    obj.useInterp = 1;
                else
                    G = gpuDevice;
                    avMem = G.AvailableMemory;
                    numbytes = numel(b0)*numel(sampTimes)*8;
                    if numbytes < 0.2*avMem
                        obj.useGPU = 1;
                        obj.useSingle = 0;
                        obj.useInterp = 0;
                    elseif numbytes < 0.4*avMem
                        obj.useGPU = 1;
                        obj.useSingle = 1;
                        obj.useInterp = 0;
                    else
                        obj.useGPU = 1;
                        obj.useSingle = 1;
                        obj.useInterp = 0;
                        obj.noPrecomp = 1;
                        obj.noPrecompNdiv = ceil(numbytes / (0.2*avMem)); % 0.2 factor was determine heuristically
                        if obj.noPrecompNdiv > max(numel(b0),numel(sampTimes))/4
                            % Unlikely, but possible. Still use single,
                            % because matrix size must be massive.
                            obj.useGPU = 0;
                            obj.useSingle = 1;
                            obj.useInterp = 1;
                            obj.noPrecomp = 0;
                        end
                    end
                end
            else
                if nargin>6 && ~isempty(useGPU)
                    obj.useGPU = useGPU;
                end
                if nargin>7 && ~isempty(useSingle)
                    obj.useSingle = useSingle;
                end
                if nargin>8 && ~isempty(useInterp)
                    obj.useInterp = useInterp;
                end
                if nargin>9 && ~isempty(svdThresh)
                    obj.svdThresh = svdThresh;
                end
                if nargin>10 && ~isempty(subFact)
                    obj.subFact = subFact;
                end
                if nargin>11 && ~isempty(noPrecomp)
                    obj.noPrecomp = noPrecomp;
                    obj.noPrecompNdiv = 2;
                end
            end
            if obj.subFact < 1
                obj.subFact = 1;
            end
            % Check for GPU support
            if ~(gpuDeviceCount>0) && obj.useGPU
                warning('No GPU detected or GPU not supported. Using CPU.')
                obj.useGPU = 0;
            end
			obj.NDim = 2; % This class currently coded/optimized for 2D only
			obj.b0 = b0;    
			obj.imSize = size(b0);   
			obj.sampTimes = sampTimes;
            obj.kSize = size(obj.sampTimes);
			obj.phs_spha = phs_spha;
			sz_a = size(phs_spha);
            if length(obj.kSize) == 2 && obj.kSize(2) == 1
                % Account for 1D sampTimes
                sz_a = [sz_a, 1];
            end
			if any(sz_a(2:end) ~= size(sampTimes))
				error('Dimension mismatch between phs_spha and sampTimes')
			end
			if length(obj.kSize)>obj.NDim || length(obj.imSize)>obj.NDim
				error('Only up to %d dimensions allowed',obj.NDim);
            end
            % Check if gpu is possible
            if ~(gpuDeviceCount>0) && obj.useGPU
                warning('No GPU detected. Using CPU. To disable this warning, set option useGPU to 0')
                obj.useGPU = 0;
            end
			obj.phs_conc = phs_conc;
            obj.phs_grid = phs_grid;
            if obj.NDim < 2
                % This is needed because size() always outputs at least two
                % values, even for vectors. Can hack 1D by having all size
                % of b0 = 1 x N.
                error('need at least 2 dimensions')
            end
			% Use single precision if requested
			if obj.useSingle
                obj.sampTimes = single(obj.sampTimes);
				obj.phs_spha = single(obj.phs_spha);
				obj.phs_conc = single(obj.phs_conc);
				obj.b0 = single(obj.b0);
				obj.phs_grid.x = single(obj.phs_grid.x);
				obj.phs_grid.y = single(obj.phs_grid.y);
				obj.phs_grid.z = single(obj.phs_grid.z);
                obj.b0mask = single(obj.b0mask);
			end
            % Move variables to GPU
			if obj.useGPU && ~obj.useInterp
				obj.sampTimes = gpuArray(obj.sampTimes);
				obj.phs_spha = gpuArray(obj.phs_spha);
				obj.phs_conc = gpuArray(obj.phs_conc);
				obj.b0 = gpuArray(obj.b0);
				obj.phs_grid.x = gpuArray(obj.phs_grid.x);
				obj.phs_grid.y = gpuArray(obj.phs_grid.y);
				obj.phs_grid.z = gpuArray(obj.phs_grid.z);
                obj.b0mask = gpuArray(obj.b0mask);
			end
			if obj.useInterp
                [obj.svdTime,obj.svdSpace,obj.traj,obj.phsShft] = prepForInterp(obj);
                if obj.useGPU
                    obj.svdTime = gpuArray(obj.svdTime);
                    obj.svdSpace = gpuArray(obj.svdSpace);
                    obj.phsShft = gpuArray(obj.phsShft);
                end
            elseif ~obj.noPrecomp
				obj.kbase = prepForDirect(obj,obj.phs_spha,obj.phs_conc,obj.sampTimes);
			end
		end

		function y = mtimes(obj,x)
			szx = [size(x), 1, 1];
            if obj.useSingle
                if (isa(x,'gpuArray') && ~isaUnderlying(x,'single')) || ~isa(x,'single')
                    %warning('useSingle specified, but input is not single. Forcing to be single.')
                    x = single(x);
                end
            end

			if obj.adjoint
				y = zeros([obj.imSize, szx(3:end)], 'like', x);
			else
				y = zeros([obj.kSize, szx(3:end)], 'like', x);
			end
			for n=1:prod(szx(3:end))
				x_a = x(:,:,n);
				if obj.useGPU && ~isa(x_a, 'gpuArray')
					x_a = gpuArray(x_a);
				end
				y_a = applyModel(obj,x_a);
				if ~isa(x,'gpuArray')
					y_a = gather(y_a);
				end
				y(:,:,n) = y_a;
            end
        end
		
		function y = applyModel(obj,x)
			if obj.useInterp
				% Use nufft's and interpolation
				y = 0;
				if obj.adjoint
					for l = 1:size(obj.svdSpace,2)
                        y_a = x.*conj(reshape(obj.svdTime(:,l), size(obj.sampTimes))); 
                        y_a = conj(obj.phsShft).*y_a;
						y_a = obj.traj'*y_a;
                        y = y + y_a.*conj(reshape(obj.svdSpace(:,l),size(obj.b0)));
					end
				else
                    if ~isempty(obj.ksphaDiv)
                        error('phiDiv not yet implemented for interpolated approach')
                    end
					for l = 1:size(obj.svdSpace,2)
                        y_a = x.*reshape(obj.svdSpace(:,l),size(obj.b0));
                        y_a = obj.traj*y_a;
                        y_a = obj.phsShft.*y_a;
                        y = y + y_a.*reshape(obj.svdTime(:,l), size(obj.sampTimes)); 
					end
				end
			else
				% Use direct model
                if obj.noPrecomp
                    y = noPrecompWorker(obj,x);
                else
                    if obj.adjoint
                        x = permute(x, [obj.NDim+1:obj.NDim+length(obj.kSize) 1:obj.NDim]);
                        y = x.*exp(-1i*obj.kbase);
                        y = sum(y,3);
                        y = sum(y,4);
                    else
                        y = x.*exp(1i*obj.kbase);
                        if ~isempty(obj.ksphaDiv)
                            y = 1i*y;
                            phiDiv_a = prepForDirect(obj,obj.ksphaDiv,obj.kconcDiv,1);
                            y = phiDiv_a.*y;
                        end
                        y = sum(y,1);
                        y = sum(y,2);
                        y = reshape(y, obj.kSize);
                    end
                end
                
                % Normalization so that a cartesian Fourier tranform would have
                % the adjoint equal to the inverse
                sz = size(obj.phs_grid.x);
                y = y/sqrt(prod(sz(1:obj.NDim)));
			end
		end

        function y = noPrecompWorker(obj,x)
            y = 0;
            Npnts = ceil(numel(x)/obj.noPrecompNdiv);
            for nc = 1:obj.noPrecompNdiv
                subinds1 = 1 + (nc-1)*Npnts;
                subinds2 = min(Npnts + (nc-1)*Npnts, numel(x));
                if subinds1 <= numel(x)
                    x_sub = x(subinds1:subinds2);
                    if obj.adjoint
                        phs = prepForDirect(obj,obj.phs_spha,obj.phs_conc,...
                            obj.sampTimes,[],[],[subinds1,subinds2]);
                        x_sub = reshape(x_sub, [ones(1,length(obj.imSize)), numel(x_sub)]);
                        tmp = x_sub .* exp(-1i*phs);
                        y = y + sum(tmp,length(obj.imSize)+1);
                    else
                        phs = prepForDirect(obj,obj.phs_spha,obj.phs_conc,...
                            obj.sampTimes,[],[subinds1,subinds2],[]);
                        phs = reshape(phs, numel(x_sub), []);
                        tmp = x_sub(:) .* exp(1i*phs);
                        if ~isempty(obj.ksphaDiv)
                            phiDiv_a = prepForDirect(obj,obj.ksphaDiv,obj.kconcDiv,...
                                1,[],[subinds1,subinds2],[]);
                            phiDiv_a = reshape(phiDiv_a, numel(x_sub), []);
                            tmp = 1i*tmp;
                            tmp = phiDiv_a.*tmp;
                        end
                        y = y + reshape(sum(tmp,1), obj.kSize);
                    end
                end
            end
        end

		function kbase = prepForDirect(obj,phs_spha_a,phs_conc_a,sampTimes_a,sphaInds,subIndsSpace,subIndsTime)
            if nargin<5 || isempty(sphaInds)
                sphaInds = 1:size(obj.phs_spha,1);
            end
            if nargin<6 || isempty(subIndsSpace)
                subIndsSpace = [];
            end
            if nargin<7 || isempty(subIndsTime)
                subIndsTime = [];
            end
            if ~isempty(subIndsTime)
                phs_spha_a = phs_spha_a(:,subIndsTime(1):subIndsTime(2));
                phs_conc_a = phs_conc_a(:,subIndsTime(1):subIndsTime(2));
                if numel(sampTimes_a)>1
                    sampTimes_a = sampTimes_a(subIndsTime(1):subIndsTime(2));
                end
            end
            if ~isempty(subIndsSpace) 
                if ~isempty(obj.b0mask)
                    b0mask_a = obj.b0mask(subIndsSpace(1):subIndsSpace(2));
                end
                b0_a = obj.b0(subIndsSpace(1):subIndsSpace(2));
            else
                if ~isempty(obj.b0mask)
                    b0mask_a = obj.b0mask;
                end
                b0_a = obj.b0;
            end
			% Move all time dims to enable implicit replication when multiplying spatial dims by time dims
            if numel(sampTimes_a) > 1
			    sampTimes_a = permute(sampTimes_a, [length(obj.kSize)+1:length(obj.kSize)+obj.NDim 1:length(obj.kSize)]);
            end
			phs_spha_a = permute(phs_spha_a, [1 length(obj.kSize)+2:length(obj.kSize)+obj.NDim 2:obj.NDim+1]);
			phs_conc_a = permute(phs_conc_a, [1 length(obj.kSize)+2:length(obj.kSize)+obj.NDim 2:obj.NDim+1]);
			% Precompute spatial variation at all times (high memory demand, but very fast)
            phs = 0;
            for n=sphaInds
				% Add all spatially varying spherical harmonic terms
				bfunc = obj.basisFuncSphHarm(n,subIndsSpace);
                phs_a = bfunc .* phs_spha_a(n, 1, :, :);
                if n>4 && ~isempty(obj.b0mask)
                    phs_a = phs_a.*b0mask_a;
                end
				phs = phs + phs_a;
            end
            for n=1:size(phs_conc_a,1)
				% Add all concomitant gradient terms
				bfunc = obj.basisFuncConcGrad(n,subIndsSpace);
                phs_a = bfunc .* phs_conc_a(n, 1, :, :);
                if ~isempty(obj.b0mask)
                    phs_a = phs_a.*b0mask_a;
                end
				phs = phs + phs_a;
            end
            if numel(sampTimes_a)>1
                phs_a = b0_a.*sampTimes_a;
                if ~isempty(obj.b0mask)
                    phs_a = phs_a.*b0mask_a;
                end
                phs = phs + phs_a;
            end
			kbase = phs;
		end

        function obj = setPhiDiv(obj,ksphaDiv,kconcDiv)
            % Check memory again, since PhiDiv greatly increases memory
            % demand.
            G = gpuDevice;
            avMem = G.AvailableMemory;
            numbytes = numel(obj.b0)*numel(obj.sampTimes)*8;
            if numbytes < 0.15*avMem
                obj.useGPU = 1;
                obj.useSingle = 0;
                obj.useInterp = 0;
            elseif numbytes < 0.3*avMem
                obj.useGPU = 1;
                obj.useSingle = 1;
                obj.useInterp = 0;
            else
                obj.useGPU = 1;
                obj.useSingle = 1;
                obj.useInterp = 0;
                obj.noPrecomp = 1;
                obj.noPrecompNdiv = ceil(numbytes / (0.1*avMem));
                obj.kbase = [];
            end
            % ksphaDiv and kconcDiv are temporal derivatives of obj.phs_spha and obj.phs_conc
            if obj.useSingle
                ksphaDiv = single(ksphaDiv);
				kconcDiv = single(kconcDiv);
            end
            % Move variables to GPU
			if obj.useGPU
				ksphaDiv = gpuArray(ksphaDiv);
				kconcDiv = gpuArray(kconcDiv);
            end
            %if obj.noPrecomp
                obj.ksphaDiv = ksphaDiv;
                obj.kconcDiv = kconcDiv;
            %else
            %    obj.phiDiv = prepForDirect(obj,ksphaDiv,kconcDiv,1);
            %end
        end
		
		function [svdTime,svdSpace,traj,phsShft] = prepForInterp(obj)
			% Create nufft object
			% TODO: below assumes perfectly axial slices. To do this properly, need to:
			% 1. have normal vector to slice as an optional input
			% 2. find linear combination of terms 2:4 in kspha for in-plane to slice
			% 3. set that to kloc, and substract from kspha. Then can still have all kspha terms in sum below
            kloc = zeros(2,size(obj.phs_spha,2), 'like', obj.phs_spha);
            if abs(obj.phs_grid.z(2,2)-obj.phs_grid.z(1,1)) ~= 0
                error('non-axial slices not supported for interpolated approach')
            end
            if abs(obj.phs_grid.x(1,2)-obj.phs_grid.x(1,1)) > eps
                if abs(obj.phs_grid.y(1,2)-obj.phs_grid.y(1,1)) ~= 0
                    error('non-axial slices not supported for interpolated approach')
                end
                res1 = abs(obj.phs_grid.y(2,1)-obj.phs_grid.y(1,1));
                res2 = abs(obj.phs_grid.x(1,2)-obj.phs_grid.x(1,1));
                if obj.phs_grid.y(2,1)-obj.phs_grid.y(1,1) > 0
                    kloc(1,:) = obj.phs_spha(3,:)/2/pi*res1;
                else
                    kloc(1,:) = -obj.phs_spha(3,:)/2/pi*res1;
                end
                if obj.phs_grid.x(1,2)-obj.phs_grid.x(1,1) > 0
                    kloc(2,:) = obj.phs_spha(2,:)/2/pi*res2;
                else
                    kloc(2,:) = -obj.phs_spha(2,:)/2/pi*res2;
                end
                % Create phase ramp to center object domain properly, since
                % nufft assumes isocenter is at matrix center
                phsShft = obj.phs_grid.y(end/2+1,1)*obj.phs_spha(3,:) + ...
                    obj.phs_grid.x(1,end/2+1)*obj.phs_spha(2,:);
            else
                res1 = abs(obj.phs_grid.x(2,1)-obj.phs_grid.x(1,1));
                res2 = abs(obj.phs_grid.y(1,2)-obj.phs_grid.y(1,1));
                if obj.phs_grid.x(2,1)-obj.phs_grid.x(1,1) > 0
                    kloc(1,:) = obj.phs_spha(2,:)/2/pi*res1;
                else
                    kloc(1,:) = -obj.phs_spha(2,:)/2/pi*res1;
                end
                if obj.phs_grid.y(1,2)-obj.phs_grid.y(1,1) > 0
                    kloc(2,:) = obj.phs_spha(3,:)/2/pi*res2;
                else
                    kloc(2,:) = -obj.phs_spha(3,:)/2/pi*res2;
                end
                % Create phase ramp to center object domain properly, since
                % nufft assumes isocenter is at matrix center
                phsShft = obj.phs_grid.x(end/2+1,1)*obj.phs_spha(2,:) + ...
                    obj.phs_grid.y(1,end/2+1)*obj.phs_spha(3,:);
            end
            phsShft = reshape(exp(1i*phsShft), size(obj.sampTimes));
			traj = nufftOp(size(obj.b0), kloc',[],obj.useGPU);
			clear kloc
			% Determine full non-linear encoding matrix
			b = prepForDirect(obj,obj.phs_spha,obj.phs_conc,obj.sampTimes,[1,4:size(obj.phs_spha,1)]);
			b = reshape(b, numel(obj.b0), numel(obj.sampTimes));
            % Sub-sample b along time dimension to speed up svd, since
            % phase is slowly varying in time. We do not subsample in
            % space, since high resolution is important for the B0 map
            if obj.subFact>1
                inds = 1:obj.subFact:size(b,2);
                if inds(end) ~= size(b,2)
                    inds = [inds, size(b,2)]';
                end
                b = b(:,inds);
            end
            b = exp(1i*b);
            % Find largest singular values and vectors
            S = 1;
            ntry = 0;
            delTry = 30;
            subspcFact = 3;
            while (min(diag(S))/max(S(:)) > obj.svdThresh) 
                if ntry > 200
                    error('many large singular values. Try providing mask or adjusting svdThresh.')
                end
                ntry = ntry + delTry;
                [U,S,V,FLAG] = svds(b,ntry,'largest','SubspaceDimension',subspcFact*ntry);
                if FLAG
                    warning('svd failure to converge. Increasing subspace.')
                    ntry = ntry - delTry;
                    subspcFact = subspcFact+1;
                end
            end
            Ns = find(diag(S)/max(S(:))<obj.svdThresh,1,'first');
            svdTime = conj(V(:,1:Ns)*S(1:Ns,1:Ns));
            if obj.subFact
                svdTime = interp1(inds,gather(svdTime),1:inds(end),'pchip');
            end
            svdSpace = U(:,1:Ns);
            clear U V
			% Compute error wrt direct approach
			if (0)
				erVal = gather(svdSpace)*gather(svdTime.');
				b = prepForDirect(obj,obj.phs_spha,obj.phs_conc,obj.sampTimes,[1,4:size(obj.phs_spha,1)]);
				b = reshape(b, numel(obj.b0), numel(obj.sampTimes));
                b = exp(1i*b);
                erVal = erVal(:) - gather(b(:));
                erVal = norm(erVal)/norm(b(:))
			end
        end

        function  res = ctranspose(obj)
            obj.adjoint = xor(obj.adjoint,1);
            res = obj;
        end
        
	end
end


