classdef sampHighOrder
%
%  sampDirectMat(NDim,b0,sampTimes,phs_spha,phs_conc,phs_grid)
%    - Object for 2D direct matrix based sampling on a single slice. Fairly high memory
%    overhead due to precomputations, but fast.
%    - input array can by higher than 2D. In this case the multiplication is just looped over all the higher dims.
% TODO: put useful help info here
% TOOD: perhaps the new sampWithRateMap method I developed should actually be implemented in here. 
%   useMeth could be an input option.

	properties
		NDim = 2;   % Number of dims to do sums over in both image domain and k-space. 
        kfull = []; % total kspace size after obj*image. Allows for vector output, which lsqr expects.
		imgNfull = []; % total image size after adjoint(obj)*k. Allows for vector output on adjoint, which lsqr expects.
	end

	properties (SetAccess = protected)
		adjoint = 0;
		b0 = []; % rad/sec
        phs_spha = []; % 0th, 1st, 2nd and 3rd order spherical harmonic terms for phase over time. Must have dimensions 16 x size(sampTimes). 
        phs_conc = []; % conc grad terms for phase over time. Must have dimensions 4 x size(sampTimes). 
        phs_grid = []; % Struct with fields phs_grid.X, phs_grid.Y, phs_grid.Z. Each must have dimensions equivalent to b0. MUST be in magnet frame.
		sampTimes = []; % sec
		kbase = []; % spatial part of spherical harmonics that is multiplied with phs_spha or phs_conc terms
		kSize = [];
		imSize = [];
		useGPU = 1;
	end

	methods

		function obj = sampHighOrder(b0,sampTimes,phs_spha,phs_conc,phs_grid,useGPU)
			if nargin == 0
				obj.tests;
				return;
			end
			if nargin>6
				obj.useGPU = useGPU;
			end
			obj.NDim = 2; % This class coded/optimized for 2D only
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
				error('Only up to 2 dimensions allowed');
			end
			obj.phs_conc = phs_conc;
            obj.phs_grid = phs_grid;
            if obj.NDim < 2
                % This is needed because size() always outputs at least two
                % values, even for vectors. Can hack 1D by having all size
                % of b0 = 1 x N.
                error('need at least 2 dimensions')
            end
			% Move all time dims to enable implicit replication when multiplying spatial dims by time dims
			obj.sampTimes = permute(obj.sampTimes, [length(obj.kSize)+1:length(obj.kSize)+obj.NDim 1:length(obj.kSize)]);
			obj.phs_spha = permute(obj.phs_spha, [1 length(obj.kSize)+2:length(obj.kSize)+obj.NDim 2:obj.NDim+1]);
			obj.phs_conc = permute(obj.phs_conc, [1 length(obj.kSize)+2:length(obj.kSize)+obj.NDim 2:obj.NDim+1]);
			% Move variables to GPU
			if obj.useGPU
				obj.sampTimes = gpuArray(obj.sampTimes);
				obj.phs_spha = gpuArray(obj.phs_spha);
				obj.phs_conc = gpuArray(obj.phs_conc);
				obj.b0 = gpuArray(obj.b0);
				obj.phs_grid.x = gpuArray(obj.phs_grid.x);
				obj.phs_grid.y = gpuArray(obj.phs_grid.y);
				obj.phs_grid.z = gpuArray(obj.phs_grid.z);
			end
			% Precompute spatial variation at all times (high memory demand, but very fast)
            phs = 0;
			for n=1:size(obj.phs_spha,1)
				% Add all spatially varying spherical harmonic terms
				bfunc = obj.basisFuncSphHarm(n);
				phs = phs + bfunc .* obj.phs_spha(n, 1, :, :);
			end
			for n=1:size(obj.phs_conc,1)
				% Add all concomitant gradient terms
				bfunc = obj.basisFuncConcGrad(n);
				phs = phs + bfunc .* obj.phs_conc(n, 1, :, :);
			end
			phs = phs + obj.b0.*obj.sampTimes;
			obj.kbase = phs;
		end

		function y = mtimes(obj,x)
			if ~obj.adjoint && ~isempty(obj.imgNfull)
				x = reshape(x,obj.imgNfull);
			elseif obj.adjoint && ~isempty(obj.kfull)
				x = reshape(x,obj.kfull);
			end
			szx = [size(x), 1, 1];

			if obj.adjoint
				% TODO: might have to make this single precision for large arrays
				y = zeros([obj.imSize, szx(3:end)], 'like', x);
				for n=1:prod(szx(3:end))
					x_a = x(:,:,n);
					x_a = permute(x_a, [obj.NDim+1:obj.NDim+length(obj.kSize) 1:obj.NDim]);
					if obj.useGPU && ~isa(x_a, 'gpuArray')
						x_a = gpuArray(x_a);
					end
					y_a = x_a.*exp(-1i*obj.kbase);
					y_a = sum(y_a,3);
					y_a = sum(y_a,4);
					if ~isa(x,'gpuArray')
						y_a = gather(y_a);
					end
					y(:,:,n) = y_a;
				end
            else
				% TODO: might have to make this single precision for large arrays
				y = zeros([obj.kSize, szx(3:end)], 'like', x);
				for n=1:prod(szx(3:end))
					x_a = x(:,:,n);
					if obj.useGPU && ~isa(x_a, 'gpuArray')
						x_a = gpuArray(x_a);
					end
					y_a = x_a.*exp(1i*obj.kbase);
					y_a = sum(y_a,1);
					y_a = sum(y_a,2);
					if ~isa(x,'gpuArray')
						y_a = gather(y_a);
					end
					y(:,:,n) = reshape(y_a, obj.kSize);
				end
            end
            
            % Normalization so that a cartesian Fourier tranform would have
            % the adjoint equal to the inverse
            sz = size(obj.phs_grid.x);
            y = y/sqrt(prod(sz(1:obj.NDim)));

			if ~obj.adjoint && ~isempty(obj.kfull)
				y = y(:);
			elseif obj.adjoint && ~isempty(obj.imgNfull)
				y = y(:);
			end

        end
        
        function  res = ctranspose(obj)
            obj.adjoint = xor(obj.adjoint,1);
            res = obj;
        end
        
	end
end


