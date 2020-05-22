classdef sampHighOrder
% Perform MRI sampling using direct summation of complex exponentials.
% Supports B0 map, time-varying spherical harmonic distributions of phase,
% and time-varying distributions of phase that are typical for concomitant
% gradients. Currently 2D only.
%
% Usage: 
%   Define sampHighOrder object using 
%       S = sampHighOrder(b0,sampTimes,phs_spha,phs_conc,phs_grid,useGPU,useSingle)
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
%   useGPU: whether to use GPU. Default = true.
%   useSingle: if true, use single instead of double precision. Default =
%       false. If set to true, memory savings will only be realized if the
%       array that S operates on is also single.
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
	end

	methods

		function obj = sampHighOrder(b0,sampTimes,phs_spha,phs_conc,phs_grid,useGPU,useSingle)
			if nargin == 0
				obj.tests;
				return;
			end
			if nargin>6
				obj.useGPU = useGPU;
            end
            if nargin>6
				obj.useSingle = useSingle;
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
			% Use single precision if requested
			if obj.useSingle
                obj.sampTimes = single(obj.sampTimes);
				obj.phs_spha = single(obj.phs_spha);
				obj.phs_conc = single(obj.phs_conc);
				obj.b0 = single(obj.b0);
				obj.phs_grid.x = single(obj.phs_grid.x);
				obj.phs_grid.y = single(obj.phs_grid.y);
				obj.phs_grid.z = single(obj.phs_grid.z);
			end
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
			szx = [size(x), 1, 1];
            if obj.useSingle
                if (isa(x,'gpuArray') && ~isaUnderlying(x,'single')) || ~isa(x,'single')
                    warning('useSingle specified, but input is not single.')
                end
            end

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
        end
        
        function  res = ctranspose(obj)
            obj.adjoint = xor(obj.adjoint,1);
            res = obj;
        end
        
	end
end


