classdef fftOp < handle
	properties
		trueAdj = 1;
		doAtA = 0;
		imgN = [];
		kfull = []; % kspace size after fftOp*image. Allows for vector output, which lsqr expects.
		imgNfull = [];
		phaseMap = []; % Applies the phase map before fft in forward operator
		realAfterConj = 0;
	end

	properties (SetAccess = protected)
		adjoint = 0;
		sampMask = [];
		shiftedMask = [];
		ND = [];
	end

	methods

		function obj = fftOp(sampMask,ND)
			if nargin == 0
				sz = 64*[1 1 1];
				sampMask = rand(sz) > 0.5;
				obj = fftOp(sampMask);
				res = testFftOp(obj);
				return;
			end
			if nargin < 2
				obj.ND = ndims(sampMask);
			else
				obj.ND = ND;
			end
			szSamp = size(sampMask);
			shiftedMask = sampMask;
			for n=1:obj.ND
				shiftedMask = fftshift(shiftedMask,n);
			end
			obj.shiftedMask = shiftedMask>0;
			obj.sampMask = sampMask>0;
			obj.imgN = size(sampMask);
		end

		function y = mtimes(obj,x)
			if obj.doAtA
				y = obj.fcnAtA_a(x);
				return;
			end

			if isa(x,'gpuArray') && ~isa(obj.sampMask,'gpuArray')
				obj.sampMask = gpuArray(obj.sampMask);
			end

			sz = size(x); sz = [sz ones(1,10)];
			ND = obj.ND;

			if obj.adjoint
				if isempty(obj.kfull)
					y = x;
				else
					y = reshape(x,obj.kfull);
					sz = obj.kfull;
				end
				y = y.*repmat(obj.sampMask, [ones(1,ndims(obj.sampMask)) sz(ndims(obj.sampMask)+1:end)]);
				for n=1:ND
					y = ifftshift(ifft(ifftshift(y,n),[],n),n);
				end
				if obj.trueAdj
					y = y*sqrt(prod(sz(1:ND)));
				end
				if ~isempty(obj.phaseMap)
					y = y.*conj(obj.phaseMap);
				end
				if obj.realAfterConj
					y = real(y);
				end
				if ~isempty(obj.imgNfull)
					y = y(:);
				end
				obj.adjoint = 0; % Because this is a handle class
			else
				if isempty(obj.imgNfull)
					y = x;
				else
					y = reshape(x,obj.imgNfull);
					sz = obj.imgNfull;
				end
				if ~isempty(obj.phaseMap)
					y = y.*obj.phaseMap;
				end
				for n=1:ND
					y = fftshift(fft(fftshift(y,n),[],n),n);
				end
				if obj.trueAdj
					y = y/sqrt(prod(sz(1:ND)));
				end
				y = y.*repmat(obj.sampMask, [ones(1,ndims(obj.sampMask)) sz(ndims(obj.sampMask)+1:end)]);
				if ~isempty(obj.kfull)
					y = y(:);
				end
			end
		end

		function y = fcnAtA_a(obj, x)
			if isa(x,'gpuArray') && ~isa(obj.shiftedMask,'gpuArray')
				shiftedMask = gpuArray(obj.shiftedMask);
			else
				shiftedMask = obj.shiftedMask;
			end
			sz = size(x); sz = [sz ones(1,10)];
			y = x;
			for n=1:obj.ND
				y = fft(y,[],n);
			end
			y = y.*repmat(shiftedMask, [ones(1,ndims(shiftedMask)) sz(ndims(shiftedMask)+1:end)]);
			for n=1:obj.ND
				y = ifft(y,[],n);
			end
		end

		function res = testFftOp(obj)
			fprintf('Testing fftOp...')
			% Test that adjoint is truly adoint.
			sz = size(obj.sampMask);
			a = randn(sz) + 1i*randn(sz);
			b = randn(sz) + 1i*randn(sz);
			Aa = obj*a;
			ATb = obj'*b;
			testVal = dot(a(:),ATb(:)) - dot(Aa(:),b(:));
			if abs(testVal) > 1e-8
				error('FAILURE: fftOp adjoint test.')
			end
			% Test that AtA fcn gives expected result
			a = phantom3dm(sz);
			ATA1 = obj' * (obj*a);
			obj.doAtA = 1;
			ATA2 = obj*a;
			testVal = ATA1(:) - ATA2(:);
			testVal = sum(testVal.*conj(testVal));
			if abs(testVal) > 1e-8
				error('FAILURE: fftOp AtA test.')
			end
			res = 1;
			fprintf('success\n')
		end

	end
end
