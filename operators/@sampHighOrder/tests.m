function tests(obj)

N0 = 64;

%% Check adjoint for basic case
imN = [N0 N0];
b0 = randn(imN) + rand;
imk = [N0 1];
sampTimes = randn(imk)/sqrt(prod(imk));
phs_spha = randn([16, size(sampTimes)]);
phs_coco = randn([4, size(sampTimes)]);
phs_grid.x = randn(size(b0));
phs_grid.y = randn(size(b0));
phs_grid.z = randn(size(b0));
S = sampHighOrder(b0,sampTimes,phs_spha,phs_coco,phs_grid);
x = randn(imN) + randn;
y = randn(imk) + randn;
Sx = S*x;
Sy = S'*y;
d1 = dot(x(:),Sy(:));
d2 = dot(Sx(:),y(:));
assert(abs(d1-d2)/min(abs(d1),abs(d2)) < 1e-8, 'Adjoint test failed.')

%% Check adjoint for larger im dims
imN = [N0 N0];
nExtra = [2 2];
b0 = randn(imN) + rand;
imk = [N0 2];
sampTimes = randn(imk)/sqrt(prod(imk));
phs_spha = randn([16, size(sampTimes)]);
phs_coco = randn([4, size(sampTimes)]);
phs_grid.x = randn(size(b0));
phs_grid.y = randn(size(b0));
phs_grid.z = randn(size(b0));
S = sampHighOrder(b0,sampTimes,phs_spha,phs_coco,phs_grid);
x = randn([imN, nExtra]) + randn;
y = randn([imk, nExtra]) + randn;
Sx = S*x;
Sy = S'*y;
d1 = dot(x(:),Sy(:));
d2 = dot(Sx(:),y(:));
assert(abs(d1-d2)/min(abs(d1),abs(d2)) < 1e-8, 'Adjoint test failed.')

%% Test single precision option
S = sampHighOrder(b0,sampTimes,phs_spha,phs_coco,phs_grid,1,1);
x = single(randn(imN) + randn);
y = single(randn(imk) + randn);
Sx = S*x;
Sy = S'*y;
d1 = dot(x(:),Sy(:));
d2 = dot(Sx(:),y(:));
assert(abs(d1-d2)/min(abs(d1),abs(d2)) < 1e-5, 'Single test failed.')

%% Test vs Fourier transform
N=64;
[phs_grid.x, phs_grid.y] = meshgrid(-N/2:N/2-1,-N/2:N/2-1);
phs_grid.z = zeros(size(phs_grid.x));
phs_spha = zeros(16,N,N);
phs_coco = zeros(4,N,N);
phs_spha(3,:,:) = pi*2/N*reshape(repmat(-N/2:N/2-1, [1 N]), [N,N]);
phs_spha(2,:,:) = pi*2/N*repmat(-N/2:N/2-1, [N 1]);
sampTimes = reshape(0:N^2-1, N, N);
b0 = zeros(size(phs_grid.x));
im0 = phantom(N);
S = sampHighOrder(b0,sampTimes,phs_spha,phs_coco,phs_grid);
k0 = S*im0;
k1 = fftnc(im0,2);
cost = abs(k0)-abs(k1); cost = sqrt(sum(cost(:).*conj(cost(:))))/numel(cost);
assert(cost<1e-10, 'fft comparison failed');
im1 = S'*k0;
cost = im1-im0; cost = sqrt(sum(cost(:).*conj(cost(:))))/numel(cost);
assert(cost<1e-10, 'Inverse test failed');

fprintf('Unit test success!\n')

end
