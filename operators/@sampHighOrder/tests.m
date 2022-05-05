function tests(obj)

N0 = 128;

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
S = sampHighOrder(b0,sampTimes,phs_spha,phs_coco,phs_grid,[],[],[],0);
x = randn(imN) + randn;
y = randn(imk) + randn;
Sx = S*x;
Sy = S'*y;
d1 = dot(x(:),Sy(:));
d2 = dot(Sx(:),y(:));
assert(abs(d1-d2)/min(abs(d1),abs(d2)) < 1e-8, 'Adjoint test failed.')

%% Test interpolation. Use random polynomials to simulate slow variations
kloc = projection(N0,N0-4,2);
sampTimes = (1:size(kloc,1))';
[x,y] = meshgrid(-N0/2:N0/2-1,-N0/2:N0/2-1);
dx = 0.8;
dy = 0.95;
kloc(:,1) = 2*pi*kloc(:,1)/max(abs(kloc(:,1)))/2/dx;
kloc(:,2) = 2*pi*kloc(:,2)/max(abs(kloc(:,2)))/2/dy;
x = (x+1)*dx;
y = (y+1)*dy;
scl = 0.5*max(abs(x(:)));
b0 = scl^5*(randn*x + randn*y) + ...
    scl*(randn*x.^4 + randn*(x.^3).*y + randn*(x.^2).*(y.^2) + randn*x.*(y.^3) + randn*y.^4) +...
    randn*x.^5 + randn*(x.^4).*y + randn*(x.^3).*(y.^2) + randn*(x.^2).*(y.^3) + randn*x.*y.^4 + randn*y.^5;
b0 = b0/max(abs(b0(:)))*pi/sampTimes(end); % Ensure theres a decent amount of phase accrual    
sclt = 0.5*sampTimes(end);
phs_spha = zeros(16,length(sampTimes));
phs_coco = zeros(4,length(sampTimes));
phs_spha(2:3,:) = kloc';
eddyAmp = 0.05;
for n= [1,4:16] %1:16
    phs_spha_a = randn + 2*randn/sclt*sampTimes + randn/sclt^2*sampTimes.^2 + randn/sclt^3*sampTimes.^3 + randn/sclt^4*sampTimes.^4;
    if n>1 && n<5
        phs_spha_a = eddyAmp*phs_spha_a/scl;
    elseif n>1 && n<10
        phs_spha_a = eddyAmp*phs_spha_a/scl^2;
    elseif n>=10
        phs_spha_a = eddyAmp*phs_spha_a/scl^3;
    end
    phs_spha(n,:) = phs_spha(n,:) + phs_spha_a';
end
for n=1:4
    phs_coco(n,:) = randn + 2*randn/sclt*sampTimes + randn/sclt^2*sampTimes.^2 + randn/sclt^3*sampTimes.^3 + randn/sclt^4*sampTimes.^4;
    phs_coco(n,:) = eddyAmp*phs_coco(n,:)/scl^2;
end
for orient=1:4
    switch orient
        case 1
            phs_grid.x = x;
            phs_grid.y = y;
        case 2
            phs_grid.x = x';
            phs_grid.y = y';
        case 3
            phs_grid.x = -x;
            phs_grid.y = -y;
        case 4
            phs_grid.x = -x';
            phs_grid.y = -y';
    end
    phs_grid.z = scl*ones(size(x));
    xvec = phantom(N0);
    yvec = fftnc(xvec);
    yvec = yvec(1:numel(sampTimes));
    yvec = yvec(:);
    S = sampHighOrder(b0,sampTimes,phs_spha,phs_coco,phs_grid,[],[],[],0);
    Sx = S*xvec;
    Sy = S'*yvec(:);
    %     tic1 = tic;
    %     for n=1:10
    %         Sx = S*x;
    %     end
    %     toc1 = toc(tic1)
    S = sampHighOrder(b0,sampTimes,phs_spha,phs_coco,phs_grid,[],[],[],1,0.01);
    Sx_b = S*xvec;
    Sy_b = S'*yvec(:);
    d1 = dot(xvec(:),Sy_b(:));
    d2 = dot(Sx_b(:),yvec(:));
    assert(abs(d1-d2)/min(abs(d1),abs(d2)) < 1e-8, 'Interp adjoint test failed.')
    cost = Sx(:)-Sx_b(:);
    cost = sqrt(sum(cost.*conj(cost)))/numel(cost);
    assert(cost < 1e-4, 'Interp comparison to direct failed.')
    %     tic2 = tic;
    %     for n=1:10
    %         Sx_b = S*x;
    %     end
    %     toc2 = toc(tic2)
    %     figure; 
    %     subplot(2,3,1); imagesc(abs(reshape(Sx,N0,N0))); colorbar; subplot(2,3,2); imagesc(abs(reshape(Sx_b,N0,N0))); colorbar; 
    %     subplot(2,3,3); imagesc(angle(reshape(Sx,N0,N0))); subplot(2,3,4); imagesc(angle(reshape(Sx_b,N0,N0)));
    %     subplot(2,3,5); imagesc(abs(reshape(Sx,N0,N0)-reshape(Sx_b,N0,N0))); colorbar;
end

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
S = sampHighOrder(b0,sampTimes,phs_spha,phs_coco,phs_grid,[],[],[],0);
x = randn([imN, nExtra]) + randn;
y = randn([imk, nExtra]) + randn;
Sx = S*x;
Sy = S'*y;
d1 = dot(x(:),Sy(:));
d2 = dot(Sx(:),y(:));
assert(abs(d1-d2)/min(abs(d1),abs(d2)) < 1e-8, 'Adjoint test failed.')

%% Test single precision option
S = sampHighOrder(b0,sampTimes,phs_spha,phs_coco,phs_grid,[],1,1,0);
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
S = sampHighOrder(b0,sampTimes,phs_spha,phs_coco,phs_grid,[],[],[],0);
k0 = S*im0;
k1 = fftnc(im0,2);
cost = abs(k0)-abs(k1); cost = sqrt(sum(cost(:).*conj(cost(:))))/numel(cost);
assert(cost<1e-10, 'fft comparison failed');
im1 = S'*k0;
cost = im1-im0; cost = sqrt(sum(cost(:).*conj(cost(:))))/numel(cost);
assert(cost<1e-10, 'Inverse test failed');


fprintf('sampHighOrder unit test success!\n')


end
