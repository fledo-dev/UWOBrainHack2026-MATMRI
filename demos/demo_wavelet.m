%% sampHighOrder with wavelet regularization demo
%
%  (c) Corey Baron, 2020

% To make things simple, we'll start with actual acquired data for 
% the field probe trajectories, B0 map, receiver sensitivity, and a ground
% truth image
load dataForDemonstrations.mat

% Build sampling object
grid.x = X; 
grid.y = Y;
grid.z = Z;
S = sampHighOrder(b0map,datatime,phs_spha',phs_conc',grid);

% Build receiver operator
R = rcvrOp(Crcvr,0);

% Sample image
data = S*(R*im0);

% Add some noise
ns_std = 5;
data = data + ns_std*(randn(size(data)) + 1i*randn(size(data)));


%% Normal highorder recon using conjugate gradient method
% Note that by default, cgne stops on the first iteration 
opFunc = @(x,transp) mrSampFunc(x,transp,S,R);
maxIt = 30;
opt.stopOnResInc = 1; % Stop on the first iteration where the residual increases, which turns out to be a pretty good criteria 
[im, resvec, mse, xnorm, xdiff] = cgne(opFunc,data,[],maxIt); 


%% Wavelet regularization
% Solve argmin(||Ax-b||^2_2 + lambda*||Wx||_1)
% Build discrete wavelet transform object
levelsPerDim = [2 2];
isDec = 0;
W = dwt(levelsPerDim,size(im),isDec);
wim = W*im;
figure; imagesc(wim); 

% Standard compressed sensing. Just using a lambda that looks okay "by eye"
NitMax = 100;
opt_bfista.gtruth = im0;
x0 = [];
lambda = 1; 
[imCS, resSqAll, RxAll, mseAll] = bfista(opFunc,data,W,lambda,x0,NitMax,opt_bfista);
close(1);
figure(1); 
subplot(1,3,1); plot(log10(resSqAll)); hold('all')
subplot(1,3,2); plot(log10(RxAll)); hold('all')
subplot(1,3,3); plot(log10(mseAll)); hold('all')

% Automatic lambda selection, similar to Varela-Mattatall G, Baron CA,
% Menon RS. Magn. Reson. Med. 2021;86:1403–1419. Here we'll base the lambda
% selection per level on histograms of wavelet xform of "zero filled"
% image, which is more generally just the adjoint of the encoding matrix A.
% Here we just use the location where the histogram of the magnitude image
% is maximized, which corresponds to the std of the noise if we assume that
% the the distribution is dominated by the {noise + aliasing} that we want
% to suppresss with wavelet regularization (i.e., it's roughly a Rayleigh
% distribution). We divide by sqrt(2) b/c we're dealing with magnitude
% while the bfista algorithm is using complex in softthresholding.
% This assumption is most valid for the 1st level high freq transform, so
% we use this lambda for all the other levels (and set the lambda for the
% low frequencies to 0, since it is not sparse).
im_zf = opFunc(data,'transp');
wim_zf = W*im_zf;
figure(2); 
subplot(1,2,1);
imagesc(wim_zf);
subplot(1,2,2)
h = histogram(gather(abs(wim_zf.high{n}(:))),round(numel(wim_zf.high{n}(:))/100));
xlim([0 20]);
vals = h.Values;
vals(1) = 0; % The first bin can be inflated if there was masking in im0
[~,threshInd] = max(vals);
lamVal = 0.5*(h.BinEdges(threshInd) + h.BinEdges(threshInd+1))/sqrt(2);
lambda2 = waveletObj(lamVal,wim_zf); % This replicates the scalar lamVal into all wavelet fields using the template wim_zf
lambda2.low = 0;

% Run the recon
[imCS2, resSqAll2, RxAll2, mseAll2] = bfista(opFunc,data,W,lambda2,x0,NitMax,opt_bfista);
figure(1);
subplot(1,3,1); plot(log10(resSqAll2)); ylabel('Data consistency'); xlabel('iterations');
subplot(1,3,2); plot(log10(RxAll2)); ylabel('Wavelet l1 norm'); xlabel('iterations');
subplot(1,3,3); plot(log10(mseAll2)); ylabel('Mean squared error'); xlabel('iterations');
legend('imCS','imCS2')








