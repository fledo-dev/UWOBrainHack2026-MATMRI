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
% Build discrete wavelet transform object
levelsPerDim = [2 2];
isDec = 0;
W = dwt(levelsPerDim,size(im),isDec);
wim = W*im;
figure; imagesc(wim); 