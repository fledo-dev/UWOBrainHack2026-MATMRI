function traj = computeTraj(obj,probe_positions,probe_raw,gammaProbes,gammaMRI,dt,B0,fieldOffsets,coilParams,nonLinSphHarm,fitOrder)
% Compute trajectory for high-order encoding model using raw data from
% field probes.

% Specify whether to weight probe contributions to fit based on their
% distance to isocenter
doWeights = 1;

% Specify whether to build up the trajectory spherical harmonics in steps
doSteps = 1;

% Determine number of probes
Nprobe = size(probe_positions,1);

% Set defaults
if isempty(fitOrder)
    if Nprobe >= 16
        fitOrder = 3;
    elseif Nprobe >= 9
        fitOrder = 2;
    else
        fitOrder = 1;
    end
end

% Correct for probe field offsets
times = (1:size(probe_raw,1))';
phsCor = 2*pi*(gammaProbes*dt)*fieldOffsets(:)'.*times;
probe_raw = probe_raw .* exp(-1i*phsCor);

% Prep. 
phsRaw = unwrap(angle(probe_raw));
magRaw = abs(probe_raw); 

% Weight probes based on lifetimes based on magnitude data. 
%   The magnitude profiles are normalized based on the initial point, and
%   then only the final 100 points are used to create the weight. Thus,
%   those with shorter lifetimes will have smaller weights.
W_mag = magRaw./magRaw(1,:,:,:,:);
W_mag = mean(mean(W_mag(end-100:end,:,:),1),3);
W_mag = W_mag(:)/max(W_mag);

% Have the phase start at 0 for all probes
phsRaw = phsRaw - phsRaw(1,:,:,:);

% Set probe positions
obj.phs_grid.x = probe_positions(:,1);
obj.phs_grid.y = probe_positions(:,2);
obj.phs_grid.z = probe_positions(:,3);

% Compute squared distance from isocenter, to be used for weighting
% Using net distance from isocenter
W = sum(probe_positions.^2,2);
    % Other metrics to potentially try instead:
    % Using sum of all distances from isocenter
    %W = sum(abs(probe_positions),2).^2;
    % Using max z-distance from isocenter
    %W = probe_positions(:,3).^2;

W = 1./W; % Closer probes should have larger weights
%W = 0.01^2 + max(W) - W; % Closer probes should have larger weights. scale by 1 cm. TODO: this approach not sufficiently tested

% Apply weights based on magnitude data 
%W = W.*W_mag; % This worsened quality for DWI

% Normalize weights
W = W/max(W);

% Apply weights to phase input
if doWeights
    phsRaw = W(:)'.*phsRaw;
end

% Create matrix for maxwell terms basis functions
Nc = 4;
B = zeros(Nprobe, Nc);
for l=1:Nc
    B(:,l) = obj.basisFuncConcGrad(l);
end
if doWeights
    B = diag(W)*B;
end

% Set number of spherical harmonics
if fitOrder == 3
    Nl = 16;
elseif fitOrder == 2
    Nl = 9;
else
    Nl = 4;
end

% Iterate between computing SphHarm and phase from conc grads
Nit = 3;
phs_spha = zeros(size(phsRaw,1),Nl,size(phsRaw,3),size(phsRaw,4),'like',phsRaw);
phs_conc = zeros(size(phsRaw,1),Nc,size(phsRaw,3),size(phsRaw,4),'like',phsRaw);
for nv = 1:size(phsRaw(:,:,:),3)
    phsRaw_conc = 0;
    phsRaw_a = phsRaw(:,:,nv);
    if ~doSteps
        % Fit all the orders of spherical harmonics simultaneously
        norderAll = fitOrder;
    else
        % Incrementally increase order, because solving for all terms
        % simultanously can cause error in the lower order terms.
        % NB: shouldn't start with order 0, because if the probes are net
        % offcenter it will mess up the B0 estimation.
        norderAll = 1:fitOrder;
    end
    % Initialize phase that is subtracted from raw phase after fitting
    % lower orders 
    phs_pre = 0; 
    % Loop through orders (required for "doSteps" option)
    for norder_ind = 1:length(norderAll)
        norder = norderAll(norder_ind);
        % Set matrix equation for normal spherical harmonic expansion of trajectory
        if norder == 3
            Nla = 10;
            Nl = 16;
        elseif norder == 2
            Nla = 5;
            Nl = 9;
        elseif norder == 1
            Nla = 1;
            Nl = 4;
        else
            error('Must choose order between 1 and 3')
        end
        if ~doSteps
            % Here we need to do all orders at once, so we start from the
            % first one (i.e., Nla=1) and go up to the one specified via
            % norder above (i.e., Nl)
            Nla = 1;
        end
        A = zeros(Nprobe, Nl-Nla+1);
        for l=Nla:Nl
            A(:,l-Nla+1) = obj.basisFuncSphHarm(l);
        end
        if doWeights
            A = diag(W)*A;
        end
        pA = pinv(A);
        for n=1:Nit
            phs_spha_a = pA*(phsRaw_a.' - phsRaw_conc - phs_pre);
            if (norder == 1) || ((Nla==1) && (Nl>=4))
                linInds = (2:4) - Nla + 1;
                phs_conc_a = obj.computeMaxwellPhase(phs_spha_a(linInds,:).',dt,gammaProbes,coilParams,B0,nonLinSphHarm);
                phsRaw_conc = B*(phs_conc_a.');
                resid(n) = norm(reshape(A*phs_spha_a + B*(phs_conc_a.'),[],1) - phsRaw_a(:));
            end
        end
        phs_pre = phs_pre + A*(phs_spha_a);
        % Save result into output
        phs_spha(:,Nla:Nl,nv) = phs_spha_a.';
    end
    phs_conc(:,:,nv) = phs_conc_a;
end
%figure; plot(resid)


% Format output and scale phase by gammas
traj.phs_spha = phs_spha * gammaMRI / gammaProbes;
traj.phs_conc = phs_conc * gammaMRI / gammaProbes;


end


