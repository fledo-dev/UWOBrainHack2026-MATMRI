function traj = computeTraj(obj,probe_positions,probe_raw,gammaProbes,gammaMRI,dt,B0,fieldOffsets,coilParams,nonLinSphHarm,fitOrder)
% Compute trajectory for high-order encoding model using raw data from
% field probes.

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
magRaw = abs(probe_raw); % TODO: could use this for weighted least squares
if fitOrder == 3
    Nl = 16;
elseif fitOrder == 2
    Nl = 9;
else
    Nl = 4;
end

% Have the phase start at 0 for all probes
phsRaw = phsRaw - phsRaw(1,:,:,:);

% Set matrix equation based on basis functions
obj.phs_grid.x = probe_positions(:,1);
obj.phs_grid.y = probe_positions(:,2);
obj.phs_grid.z = probe_positions(:,3);
A = zeros(Nprobe, Nl);
for l=1:Nl
    A(:,l) = obj.basisFuncSphHarm(l);
end
pA = pinv(A);

% Create similar matrix, but for phase from maxwell terms
Nc = 4;
B = zeros(Nprobe, Nc);
for l=1:Nc
    B(:,l) = obj.basisFuncConcGrad(l);
end

% Iterate between computing SphHarm and phase from conc grads
Nit = 3;
phsRaw_conc = 0;
phs_spha = zeros(size(phsRaw,1),Nl,size(phsRaw,3),size(phsRaw,4),'like',phsRaw);
phs_conc = zeros(size(phsRaw,1),Nc,size(phsRaw,3),size(phsRaw,4),'like',phsRaw);
for nv = 1:size(phsRaw(:,:,:),3)
    phsRaw_a = phsRaw(:,:,nv);
    for n=1:Nit
        phs_spha_a = pA*(phsRaw_a.'-phsRaw_conc);
        phs_conc_a = obj.computeMaxwellPhase(phs_spha_a(2:4,:).',dt,gammaProbes,coilParams,B0,nonLinSphHarm);

        phsRaw_conc = B*(phs_conc_a.');
        
        resid(n) = norm(reshape(A*phs_spha_a + B*(phs_conc_a.'),[],1) - phsRaw_a(:));
    end
    phs_spha(:,:,nv) = phs_spha_a.';
    phs_conc(:,:,nv) = phs_conc_a;
end
%figure; plot(resid)


% Format output and scale phase by gammas
traj.phs_spha = phs_spha * gammaMRI / gammaProbes;
traj.phs_conc = phs_conc * gammaMRI / gammaProbes;


end


