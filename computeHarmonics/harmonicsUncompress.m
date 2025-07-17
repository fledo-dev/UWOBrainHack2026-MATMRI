function phs_spha = harmonicsUncompress(phs_spha,compMat)
% Uncompress basis functions
% phs_spha [Nsamples, Nbasis, otherDimensions]

if ndims(phs_spha)>6
    error('not coded for more than 6 dims')
end

phs_spha = permute(phs_spha,[2 1 3:6]);
sz_a = size(phs_spha);
compMatInv = pinv(compMat);
% Scale to SI units (compMat uses units of dm)
compMatInv(1:5,:) = 10^2*compMatInv(1:5,:); % 2nd order
if size(compMatInv,1) > 5
    compMatInv(6:12,:) = 10^3*compMatInv(6:12,:); % 3nd order
end
if size(compMatInv,1) > 12
    compMatInv(13:21,:) = 10^4*compMatInv(13:21,:); % 4th order
end
if size(compMatInv,1) > 21
    compMatInv(22:32,:) = 10^5*compMatInv(22:32,:); % 5th order
end
if size(compMatInv,1) > 32
    error('higher than 5th order compression not yet coded!')
end
% Uncompress 2nd and higher orders
phs_spha_a = compMatInv*phs_spha(5:end,:);
phs_spha_a = reshape(phs_spha_a, [size(phs_spha_a,1), sz_a(2:end)]);
phs_spha = cat(1, phs_spha(1:4,:,:,:,:,:), phs_spha_a);
phs_spha = permute(phs_spha,[2 1 3:6]);

end

