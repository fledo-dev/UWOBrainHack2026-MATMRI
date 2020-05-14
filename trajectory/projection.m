function [kloc] = projection(Nproj,NperProj,ndim)
%PROJECTION Summary of this function goes here
%   Detailed explanation goes here
% (c) Corey Baron 2015
%

switch ndim
    case 2
        phi = linspace(0,pi,Nproj+1);
        phi = repmat(phi(1:end-1),[NperProj, 1]);
        r = linspace(-0.5,0.5,NperProj+1);
        r = repmat(r(1:end-1)', [1 Nproj]);
        kloc = cat(2,r(:).*cos(phi(:)), r(:).*sin(phi(:))); 
    otherwise 
        error('requested dim not coded!')
end

end

