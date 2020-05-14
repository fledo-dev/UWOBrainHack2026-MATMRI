function m_sos = calcSOS(m, dim)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if (nargin<2)
  dim = ndims(m);
end

m_sos = squeeze(sqrt(sum(m.*conj(m), dim)));
