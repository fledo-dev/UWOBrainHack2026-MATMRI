%% Regridding Demo
%
%  Regridding is performed by first defining a type-2 nuFFT object,
%  updating it's density compensation function (dcf), and then applying the
%  adjoint of the nuFFT.
%


%% Generate Simulation k-space data using nuFFT
% Create spiral trajectory


im = phantom;



