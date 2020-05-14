% Note: we avoid genpath on top repo dir as it will add the whole .git
% folder and its subdirectories too.

% Delete any parallel pools, because the workers won't get the new path otherwise
if(exist('gcp'))
  poolobj = gcp('nocreate');
  delete(poolobj);
end

topPath = mfilename('fullpath');
topPath = topPath(1:end-length(mfilename));

% Add paths
addpath([topPath, 'bview']);
addpath([topPath, 'demos']);
addpath([topPath, 'general']);
addpath([topPath, 'operators']);
addpath([topPath, 'trajectory']);
