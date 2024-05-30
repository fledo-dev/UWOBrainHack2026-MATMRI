function extractLTE(nii_in,bmatFile, bvalFile, bvecFile)
% Function to extract LTE scans from a b-tensor encoding scan nifti file

% Load data
info = niftiinfo(nii_in);
im = niftiread(nii_in);
brank = bmatRank(load(bmatFile));
bval = load(bvalFile);
bvec = load(bvecFile);

% Get filename
fname = nii_in;
while contains(fname,'.')
    [~,fname] = fileparts(fname);
end
fname = [fname, '_LTE'];

% Extract LTE scans
im = im(:,:,:,brank < 1.5);
info.ImageSize(4) = size(im,4);
info.raw.dim(5) = size(im,4);
niftiwrite(im,[fname,'.nii'],info,'Compressed',true);

if nargin > 2
    bval = bval(brank < 1.5);
    fid = fopen([fname, '.bval'], 'w');
    fprintf(fid, '%d ', round(bval(1:end-1)));
    fprintf(fid, '%d', round(bval(end)));
    fclose(fid);
end
if nargin > 3
    bvec = bvec(:, brank < 1.5);
    fid = fopen([fname, '.bvec'], 'w');
    for n=1:3
        fprintf(fid, '%.6f ', bvec(n,1:end-1));
        fprintf(fid, '%.6f\n', bvec(n,end));
    end
    fclose(fid);
end


end