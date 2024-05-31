function extractLTE(nii_in,bmatFile, bvalFile, bvecFile, saveMode, bthresh)
% Function to extract LTE scans from a b-tensor encoding scan nifti file

if nargin < 3
    bvalFile = [];
end
if nargin < 4
    bvecFile = [];
end
if (nargin < 5) || isempty(saveMode)
    % 0: save LTE + b0. bthresh is ignored for this case, and b0 is
    %    included based on code in bmatRank.m
    % 1: save LTE + b0 in one file, and STE in another file. bthresh is
    %    ignored for this case, and b0 is included based on code in bmatRank.m
    % 2: save b0, LTE, and STE in three separate files. Useful for doing
    %    fsl eddy separately on b0+LTE and b0+STE. Uses bthresh.
    saveMode = 0;
end
if (nargin < 6) || isempty(bthresh)
    % Threshold for considering scan to be b ~ 0
    bthresh = 20;
end

% Load data
info = niftiinfo(nii_in);
im = niftiread(nii_in);
brank = bmatRank(load(bmatFile));
if ~isempty(bvalFile)
    bval = load(bvalFile);
end
if ~isempty(bvecFile)
    bvec = load(bvecFile);
end

% Get filename
fname = nii_in;
while contains(fname,'.')
    [~,fname] = fileparts(fname);
end

% Extract and save subsets
if saveMode < 2
    inds = brank < 1.5;
    save_nii(inds,im,info,[fname,'_LTE'],bval,bvec);
    if (saveMode == 1) && any(brank > 2.5)
        inds = brank > 2.5;
        save_nii(inds,im,info,[fname,'_STE'],bval,bvec);
        if any(and(brank>=1.5,brank<=2.5))
            inds = and(brank>=1.5,brank<=2.5);
            save_nii(inds,im,info,[fname,'_PTE'],bval,bvec);
        end
    end
elseif saveMode == 2
    inds = bval < bthresh;
    save_nii(inds,im,info,[fname,'_b0'],bval,bvec);
    %
    inds = and(bval>=bthresh,brank<1.5);
    save_nii(inds,im,info,[fname,'_LTE'],bval,bvec);
    %
    if any(brank > 2.5)
        inds = and(bval>=bthresh,brank>2.5);
        save_nii(inds,im,info,[fname,'_STE'],bval,bvec);
    end
    %
    if any(and(brank>=1.5,brank<=2.5))
        inds = and(bval>=bthresh,and(brank>=1.5,brank<=2.5));
        save_nii(inds,im,info,[fname,'_PTE'],bval,bvec);
    end
else
    error('unknown saveMode')
end

end

function save_nii(inds,im,info,fname,bval,bvec)
    im_a = im(:,:,:,inds);
    info.ImageSize(4) = size(im_a,4);
    info.raw.dim(5) = size(im_a,4);
    niftiwrite(im_a,[fname, '.nii'],info,'Compressed',true);

    if ~isempty(bval)
        bval = bval(inds);
        fid = fopen([fname, '.bval'], 'w');
        fprintf(fid, '%d ', round(bval(1:end-1)));
        fprintf(fid, '%d', round(bval(end)));
        fclose(fid);
    end

    if ~isempty(bvec)
        bvec = bvec(:,inds);
        fid = fopen([fname, '.bvec'], 'w');
        for n=1:3
            fprintf(fid, '%.6f ', bvec(n,1:end-1));
            fprintf(fid, '%.6f\n', bvec(n,end));
        end
        fclose(fid);
    end
end