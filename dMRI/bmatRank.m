function rankOut = bmatRank(bmat)
% Approximates the rank of a list of b-matrices. Dim 1 should be b-matrix
% entries (s/mm2), dim 2 should be number of acquisitions. 
% %  6 matrix entries makes symmetric matrices be assumed. First three are
% diagonal (xx, yy, zz), next are cross terms (xy, xz, yz).

% Set threshold used for eigenvalues (s/mm^2)
eigThresh = 5;

rankOut = zeros(1,size(bmat,2),'like',bmat);
for n=1:size(bmat,2)
    if size(bmat,1)==6
        bmat_a = [bmat(1,n), bmat(4,n), bmat(5,n);
                  bmat(4,n), bmat(2,n), bmat(6,n);
                  bmat(5,n), bmat(6,n), bmat(3,n)];
    else
        error('unkown size of dim 1')
    end

    % Find eigenvalues
    e = eig(bmat_a);

    % Find rank
    rankOut(n) = sum(e>eigThresh);
    
    % Account for b=0 scans
    if rankOut(n) < 1
        rankOut(n) = 1;
    end
end