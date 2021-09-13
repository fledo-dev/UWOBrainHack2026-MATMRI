function maxEig = powermethod(Ain,x0,NitMax,doAtA,xdiffThresh)
% Use the power method to find the largest eigenvalue of either the matrix
% A or the matrix A'*A. The latter case is the default. The input A can be
% a function handle.

if nargin<4 || isempty(doAtA)
    doAtA = 1;
end
if nargin<5 || isempty(xdiffThresh)
    xdiffThresh = 1e-5;
end

% Account for different ways of supplying A
if isa(Ain,'function_handle')
    A = Ain;
else
    A = @(x,transp) Asub(x,transp,Ain);
end

% Initialize
finished = 0;
nit = 0;
x = x0;
x_prev = x0;

% Iteratively find eigenvector associated with maximum eigenvalue
while ~finished
    if doAtA
        x = A(A(x,'transp'),'notransp');
    else
        x = A(x,'notransp');
    end
    x = x/norm(x(:));
    
    % Track progress
    xdiff = norm(x(:) - x_prev(:));
    x_prev = x;
    
    nit = nit+1;
    if xdiff < xdiffThresh
        finished = 1;
    end
    if nit >= NitMax
        finished = 1;
    end
end

% Use Rayleigh quotient to find eigenvalue
if doAtA
    AtAx = A(A(x,'transp'),'notransp');
else
    AtAx = A(x,'notransp');
end
maxEig = x(:)' * AtAx(:) / (x(:)'*x(:));

end