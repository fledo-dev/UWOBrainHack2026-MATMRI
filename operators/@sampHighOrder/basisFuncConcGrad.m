function bfunc = basisFuncConcGrad(obj,n)
    % Outputs basis function for index n of concomitant gradient terms
    % n = index of phs_conc in sampHighOrder.m
    x = obj.phs_grid.x;
    y = obj.phs_grid.y;
    z = obj.phs_grid.z;
    switch n
    case 1
        bfunc = z.*z;
    case 2
        bfunc = x.*x + y.*y;
    case 3                
        bfunc = x.*z;
    case 4                
        bfunc = y.*z;
    otherwise
        error('Basis function undefined.')
    end
end