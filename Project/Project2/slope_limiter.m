function dU = slope_limiter(U, h, M, limiter)

dUL = U(2:end-1) - U(1:end-2);
dUR = U(3:end) - U(2:end-1);

switch limiter
    
    case 'None'
        dU = 0*dUL;
        
    case 'MINMOD'
        
        dU = MINMOD([dUL', dUR'])';
        
    case 'MUSCL'
        dU = MINMOD([0.5*(dUL'+dUR'), 2*dUL', 2*dUR'])';
        
    case 'TVB'
        dU = TVB([0.5*(dUL'+dUR'), 2*dUL', 2*dUR'], M, h)';
end