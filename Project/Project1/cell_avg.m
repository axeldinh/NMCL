function U_avg = cell_avg(U_mid, bc)

switch bc
    
    case 'Periodic'
        U_avg = ([U_mid(:,end), U_mid] + [U_mid, U_mid(:,1)])/2;
        
    case 'Open'
        U_avg = ([U_mid(:,1), U_mid] + [U_mid, U_mid(:,end)])/2;
        
end