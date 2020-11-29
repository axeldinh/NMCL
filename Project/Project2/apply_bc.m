function  U_ext = apply_bc(U, bc, n)
% Add n ghosts nodes

switch bc
    case 'Periodic'
        U_ext = [U(:,end-n+1:end), U, U(:,1:n)];
    case 'Open'
        U_ext = [repmat(U(:,1),1,n), U, repmat(U(:,end),1,n)];
end