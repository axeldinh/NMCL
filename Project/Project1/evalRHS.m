function rhs = evalRHS(U, S, h, time, xc, bc, num_flux)
% Computes the right hand side of the equation
% Here the temporal ode is solved with forward Euler

% Make 1 ghost node
U_ext = apply_bc(U, bc, 1);

% Compute the numerical_flux
UL = U_ext(:,1:end-1);
UR = U_ext(:,2:end);
Flux = numerical_flux(UL, UR, num_flux);

% Compute the Rhs
rhs = -1/h * (Flux(:,2:end) - Flux(:,1:end-1));

% Add the source term
rhs = rhs + S(xc, time);
end
