function rhs = evalRHS(U, S, x, k, h, time, bc, num_flux, M, limiter)

N = length(U(1,:));

% Make 2 ghost nodes

U_ext = apply_bc(U, bc, 2);

% Compute the slopes
dU = zeros(2, N+2);
dU(1,:) = slope_limiter(U_ext(1,:), h, M, limiter);
dU(2,:) = slope_limiter(U_ext(2,:), h, M, limiter);

% Compute the approximation at time time+k/2
U_half = zeros(2, N+2);
for i = 1:N+2
    A = compute_grad(U_ext(:,1+i));
    U_half(:,i) = U_ext(:,1+i) - k/h * A * dU(:,i);
end

% Approximate the cell interfaces
UL = U_half(:,1:end-1) + 0.5*dU(:,1:end-1);
UR = U_half(:,2:end) - 0.5*dU(:,2:end);

% Compute the numerical flux
Flux = numerical_flux(UL, UR, num_flux);

rhs = -(1/h) * (Flux(:,2:end) - Flux(:,1:end-1));

% Add the source term
rhs = rhs + S(x,time);