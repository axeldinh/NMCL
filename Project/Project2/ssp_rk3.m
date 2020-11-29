function U1 = ssp_rk3(U, S, x, k, h, time, bc, num_flux, M, limiter)

k1 = U + k*evalRHS(U, S, x, k, h, time, bc, num_flux, M, limiter);
k2 = 0.75*U + 0.25*(k1 + k*evalRHS(k1, S, x, k, h, time+k, bc, num_flux, M, limiter));
U1 = 1/3 * U + 2/3 * (k2 + k*evalRHS(k2, S, x, k, h, time+k/2, bc, num_flux, M, limiter));