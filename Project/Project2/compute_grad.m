function grad = compute_grad(U)
% compute the gradient of the flux
g = 1;
grad = [0,1;-(U(2)/U(1))^2 + g*U(1), 2*(U(2)/U(1))];