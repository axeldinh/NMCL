function flux = GodunovFlux(UL, UR,f,g)
% Computes the Godunov flux with Roe linearization

% Compute the Roe linearization
zl = UL/sqrt(UL(1));
zr = UR/sqrt(UR(1));
z_bar = (zr + zl)/2;
u_bar = z_bar(2)/z_bar(1);
c_bar = ((zl(1)^2 + zr(1)^2)/2) * g;

D = [u_bar + sqrt(c_bar), 0; 0, u_bar - sqrt(c_bar)];
S = [1,1; u_bar + sqrt(c_bar), u_bar - sqrt(c_bar)];
Sinv = 1/(-2*sqrt(c_bar)) * [u_bar-sqrt(c_bar), -1; -u_bar-sqrt(c_bar), 1];
abs_A = S*abs(D)*Sinv;

flux = 0.5 * (f(UL) + f(UR)) - 0.5*abs_A*(UR-UL);

end