function flux = LFflux(u,v,f)

g = 1;
alpha = max(abs(u(2,:)./u(1,:)) + sqrt(g*u(1,:)), abs(v(2,:)./v(1,:)) + sqrt(g*v(1,:)));

flux = 0.5*(f(u) + f(v)) - 0.5*alpha.*(v-u);