function U = solver(U0,S,a,b,N,T,CFL,bc,num_flux)
% Solves the shallow water equations, with finite volume scheme and forward euler time integration
%
% Parameters
%   -U: Initial conditions (Integrated over cells)
%   -S: Source term
%   -a,b: Space limits [a,b]
%   -N: number of cells
%   -T: Final time
%   -CFL: Courant Friedrich Lewy condition
%   -bc: Bondary condition, can be 'Periodic', 'Open'
%   -num_flux: numerical flux, can be 'LF', 'Roe'

h = (b-a)/N;
xf = a:h:b;
xc = a+0.5*h:h:b-0.5*h;

% Compute cell averages

U = zeros(2,N);
for i = 1:N
    U(:,i) = integral(U0, xf(i), xf(i+1), 'ArrayValued', true, 'AbsTol', 1e-14)/h;
end

time = 0;

while time < T
    
    s = max(abs(U(2,:)./U(1,:)) + sqrt(U(1,:)));
    k = CFL*h / s;
        
    if time + k > T
        k = T-time;
    end
    
    U = U + k*evalRHS(U,S,h,time,xc,bc,num_flux);
    
    time = time +k;
end
    