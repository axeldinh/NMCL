%Plots solutions project 2-1-b

clc
close all
clear all

[U0, S, a, b, bc,g] = Initial_conditions(1);
CFL = 0.5; T = 2; M = 1;
N = 500;

% Compute the solutions
U_LF_None = solver(U0,S,a,b,N,T,CFL,bc,'LF',M,'None');
U_LF_Minmod = solver(U0,S,a,b,N,T,CFL,bc,'LF',M,'MINMOD');
U_LF_MUSCL = solver(U0,S,a,b,N,T,CFL,bc,'LF',M,'MUSCL');
U_LF_TVB = solver(U0,S,a,b,N,T,CFL,bc,'LF',M,'TVB');

U_Roe_None = solver(U0,S,a,b,N,T,CFL,bc,'Roe',M,'None');
U_Roe_Minmod = solver(U0,S,a,b,N,T,CFL,bc,'Roe',M,'MINMOD');
U_Roe_MUSCL = solver(U0,S,a,b,N,T,CFL,bc,'Roe',M,'MUSCL');
U_Roe_TVB = solver(U0,S,a,b,N,T,CFL,bc,'Roe',M,'TVB');

%  Compute exact solutions
h = (b-a)/N;
xf = a:h:b;
xc = a+0.5*h:h:b-0.5*h;
U_ex = @(x) U0(x-T);
U_exact = zeros(2,N);
for j = 1:N
    U_exact(:,j) = integral(U_ex, xf(j), xf(j+1), 'AbsTol', 1e-14, 'ArrayValued', true)/h;
end

figure()
suptitle('Lax-friedrich flux')
subplot(2,1,1)
plot(xc,U_exact(1,:), '-k', 'linewidth', 2)
hold on
plot(xc,U_LF_None(1,:), '--', 'linewidth', 2)
plot(xc,U_LF_Minmod(1,:), '--', 'linewidth', 2)
plot(xc,U_LF_MUSCL(1,:), '--', 'linewidth', 2)
plot(xc,U_LF_TVB(1,:), '--', 'linewidth', 2)

legend('Exact', 'None', 'minmod', 'muscl', 'TVB', 'Location', 'best')

subplot(2,1,2)
plot(xc,U_exact(2,:), '-k', 'linewidth', 2)
hold on
plot(xc,U_LF_None(2,:), '--', 'linewidth', 2)
plot(xc,U_LF_Minmod(2,:), '--', 'linewidth', 2)
plot(xc,U_LF_MUSCL(2,:), '--', 'linewidth', 2)
plot(xc,U_LF_TVB(2,:), '--', 'linewidth', 2)

legend('Exact', 'None', 'minmod', 'muscl', 'TVB', 'Location', 'best')

figure()
suptitle('Roe flux')
subplot(2,1,1)
plot(xc,U_exact(1,:), '-k', 'linewidth', 2)
hold on
plot(xc,U_Roe_None(1,:), '--', 'linewidth', 2)
plot(xc,U_Roe_Minmod(1,:), '--', 'linewidth', 2)
plot(xc,U_Roe_MUSCL(1,:), '--', 'linewidth', 2)
plot(xc,U_Roe_TVB(1,:), '--', 'linewidth', 2)

legend('Exact', 'None', 'minmod', 'muscl', 'TVB', 'Location', 'best')

subplot(2,1,2)
plot(xc,U_exact(2,:), '-k', 'linewidth', 2)
hold on
plot(xc,U_Roe_None(2,:), '--', 'linewidth', 2)
plot(xc,U_Roe_Minmod(2,:), '--', 'linewidth', 2)
plot(xc,U_Roe_MUSCL(2,:), '--', 'linewidth', 2)
plot(xc,U_Roe_TVB(2,:), '--', 'linewidth', 2)

legend('Exact', 'None', 'minmod', 'muscl', 'TVB', 'Location', 'best')