% Plot project 2_2_3

close all
clear all
clc

%% Initialization

[U0, S, a, b, bc, g] = Initial_conditions(4);
T = 0.5; CFL = 0.5; M = 500;
N = 500; N_ref = 5000;

h = (b-a)/N;
xc = a+0.5*h:h:b-0.5*h;

h_ref = (b-a)/N_ref;
xc_ref = a+0.5*h_ref:h_ref:b-0.5*h_ref;

%% Create reference solution and rescale it

U_ref = solver(U0,S,a,b,N_ref,T,CFL,bc,'LF',M,'None');
U_ref = ref_to_current(U_ref,xc_ref,xc);

%% Solve for the different numerical fluxes and slope limiters

U_LF_None = solver(U0,S,a,b,N,T,CFL,bc,'LF',M,'None');
U_LF_Minmod = solver(U0,S,a,b,N,T,CFL,bc,'LF',M,'MINMOD');
U_LF_MUSCL = solver(U0,S,a,b,N,T,CFL,bc,'LF',M,'MUSCL');
U_LF_TVB = solver(U0,S,a,b,N,T,CFL,bc,'LF',M,'TVB');

U_Roe_None = solver(U0,S,a,b,N,T,CFL,bc,'Roe',M,'None');
U_Roe_Minmod = solver(U0,S,a,b,N,T,CFL,bc,'Roe',M,'MINMOD');
U_Roe_MUSCL = solver(U0,S,a,b,N,T,CFL,bc,'Roe',M,'MUSCL');
U_Roe_TVB = solver(U0,S,a,b,N,T,CFL,bc,'Roe',M,'TVB');

%% Plot the solutions

figure()
suptitle('Lax-friedrich flux')
subplot(2,1,1)
plot(xc,U_ref(1,:), '-k', 'linewidth', 2)
hold on
plot(xc,U_LF_None(1,:), '--', 'linewidth', 2)
plot(xc,U_LF_Minmod(1,:), '--', 'linewidth', 2)
plot(xc,U_LF_MUSCL(1,:), '--', 'linewidth', 2)
plot(xc,U_LF_TVB(1,:), '--', 'linewidth', 2)

legend('Exact', 'None', 'minmod', 'muscl', 'TVB', 'Location', 'best')

subplot(2,1,2)
plot(xc,U_ref(2,:), '-k', 'linewidth', 2)
hold on
plot(xc,U_LF_None(2,:), '--', 'linewidth', 2)
plot(xc,U_LF_Minmod(2,:), '--', 'linewidth', 2)
plot(xc,U_LF_MUSCL(2,:), '--', 'linewidth', 2)
plot(xc,U_LF_TVB(2,:), '--', 'linewidth', 2)

legend('Exact', 'None', 'minmod', 'muscl', 'TVB', 'Location', 'best')

figure()
suptitle('Roe flux')
subplot(2,1,1)
plot(xc,U_ref(1,:), '-k', 'linewidth', 2)
hold on
plot(xc,U_Roe_None(1,:), '--', 'linewidth', 2)
plot(xc,U_Roe_Minmod(1,:), '--', 'linewidth', 2)
plot(xc,U_Roe_MUSCL(1,:), '--', 'linewidth', 2)
plot(xc,U_Roe_TVB(1,:), '--', 'linewidth', 2)

legend('Exact', 'None', 'minmod', 'muscl', 'TVB', 'Location', 'best')

subplot(2,1,2)
plot(xc,U_ref(2,:), '-k', 'linewidth', 2)
hold on
plot(xc,U_Roe_None(2,:), '--', 'linewidth', 2)
plot(xc,U_Roe_Minmod(2,:), '--', 'linewidth', 2)
plot(xc,U_Roe_MUSCL(2,:), '--', 'linewidth', 2)
plot(xc,U_Roe_TVB(2,:), '--', 'linewidth', 2)

legend('Exact', 'None', 'minmod', 'muscl', 'TVB', 'Location', 'best')