% Project 1_4

close all
clear all
clc

%% Initialization
[U0, S, a, b, bc, g] = Initial_conditions(4);
CFL = 0.5; T = 0.1;

N_ref = 10000;
N = 1000;

%% Compute reference solution

U_ref = solver(U0,S,a,b,N_ref,T,CFL,bc,'LF');

%% Compute solutions

U_LF = solver(U0,S,a,b,N,T,CFL,bc,'LF');
U_Roe = solver(U0,S,a,b,N,T,CFL,bc,'Roe');

%% Plot the solutions

h_ref = (b-a)./N_ref; h = (b-a)./N;
xc_ref = a+0.5*h_ref:h_ref:b-0.5*h_ref;
xc = a+0.5*h:h:b-0.5*h;

figure()

subplot(2,1,1)
plot(xc_ref, U_ref(1,:), '-k', 'LineWidth', 2)
hold on
plot(xc, U_LF(1,:), '--', 'LineWidth', 2)
plot(xc, U_Roe(1,:), '--', 'LineWidth', 2)
legend('Exact', 'LF', 'Roe')
title('Height')

subplot(2,1,2)
plot(xc_ref, U_ref(2,:), '-k', 'LineWidth', 2)
hold on
plot(xc, U_LF(2,:), '--', 'LineWidth', 2)
plot(xc, U_Roe(2,:), '--', 'LineWidth', 2)
legend('Exact', 'LF', 'Roe')
title('Discharge')

