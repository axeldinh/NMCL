% Project 1_2_a, plot solution for (5) = 2, (6) = 3
close all
clear all
clc

%% Initialize

% num_flux = 'LF'; % for project 1_2_a
num_flux = 'Roe'; % for project 1_3_c

[U0, S, a, b, bc, g] = Initial_conditions(2);
CFL = 0.5; T = 2;
N = 500;

%% Compute the solution
U = solver(U0,S,a,b,N,T,CFL,bc,num_flux);

%% Plot the solution

h = (b-a)/N;
xc = a+0.5*h:h:b-0.5*h;

figure()
subplot(2,1,1)
plot(xc, U(1,:), '-', 'LineWidth', 2)

subplot(2,1,2)
plot(xc, U(2,:), '-', 'LineWidth', 2)