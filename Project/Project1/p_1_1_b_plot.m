% Project 1_1_b / 1_3_c plot

close all
clear all
clc

%% Initialize

%num_flux = 'LF'; % for project 1_1_b
num_flux = 'Roe'; % for project 1_3_c

[U0, S, a, b, bc, g] = Initial_conditions(1);

CFL = 0.5; T = 2;
N = 1000;
U_ex = @(x) U0(x-a*T);

%% Compute the solution
U = solver(U0,S,a,b,N,T,CFL,bc,num_flux);

%% Compute the exact solution
h = (b-a)/N;
xf = a:h:b;
xc = a+0.5*h:h:b-0.5*h;

U_exact = zeros(2,N);
for i = 1:N
    U_exact(:,i) = integral(U_ex, xf(i), xf(i+1), 'AbsTol', 1e-14, 'ArrayValued', true)/h;
end

%% Plot the solutions

figure()
suptitle("Lax-Friedrich Flux, " + num2str(N) + "cells")

subplot(2,1,1)
plot(xc, U_exact(1,:), '-k', 'LineWidth', 2)
hold on
plot(xc, U(1,:), '--', 'LineWidth', 2)
ylabel('Height')

subplot(2,1,2)
plot(xc, U_exact(2,:), '-k', 'LineWidth', 2)
hold on
plot(xc, U(2,:), '--', 'LineWidth', 2)
ylabel('Discharge')



