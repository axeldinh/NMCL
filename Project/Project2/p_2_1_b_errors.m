% Compute and plot errors project 2-1-b

clc
close all
clear all

[U0, S, a, b, bc,g] = Initial_conditions(1);
CFL = 0.5; T = 2; M = 100;
N_set = [10,20,50,100,150,300,500];

U_ex  =@(x) U0(x-a*T);
errors = zeros(16, length(N_set));
p = 2; % norm to compute error

for i = 1:length(N_set)
    N = N_set(i);
    % Compute the solutions
    U_LF_None = solver(U0,S,a,b,N,T,CFL,bc,'LF',M,'None');
    U_LF_Minmod = solver(U0,S,a,b,N,T,CFL,bc,'LF',M,'MINMOD');
    U_LF_MUSCL = solver(U0,S,a,b,N,T,CFL,bc,'LF',M,'MUSCL');
    U_LF_TVB = solver(U0,S,a,b,N,T,CFL,bc,'LF',M,'TVB');

    U_Roe_None = solver(U0,S,a,b,N,T,CFL,bc,'Roe',M,'None');
    U_Roe_Minmod = solver(U0,S,a,b,N,T,CFL,bc,'Roe',M,'MINMOD');
    U_Roe_MUSCL = solver(U0,S,a,b,N,T,CFL,bc,'Roe',M,'MUSCL');
    U_Roe_TVB = solver(U0,S,a,b,N,T,CFL,bc,'Roe',M,'TVB');
    
    h = (b-a)/N;
    xf = a:h:b;
    xc = a+0.5*h:h:b-0.5*h;
    
    % Compute the exact averages
    U_exact = zeros(2,N);
    for j = 1:N
        U_exact(:,j) = integral(U_ex, xf(j), xf(j+1), 'ArrayValued', true, 'AbsTol', 1e-14)/h;
    end
    
    % Compute the errors
    errors(1:2,i) = p_error(U_LF_None, U_exact, h, p);
    errors(3:4,i) = p_error(U_LF_Minmod, U_exact, h, p);
    errors(5:6,i) = p_error(U_LF_MUSCL, U_exact, h, p);
    errors(7:8,i) = p_error(U_LF_TVB, U_exact, h, p);
    
    errors(9:10,i) = p_error(U_Roe_None, U_exact, h, p);
    errors(11:12,i) = p_error(U_Roe_Minmod, U_exact, h, p);
    errors(13:14,i) = p_error(U_Roe_MUSCL, U_exact, h, p);
    errors(15:16,i) = p_error(U_Roe_TVB, U_exact, h, p);
end

% Plot errors
% Lax-Friedrich Flux
figure()
suptitle('Lax-Friedrich flux errors')

subplot(2,1,1)
title('Height')
loglog(N_set, N_set.^(-1), '-k', 'linewidth', 2)
hold on
loglog(N_set, errors(1,:), '--', 'linewidth', 2)
loglog(N_set, errors(3,:), '--', 'linewidth', 2)
loglog(N_set, errors(5,:), '--', 'linewidth', 2)
loglog(N_set, errors(7,:), '--', 'linewidth', 2)
xlabel('N')
ylabel('Error')
legend('O(h)','None', 'minmod', 'muscl', 'TVB', 'Location', 'best')

subplot(2,1,2)
title('Discharge')
loglog(N_set, N_set.^(-1), '-k', 'linewidth', 2)
hold on
loglog(N_set, errors(2,:), '--', 'linewidth', 2)
loglog(N_set, errors(4,:), '--', 'linewidth', 2)
loglog(N_set, errors(6,:), '--', 'linewidth', 2)
loglog(N_set, errors(8,:), '--', 'linewidth', 2)
xlabel('N')
ylabel('Error')
legend('O(h)','None', 'minmod', 'muscl', 'TVB', 'Location', 'best')


% Roe Flux
figure()
suptitle('Roe flux errors')

subplot(2,1,1)
title('Height')
loglog(N_set, N_set.^(-1), '-k', 'linewidth', 2)
hold on
loglog(N_set, errors(9,:), '--', 'linewidth', 2)
loglog(N_set, errors(11,:), '--', 'linewidth', 2)
loglog(N_set, errors(13,:), '--', 'linewidth', 2)
loglog(N_set, errors(15,:), '--', 'linewidth', 2)
xlabel('N')
ylabel('Error')
legend('O(h)','None', 'minmod', 'muscl', 'TVB', 'Location', 'best')

subplot(2,1,2)
title('Discharge')
loglog(N_set, N_set.^(-1), '-k', 'linewidth', 2)
hold on
loglog(N_set, errors(10,:), '--', 'linewidth', 2)
loglog(N_set, errors(12,:), '--', 'linewidth', 2)
loglog(N_set, errors(14,:), '--', 'linewidth', 2)
loglog(N_set, errors(16,:), '--', 'linewidth', 2)
xlabel('N')
ylabel('Error')
legend('O(h)','None', 'minmod', 'muscl', 'TVB', 'Location', 'best')