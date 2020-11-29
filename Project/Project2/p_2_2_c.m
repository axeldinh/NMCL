% Errors for project 2_2_c, for initial conditions 2=(6), 3=(7)
close all
clear all
clc

[U0, S, a, b, bc,g] = Initial_conditions(2);
CFL = 0.5; T = 2; M = 1;
N_set = [10,20,50,100,150,300,500];
N_ref = 1000;

errors = zeros(16, length(N_set));
p = 2;

%Compute Reference solution
U_ref = solver(U0,S,a,b,N_ref,T,CFL,bc,'LF',M,'None');
h_ref = (b-a)/N_ref;
xc_ref = a+0.5*h_ref:h_ref:b-0.5*h_ref;

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
    
    % Compute the errors (need to rescale the reference solution)
    h = (b-a)/N;
    xc = a+0.5*h:h:b-0.5*h;
    
    errors(1:2,i) = p_error(U_LF_None, ref_to_current(U_ref,xc_ref,xc), h, p);
    errors(3:4,i) = p_error(U_LF_Minmod, ref_to_current(U_ref,xc_ref,xc), h, p);
    errors(5:6,i) = p_error(U_LF_MUSCL, ref_to_current(U_ref,xc_ref,xc), h, p);
    errors(7:8,i) = p_error(U_LF_TVB, ref_to_current(U_ref,xc_ref,xc), h, p);
    
    errors(9:10,i) = p_error(U_Roe_None, ref_to_current(U_ref,xc_ref,xc), h, p);
    errors(11:12,i) = p_error(U_Roe_Minmod, ref_to_current(U_ref,xc_ref,xc), h, p);
    errors(13:14,i) = p_error(U_Roe_MUSCL, ref_to_current(U_ref,xc_ref,xc), h, p);
    errors(15:16,i) = p_error(U_Roe_TVB, ref_to_current(U_ref,xc_ref,xc), h, p);
    
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