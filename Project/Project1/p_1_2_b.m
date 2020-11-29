% Project 1_2_b/1_3_c

%% Initialization

%num_flux = 'LF'; % for project 1_2_b
num_flux = 'Roe'; % for project 1_3_c

[U0, S, a, b, g] = Initial_conditions(2);
CFL = 0.5; T = 2;
N_set = [10,20,50,100,500, 1000];
N_ref = 5000; h_ref = (b-a)/N_ref; xc_ref = a+0.5*h_ref:h_ref:b-0.5*h_ref;

errors = zeros(2,length(N_set));
p = 2; % to compute p norm

%% Compute reference solution

U_ref = solver(U0,S,a,b,N_ref,T,CFL,bc,num_flux);

%% Compute the solutions and errors

for i = 1:length(N_set)
    
    N = N_set(i);
    
    U = solver(U0,S,a,b,N,T,CFL,bc,num_flux);
    
    h = (b-a)/N;
    xc = a+0.5*h:h:b-0.5*h;
    errors(:,i) = p_error(U, ref_to_current(U_ref,xc_ref,xc), h, p);
    
end

%% Plot Errors

figure()
suptitle('Errors')

subplot(2,1,1)
loglog(N_set,N_set.^(-1), '-k', 'LineWidth', 2)
hold on
loglog(N_set, errors(1,:), '--', 'LineWidth', 2)
legend('Exact', 'LF')
ylabel('Height')

subplot(2,1,2)
loglog(N_set,N_set.^(-1), '-k', 'LineWidth', 2)
hold on
loglog(N_set, errors(2,:), '--', 'LineWidth', 2)
legend('Exact', 'LF')
ylabel('Discharge')



