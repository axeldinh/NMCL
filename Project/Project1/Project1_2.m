clc
close all
clear all

%% Project 1- We want to solve the shallow water equations d[h,m]/dt + d[m, m^2 / h + 0.5*g*h^2]/dx = S(x,t) (h*u = m)
% Eigenvalues of the flux term: lambda_i = u_i +- sqrt(h_ig)
% CFL condition k <= CFL*(dx/max(|lambda_i|))

%% Initialization

data = 1;

Plot = 0;
Plot_fine = 1;
Error = 0;

switch data
        
    case 1
        u = 0.25;
        g = 1;
        hIC = @(x) 1 - 0.1*sin(pi*x);
        mIC = @(x) 0*x;
        SIC = @(x,t) [0*x; 0*x];
        bc = 'Periodic';
        
    case 2
        u = 0.25;
        g = 1;
        hIC = @(x) 1 - 0.2*sin(pi*x);
        mIC = @(x) 0*x + 0.5;
        SIC = @(x,t) [0*x; 0*x];
        bc = 'Periodic';
end

f = @(u) [u(2,:); u(2,:).^2./u(1,:) + 0.5*g*u(1,:).^2];
CFL = 0.5;

a = 0;
b = 2;

T = 2;

err = [];

if Error
    Dx = logspace(-1,-4,10);
else
    Dx = 0.1;
end
    

%% Discretization very fine mesh

Dx_fine = 0.00001;
for dx = Dx_fine
    dt = 1;

    N = length(a:dx:b);

    x = a:dx:b;
    x_mid = a+dx/2:dx:b-dx/2;

    %% Solving

    % We take the cell averages (depends on the boundary conditions)
    U_mid = [hIC(x_mid); mIC(x_mid)];
    U_fine = cell_avg(U_mid, bc);
    S_mid = SIC(x_mid,0);
    S = cell_avg(S_mid, bc);

    time = 0;
    s = max(abs(U_fine(2)./U_fine(1)) + sqrt(g*U_fine(1)));
    
    iter = 0;

    while time < T

        dt = CFL*dx/s;

        if time + dt > T
            dt = T - time;
        end

        U_ext = apply_bc(U_fine, bc);
        Flux = LFflux(U_ext(:,1:end-1), U_ext(:,2:end),f);
        U_fine = U_fine - dt/dx * (Flux(:,2:end) - Flux(:,1:end-1));
        U_fine = U_fine + dt*S;

        time = time + dt;

        S = cell_avg(SIC(x_mid, time),bc);
        
        if Plot_fine && mod(iter,100) == 0
            subplot(2,1,1)
            plot(x,U_fine(1,:), 'DisplayName', 'Approximation')
            title("Height, time = " + num2str(time))
            pause(.01)
            
            subplot(2,1,2)
            plot(x,U_fine(2,:), 'DisplayName', 'Approximation')
            title("Discharge, time = " + num2str(time))
            pause(.01)
        end
    end
end


%% Discretization


for dx = Dx
    dt = 1;

    N = length(a:dx:b);

    x = a:dx:b;
    x_mid = a+dx/2:dx:b-dx/2;

    %% Solving

    % We take the cell averages (depends on the boundary conditions)
    U_mid = [hIC(x_mid); mIC(x_mid)];
    U = cell_avg(U_mid, bc);
    S_mid = SIC(x_mid,0);
    S = cell_avg(S_mid, bc);

    time = 0;
    s = max(abs(U(2)./U(1)) + sqrt(g*U(1)));
    
    iter = 0;

    while time < T

        dt = CFL*dx/s;

        if time + dt > T
            dt = T - time;
        end

        U_ext = apply_bc(U, bc);
        Flux = LFflux(U_ext(:,1:end-1), U_ext(:,2:end),f);
        U = U - dt/dx * (Flux(:,2:end) - Flux(:,1:end-1));
        U = U + dt*S;

        time = time + dt;

        S = cell_avg(SIC(x_mid, time),bc);
        
        if Plot && mod(iter,100) == 0
            subplot(2,1,1)
            plot(x,U(1,:), 'DisplayName', 'Approximation')
            title("Height, time = " + num2str(time))
            pause(.01)
            
            subplot(2,1,2)
            plot(x,U(2,:), 'DisplayName', 'Approximation')
            title("Discharge, time = " + num2str(time))
            pause(.01)
        end
    end
    err = [err, max(max(U-fine_to_current(U_fine, Dx_fine, dx, a, a)))];
end

if Error
    
    figure()
    loglog(Dx, err, 'DisplayName', 'Error')
    hold on
    loglog(err, err, 'DisplayName', 'o(dx)')
    xlabel('dx')
    ylabel('L inf error')
    legend show
    
end