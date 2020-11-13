clc
close all
clear all

%% Project 1- We want to solve the shallow water equations d[h,m]/dt + d[m, m^2 / h + 0.5*g*h^2]/dx = S(x,t) (h*u = m)
% Eigenvalues of the flux term: lambda_i = u_i +- sqrt(h_ig)
% CFL condition k <= CFL*(dx/max(|lambda_i|))

%% Initialization

data = 1;

Plot = 0;
Error = 1;

switch data
    
    case 1
        u = 0.25;
        g = 1;
        hIC = @(x) 1 + 0.5*sin(pi*x);
        mIC = @(x) u*hIC(x);
        SIC = @(x,t) [pi/2 * (u-1) * cos(pi*(x-t)); pi/2 * cos(pi*(x-t)) .* (-u + u^2 + g*hIC(x-t))];
        bc = 'Periodic';
        
    case 2
        u = 0.25;
        g = 1;
        hIC = @(x) 1 + 0.5*sin(pi*x);
        mIC = @(x) u*hIC(x);
        SIC = @(x,t) [pi/2 * (u-1) * cos(pi*(x-t)); pi/2 * cos(pi*(x-t)) .* (-u + u^2 + g*hIC(x-t))];
        bc = 'Open';
end

f = @(u) [u(2,:); u(2,:).^2./u(1,:) + 0.5*g*u(1,:).^2];
CFL = 0.5;

U_exact = @(x,t) [hIC(x-t); u*hIC(x-t)];

a = 0;
b = 2;

T = 2;

err = [];

if Error
    Dx = [0.1, 0.01, 0.001, 0.0001];
else
    Dx = 0.01;
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
            U_ex = U_exact(x,time);
            subplot(2,1,1)
            xlim([0,2])
            ylim([0,1.5])
            plot(x,U(1,:), 'DisplayName', 'Approximation')
            hold on
            plot(x, U_ex(1,:), '--', 'DisplayName', 'Exact')
            title("Height, time = " + num2str(time))
            pause(.01)
            hold off
            
            subplot(2,1,2)
            plot(x,U(2,:), 'DisplayName', 'Approximation')
            hold on
            xlim([0,2])
            ylim([0.1,0.4])
            plot(x, U_ex(2,:), '--', 'DisplayName', 'Exact')
            title("Discharge, time = " + num2str(time))
            pause(.01)
            hold off
        end



    end
    err = [err, max(max(U-U_exact(x,time)))];
end

if Error
    
    figure()
    loglog(Dx, err, 'DisplayName', 'Error')
    hold on
    loglog(err, err, 'DisplayName', 'o(dx)')
    xlabel('dx')
    ylabel('L_inf error')
    legend show
    
end