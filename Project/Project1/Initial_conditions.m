function [U0, S, a, b, bc, g] = Initial_conditions(data)

g = 1;
a = 0;
b = 2;

switch data
    
    case 1
        u = 0.25;
        hIC = @(x) 1+0.5*sin(pi*x);
        mIC = @(x) 0.25*hIC(x);
        S = @(x,t) [pi/2 * (u-1) * cos(pi*(x-t)); pi/2 * cos(pi*(x-t)) .* (-u + u^2 + g*hIC(x-t))];
        bc = 'Periodic';
        
    case 2
        hIC = @(x) 1-0.1*sin(pi*x);
        mIC = @(x) 0*x;
        S = @(x,t) [0*x;0*x];
        bc = 'Periodic';
        
    case 3
        hIC = @(x) 1-0.2*sin(2*pi*x);
        mIC = @(x) 0*x;
        S = @(x,t) [0*x;0*x];
        bc = 'Periodic';
        
    case 4
        hIC = @(x) ones(size(x));
        mIC = @(x) -1.5*(x<1);
        S = @(x,t) [0*x;0*x];
        bc = 'Open';
end

U0 = @(x) [hIC(x);mIC(x)];