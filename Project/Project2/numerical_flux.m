function Flux = numerical_flux(UL,UR,num_flux)

g = 1;
f = @(u) [u(2,:) ; (u(2,:).^2)./u(1,:) + 0.5*g*u(1,:).^2];

switch num_flux
    
    case 'LF' % Compute the Lax-Friedrich numerical flux

        % We compute the local maximum speed
        alpha = max(abs(UL(2,:)./UL(1,:)) + sqrt(g*UL(1,:)), abs(UR(2,:)./UR(1,:)) + sqrt(g*UR(1,:)));
        Flux = 0.5*(f(UL)+f(UR)) - 0.5*alpha.*(UR-UL);
    
    case 'Roe' % Compute the Godunov flux with Roe linearization
        N = length(UL(1,:));
        Flux = zeros(2,N);
        for i = 1:N % As the Roe matrix is defined for each cell, we do a loop
            Flux(:,i) = GodunovFlux(UL(:,i), UR(:,i), f,g);
        end
end
        
        
        