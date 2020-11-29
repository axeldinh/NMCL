function plot_SW(xc, U, xc_exact, U_exact, time)

subplot(2,1,1)
% Plot the height
plot(xc, U(1,:), '-b', 'linewidth', 2)
title("time = " + num2str(time))

if length(U_exact) > 0
    hold on
    plot(xc_exact, U_exact(1,:), '--r', 'linewidth', 2)
    hold off
    legend('Height Approximation','Exact Solution', 'Location', 'best')
else
    legend('Height Approximation', 'Location', 'best')
end

subplot(2,1,2)
%Plot the discharge
plot(xc, U(2,:), '-b', 'linewidth', 2)
title("time = " + num2str(time))
legend('Discharge Approximation', 'Location', 'best')
if length(U_exact) > 0
    hold on
    plot(xc_exact, U_exact(2,:), '--r', 'linewidth', 2)
    hold off
    legend('Discharge Approximation','Exact Solution', 'Location', 'best')
else
    legend('Discharge Approximation', 'Location', 'best')
end

pause(.01)

