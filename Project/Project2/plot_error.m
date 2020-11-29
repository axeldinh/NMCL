function plot_error(err, h, p, component)

switch component
    
    case 'both'
        order_height = (log(err(1,1)) - log(err(1,end)))/(log(h(1))-log(h(end)));
        order_height = floor(order_height);
        order_discharge = (log(err(2,1)) - log(err(2,end)))/(log(h(1))-log(h(end)));
        order_discharge = floor(order_discharge);

        figure()
        subplot(2,1,1)
        loglog(h, err(1,:), 'DisplayName', 'Error on height')
        hold on
        loglog(h, h.^order_height, 'DisplayName', "o(h^"+num2str(order_height)+")")
        loglog(h, h.^(order_height+1), 'DisplayName', "o(h^"+num2str(order_height+1)+")")
        xlabel('h', 'fontsize', 12)
        ylabel('Error', 'fontsize', 12)
        title("Error in norm " + num2str(p) + " of the height", 'fontsize', 14)
        legend show

        subplot(2,1,2)
        loglog(h, err(2,:), 'DisplayName', 'Error on discharge')
        hold on
        loglog(h, h.^order_discharge, 'DisplayName', "o(h^"+num2str(order_discharge)+")")
        loglog(h, h.^(order_discharge+1), 'DisplayName', "o(h^"+num2str(order_discharge+1)+")")
        xlabel('h', 'fontsize', 12)
        ylabel('Error', 'fontsize', 12)
        title("Error in norm " + num2str(p) + " of the discharge", 'fontsize', 14)
        legend show
        
    case 'height'
        
        order_height = (log(err(1,1)) - log(err(1,end)))/(log(h(1))-log(h(end)));
        order_height = floor(order_height);
        
        figure()
        loglog(h, err(1,:), 'DisplayName', 'Error on height')
        hold on
        loglog(h, h.^order_height, 'DisplayName', "o(h^"+num2str(order_height)+")")
        loglog(h, h.^(order_height+1), 'DisplayName', "o(h^"+num2str(order_height+1)+")")
        xlabel('h', 'fontsize', 12)
        ylabel('Error', 'fontsize', 12)
        title("Error in norm " + num2str(p) + " of the height", 'fontsize', 14)
        legend show
        
    case 'Discharge'
        
        order_discharge = (log(err(2,1)) - log(err(2,end)))/(log(h(1))-log(h(end)));
        order_discharge = floor(order_discharge);
        
        figure()
        loglog(h, err(2,:), 'DisplayName', 'Error on discharge')
        hold on
        loglog(h, h.^order_discharge, 'DisplayName', "o(h^"+num2str(order_discharge)+")")
        loglog(h, h.^(order_discharge+1), 'DisplayName', "o(h^"+num2str(order_discharge+1)+")")
        xlabel('h', 'fontsize', 12)
        ylabel('Error', 'fontsize', 12)
        title("Error in norm " + num2str(p) + " of the discharge", 'fontsize', 14)
        legend show
        
end