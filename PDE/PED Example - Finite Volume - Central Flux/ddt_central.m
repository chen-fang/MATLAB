function [ du_bar_dt ] = ddt_central( t, u_bar )

    dx = 1. / length(u_bar);
    f_bar = u_bar.^2 / 2; % Burger's equation
    f_interface = (f_bar(2:end) + f_bar(1:end-1)) / 2;
    
    % Boundary condition
    % Assume flux f=0 at both ends
    f_interface = [ 0; f_interface; 0 ]; 
    
    du_bar_dt = (f_interface(1:end-1) - f_interface(2:end)) / dx;

end

