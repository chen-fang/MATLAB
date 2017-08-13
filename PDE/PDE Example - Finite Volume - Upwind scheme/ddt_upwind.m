function [ dudt ] = ddt_upwind( t, u )

    % use Burger's equation as an example
    f = u.^2 / 2;
    % shock speed
    ss = ( u(2:end) + u(1:end-1) ) / 2;
    % 
    f_interface = [0]; % left B.C
    for i = 2 : length( u )
        if ss(i-1) > 0
            f_interface = [ f_interface; f( i-1 ) ];
        else
            f_interface = [ f_interface; f( i ) ];
        end
    end
    % right B.C
    f_interface = [ f_interface; 0 ];
    %
    dx =  1 / length(u);
    %
    dudt = ( f_interface(1:end-1) - f_interface(2:end) ) / dx;
end

