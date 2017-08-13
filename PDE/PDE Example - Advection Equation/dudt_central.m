% central diffenrence
% to achieve stability. See notes.

function [ f ] = dudt_central( t, u )

    N = length(u);
    dx = 2 * pi / N;
    
    % the below vectorizing method will be faster
    % than the for-loop
    %
    % u(end) == u(0) as we assume periodic function
    %
    
    %% U > 0 with forward spatial discretization
    U = 1;
    dudx = [ ( u(2)-u(end) ) /(2*dx); ... 
        ( u(3:end)-u(1:end-2) ) /(2*dx); ...
        ( u(1) - u(end-1) ) /(2*dx) ];
    
    %%
    f = -U * dudx;end

