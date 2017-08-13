% << upwind >> difference instead of << central >> diffenrence is used
% to achieve stability. See notes.

function [ f ] = dudt_upwind( t, u )

    N = length(u);
    dx = 2 * pi / N;
    
    % the below vectorizing method will be faster
    % than the for-loop
    %
    % u(end) == u(0) as we assume periodic function
    %
    
    %% U > 0 with forward spatial discretization
    U = 1;
    dudx = [ ( u(1)-u(end) )/dx; ... 
        ( u(2:end)-u(1:end-1) ) /dx ];

    %% U < 0 with backward spatial discretization
%     U = -1;
%     dudx = [ ( u(end)-u(1) )/dx; ...
%         ( u(1:end-1)-u(2:end) ) /dx ];
%     

    
    % problem unsolved...
    % Analytical solution should move towards left
    % when U < 0, but this numerical result shows the
    % solution still moves towards right, same as the
    % case with U > 0. Why???
    
    %% Unstable case
    % U < 0 with forward spatital discretization  
%     U = -1;
%     dudx = [ ( u(1)-u(end) )/dx; ... 
%         ( u(2:end)-u(1:end-1) ) /dx ];
    
    %%
    f = -U * dudx;end

