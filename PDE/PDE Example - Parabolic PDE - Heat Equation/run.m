n = 100;

k = 1;
dx = 1/n;

global A b
A = zeros( n-1 );
b = zeros( n-1, 1 );

for i = 1 : n-1
    A(i,i) = -2*k/dx.^2;
end

for i = 1 : n-2
    A(i,i+1) = k/dx.^2;
    A(i+1,i) = k/dx.^2;
end


% boundary conditions
u0 = 1;
un = -1;


b(1) = k/dx.^2 * u0;
b(end) = k/dx.^2 * un;

% solution vector
u_0 = zeros( n-1, 1 );

% PDE solver
% ode45( function, time range, initial value )
[t,u] = ode45( @dudt, [0,1], u_0 ); 



figure( 'Name', 'Heat Equation u(x)' )
% Animation
for i = 1 : length(u)
    % plot
    plot( u(i,:) );
    % fix axis (after 'plot' to function properly)
    axis( [ 0 100 -1 1 ] );
    % show text
    txt = [ 't = ', num2str(t(i)) ];
    text( 70, 0.8, txt );
    % draw
    drawnow  
end