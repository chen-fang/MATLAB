% grid
x_interface = linspace( 0, 1, 101 );
%
x = ( x_interface(1:end-1) + x_interface(2:end) ) / 2;
% initial shape
u = sin( x *2*pi );
u0 = u;

% figure;
% plot( x, u0 );

% solve
[t,u] = ode45( @ddt_upwind, [0,1], u0 );
figure;
plot( x, u );
figure;
plot( x, u(end,:) );