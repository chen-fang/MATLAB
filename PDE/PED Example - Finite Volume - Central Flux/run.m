% grid
x_interface = linspace( 0, 1, 101 );
% initial shape
u_interface = sin( x_interface * 2*pi );
% take average value 
u = (u_interface(1:end-1) + u_interface(2:end)) / 2;

figure;
plot(u);

u0 = u;
[t,u] = ode45( @ddt_central, [0,0.1], u0 );
x_bar = ( x_interface(1:end-1) + x_interface(2:end) ) / 2;

figure;
plot( x_bar, u );

[t,u] = ode45( @ddt_central, [0,0.5], u0 );
figure;
plot( x_bar, u );

figure;
plot( x_bar, u(end,:) );