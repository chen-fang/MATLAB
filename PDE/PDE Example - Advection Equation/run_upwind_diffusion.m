%% upwind difference
%  major << diffusion >> error
%  assume periodic function [ 0, 2pi ]

%%
% N_max is chosen as a big number
% s.t the numerial solution is close to the analytical solution
N_max = 1000;  

figure( 'Name', 'Dissipation error in advection equation' );

%%
% initial solution
x0 = [0:N_max-1] / N_max * 2*pi;
u0 = sin( x0/2 ) .^16;
plot( x0, u0, 'b' );

%% obvious error appears when <N> is small s.t <dx> is big
hold on;
p2 = [];
t = [];
for N = 100 : 250
    delete( p2 ); % to update p2
    delete( t );
    
    x = [0:N-1] / N * 2*pi;
    u00 = sin( x/2 ) .^16;
    [t,u] = ode45( @dudt_upwind, [ 0, 2*pi ], u00 );

    p2 = plot( x, u(end,:), 'r' );
    
    axis( [ 0, 7, 0, 1 ] );
    
    txt = [ 'N = ', num2str(N) ];
    t = text( 5, 0.8, txt );
  
    drawnow;
  
end

hold off;


