%% central difference
%  major << dispersion >> error

N_max = 100;

figure( 'Name', 'Central difference causes dispersion error in advection equation' );


%% initial solution
x0 = [0:N_max-1] / N_max * 2*pi;
u0 = sin( x0/2 ) .^16;
p1 = plot( x0, u0, 'b' );

%% obvious error appears when <N> is small s.t <dx> is big
hold on;
p2 = [];
t = [];
for N = 1 : 30  % number of periods
    delete( p2 ); % to update p2
    delete( t );

    [t,u] = ode45( @dudt_central, [ 0, 2*N*pi ], u0 );

    p2 = plot( x0, u(end,:), 'r' );
    
    axis( [ 0, 7, 0, 1 ] );
    
    txt = [ 'N = ', num2str(N) ];
    t = text( 5, 0.8, txt );
  
    drawnow;
  
end

hold off;

