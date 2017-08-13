N = 100;
x = [0:N-1] / N * 2*pi;  % [ 0, 100 )

phi0 = sin(x/2).^16;
psi0 = phi0 * 0;

% [t, phi_psi] = ode45( @ddt_wave, [0,pi], [phi0,psi0] );
% 
% plot( x, phi_psi( end, 1:100 ) );

figure;
phi_psi0 = [phi0,psi0];
for i = 1 : 1
    [ t, phi_psi ] = ode45( @ddt_wave, [0,0.1], phi_psi0 );
    hold off;
    phi_psi0 = phi_psi( end, : ); % new I.C. for the next time step
    plot( x, phi_psi0(1:100), 'b' );
    hold on;
    plot( x, phi_psi0(101:200), 'r' );
    pause( 0.05 );
end
