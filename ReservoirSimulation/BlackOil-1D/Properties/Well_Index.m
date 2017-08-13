function [f] = WI( DX, DZ, perm )

re = 0.2 * DX;
rw = 0.5; % [ m ]

f = 2*pi*perm*DZ / log(re/rw);