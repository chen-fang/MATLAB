function [f] = ResAccum( i, OIL, WAT, Sw, dens )

f( OIL ) = dens( OIL, i ) * ( 1-Sw(i) );
f( WAT ) = dens( WAT, i ) * Sw(i);