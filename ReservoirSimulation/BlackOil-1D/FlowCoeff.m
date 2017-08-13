function [f] = FlowCoeff( OIL, WAT, perm, viscosity, DX )

f( OIL, 1 ) = perm/viscosity(OIL)/DX^2;
f( WAT, 1 ) = perm/viscosity(WAT)/DX^2;