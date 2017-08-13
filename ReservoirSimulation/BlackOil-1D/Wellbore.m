function [f] = Wellbore( i, WI, dens, rperm, viscosity, P, Pbh, DV )

f = WI * dens(:,i) .* rperm(:,i) ./viscosity' * ( P - Pbh ) / DV;