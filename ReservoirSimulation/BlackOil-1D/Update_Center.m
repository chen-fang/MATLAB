function [dens, rperm, dens_DP, rperm_DS] = Update_Center( Ncell, OIL, WAT, P, Sw, Co, constant_water_density )

%     ps( Sw );
for i = 1 : Ncell

    
    dens( OIL,i) = Density_Oil( Co, P(i) );
    dens( WAT,i) = constant_water_density;
    rperm( OIL,i) = RPerm_Oil( Sw(i) );
    rperm( WAT,i) = RPerm_Wat( Sw(i) );
    
%     p( rperm( OIL, i ) );
    
    dens_DP( OIL, i ) = Density_Oil_DP( Co, P(i) );
    dens_DP( WAT, i ) = 0;
    
    rperm_DS( OIL, i ) = RPerm_Oil_DS( Sw(i) );
    rperm_DS( WAT, i ) = RPerm_Wat_DS( Sw(i) );
end