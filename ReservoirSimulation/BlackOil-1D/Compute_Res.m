function [Res] = Compute_Res( Res, Ncell, Nface, OIL, WAT, P, Sw, Sw_old, P_old, dens, densf, rperm, rpermf, DP, flowcoeff, Co, constant_water_density, porosity, viscosity, DT, WI, wid_inj, const_Qw_inj, wid_prd, BHP_PRD, DV, print_flag )

Res = zeros( length(Res), 1 );

% if print_flag == 1
%     pp( P );
%     ps( Sw );
% end

%% accum
for i = 1 : Ncell    
    dens_old( OIL ,i ) = Density_Oil( Co, P_old(i) );
    dens_old( WAT, i ) = constant_water_density;
    
    accum_old = ResAccum( i, OIL, WAT, Sw_old, dens_old );
    accum_new = ResAccum( i, OIL, WAT, Sw, dens );
    
    Res(2*i-1,1) = Res(2*i-1,1) + porosity/DT*( accum_new(OIL) - accum_old(OIL) );
    Res(2*i,1)   = Res(2*i,1)   + porosity/DT*( accum_new(WAT) - accum_old(WAT) );
    

    if print_flag == 1 && i == 1
        tmp_oil_old = porosity/DT* accum_old(OIL);
        p( tmp_oil_old );        
        
        tmp_oil = porosity/DT* accum_new(OIL);
        p( tmp_oil );
    end
%      

    
%     p( Res( 2*i-1, 1 ) );
    %print( accum_new(OIL) )
    %print( Res(2*i-1,1) )
end

%% flow
for j = 1 : Nface
    lcid = j;
    rcid = j+1;
    
    flow = ResFlow( j, OIL, WAT, densf, rpermf, DP );   

    Res( lcid*2-1 : lcid*2, 1 ) = Res( lcid*2-1 : lcid*2, 1 ) - flowcoeff .* flow;
    Res( rcid*2-1 : rcid*2, 1 ) = Res( rcid*2-1 : rcid*2, 1 ) + flowcoeff .* flow;
    
%     p( densf(OIL,j) );
    %flow
%      tmp = flowcoeff .* flow;
%      p( tmp(WAT) );
end

%% producers with constant BHP
for well_id = 1 : length( wid_prd )
    cell_id = wid_prd(well_id);
    Q = Wellbore( cell_id, WI, dens, rperm, viscosity, P(cell_id), BHP_PRD(well_id), DV );
    Res( 2*cell_id-1, 1 ) = Res( 2*cell_id-1, 1 ) + Q( OIL );
    Res( 2*cell_id,   1 ) = Res( 2*cell_id,   1 ) + Q( WAT );
    
    %print( Q(WAT) )
end


%% Injectors with constant Qw - ONLY WATER
for well_id = 1 : length( wid_inj )
    cell_id = wid_inj( well_id ); 

    % add to cell residual
    Res( 2*cell_id, 1 ) = Res( 2*cell_id, 1 ) + const_Qw_inj / DV;
    
    % new residual equation - water only 
    new_res_id  = 2 * Ncell + well_id;
 
    Q = Wellbore( cell_id, WI, dens, rperm, viscosity, P(cell_id), P(Ncell+well_id), DV );
    Res( new_res_id, 1 ) = Q(WAT) - const_Qw_inj / DV;
end
    
if print_flag == 1
%     Res( new_res_id, 1 )
    %p( Res( new_res_id, 1 ) );
end

