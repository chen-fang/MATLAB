function [J] = Compute_Jacobian( J, Ncell, Nface, OIL, WAT, L, R, P, Sw, dens, densf, dens_DP, rperm, rpermf, rperm_DS, DP, flowcoeff, porosity, viscosity, DT, WI, wid_inj, wid_prd, BHP_PRD, DV )

J = zeros( length(J), length(J) );

%% flow
for j = 1 : Nface
    lcid = j;
    rcid = j+1;       
    
    % W.R.T [ Pressure ] ----------------------------------
    tmp = UpstreamPTN( L, R, DP(j), OIL, dens_DP(OIL,lcid), dens_DP(OIL,rcid), WAT, dens_DP(WAT,lcid), dens_DP(WAT,rcid) );

    jflowDP_toL = flowcoeff .* rpermf(:,j) .* ( tmp(:,L) * DP(j) - densf(:,j) );
    jflowDP_toR = flowcoeff .* rpermf(:,j) .* ( tmp(:,R) * DP(j) + densf(:,j) );
    
%     jflowDP_toL( OIL, 1 )
%     jflowDP_toL( WAT, 1 )
%     jflowDP_toR( OIL, 1 )
%     jflowDP_toR( WAT, 1 )

    % current cell to center
    J( JID(OIL,lcid):JID(WAT,lcid), JID(OIL,lcid) ) = J( JID(OIL,lcid):JID(WAT,lcid), JID(OIL,lcid) ) - jflowDP_toL;
    % current cell to right
    J( JID(OIL,lcid):JID(WAT,lcid), JID(OIL,rcid) ) = J( JID(OIL,lcid):JID(WAT,lcid), JID(OIL,rcid) ) - jflowDP_toR;   
    % next cell to its center
    J( JID(OIL,rcid):JID(WAT,rcid), JID(OIL,rcid) ) = J( JID(OIL,rcid):JID(WAT,rcid), JID(OIL,rcid) ) + jflowDP_toR;
    % next cell to its left
    J( JID(OIL,rcid):JID(WAT,rcid), JID(OIL,lcid) ) = J( JID(OIL,rcid):JID(WAT,rcid), JID(OIL,lcid) ) + jflowDP_toL;
    
    % W.R.T [ Saturation ] --------------------------------
    tmp = UpstreamPTN( L, R, DP(j), OIL, rperm_DS(OIL,lcid), rperm_DS(OIL,rcid), WAT, rperm_DS(WAT,lcid), rperm_DS(WAT,rcid) );
    
    jflowDS_toL = flowcoeff .* densf(:,j) .* tmp(:,L) * DP(j);
    jflowDS_toR = flowcoeff .* densf(:,j) .* tmp(:,R) * DP(j);
    
%     [ jflowDS_toL( OIL, 1 ), jflowDS_toL( WAT, 1 ) ]   
%     [ jflowDS_toR( OIL, 1 ), jflowDS_toR( WAT, 1 ) ]
    
    
    % current cell to center
    J( JID(OIL,lcid):JID(WAT,lcid), JID(WAT,lcid) ) = J( JID(OIL,lcid):JID(WAT,lcid), JID(WAT,lcid) ) - jflowDS_toL;
    % current cell to right
    J( JID(OIL,lcid):JID(WAT,lcid), JID(WAT,rcid) ) = J( JID(OIL,lcid):JID(WAT,lcid), JID(WAT,rcid) ) - jflowDS_toR;
    % next cell to its center
    J( JID(OIL,rcid):JID(WAT,rcid), JID(WAT,rcid) ) = J( JID(OIL,rcid):JID(WAT,rcid), JID(WAT,rcid) ) + jflowDS_toR;
    % next cell to its left
    J( JID(OIL,rcid):JID(WAT,rcid), JID(WAT,lcid) ) = J( JID(OIL,rcid):JID(WAT,rcid), JID(WAT,lcid) ) + jflowDS_toL;
end

%% accum
for i = 1 : Ncell
    % W.R.T [ Pressure ]    
    J( JID(OIL,i), JID(OIL,i) ) = J( JID(OIL,i), JID(OIL,i) ) + porosity/DT * (1-Sw(i)) * dens_DP(OIL,i);
    J( JID(WAT,i), JID(OIL,i) ) = J( JID(WAT,i), JID(OIL,i) ) + porosity/DT * Sw(i)     * dens_DP(WAT,i);
  
    % W.R.T [ Saturation ]
    J( JID(OIL,i), JID(WAT,i) ) = J( JID(OIL,i), JID(WAT,i) ) - porosity/DT * dens(OIL,i);
    J( JID(WAT,i), JID(WAT,i) ) = J( JID(WAT,i), JID(WAT,i) ) + porosity/DT * dens(WAT,i);
end

%% producers with constant BHP
for well_id = 1 : length( wid_prd )
    cell_id = wid_prd( well_id );
    dQ_dP = WI * rperm( :, cell_id ) ./viscosity' .* ( dens_DP( :, cell_id ) * ( P(cell_id) - BHP_PRD(well_id) ) + dens( :, cell_id ) ) / DV;
    dQ_dS = WI * dens( :, cell_id ) ./viscosity' * ( P(cell_id) - BHP_PRD(well_id) ) .* rperm_DS( :, cell_id ) / DV;

    J( JID(OIL,cell_id) : JID(WAT,cell_id), JID(OIL,cell_id) ) = J( JID(OIL,cell_id) : JID(WAT,cell_id), JID(OIL,cell_id) ) + dQ_dP;
    J( JID(OIL,cell_id) : JID(WAT,cell_id), JID(WAT,cell_id) ) = J( JID(OIL,cell_id) : JID(WAT,cell_id), JID(WAT,cell_id) ) + dQ_dS;
end



%% injectors with constant Qw
for well_id = 1 : length( wid_inj ) 
    cell_id = wid_inj( well_id );
    bhp_id = Ncell + well_id;
    new_res_id  = 2 * Ncell + well_id;
    
    dQ_dP = WI * rperm( WAT, cell_id ) /viscosity(WAT) * ( dens_DP( WAT, cell_id ) * ( P(cell_id) - P(bhp_id) ) + dens( WAT, cell_id ) ) / DV;
    dQ_dS = WI * dens( WAT, cell_id ) /viscosity(WAT) * ( P(cell_id) - P(bhp_id) ) * rperm_DS( WAT, cell_id ) / DV;
    
    
    dQ_dBHP = -1 * WI * rperm( WAT, cell_id ) * dens( WAT, cell_id ) / viscosity(WAT) / DV;    

    J( new_res_id, JID( OIL, cell_id ) ) = dQ_dP;
    J( new_res_id, JID( WAT, cell_id ) ) = dQ_dS;
   
%     'dQ-dS'
%     p( P(cell_id) );
%     p( P(bhp_id) );
%     p( dQ_dS );
    
    J( new_res_id, new_res_id ) = dQ_dBHP;
end

%% Injectors with constant BHP - ONLY WATER
% for well_id = 1 : length(wid_inj_bhp)
%     cell_id = wid_inj_bhp(well_id);
%     dQ_dP = WI * rperm( :, cell_id ) ./viscosity' .* ( dens_DP( :, cell_id ) * ( P(cell_id) - BHP_INJ(well_id) ) + dens( :, cell_id ) ) / DV;
%     dQ_dS = WI * dens( :, cell_id ) ./viscosity' * ( P(cell_id) - BHP_INJ(well_id) ) .* rperm_DS( :, cell_id ) / DV;
%     
% %     print( dQ_dP(2) )
% %     print( dQ_dS(2) )
% 
%     J( JID(WAT,cell_id), JID(OIL,cell_id) ) = J( JID(WAT,cell_id), JID(OIL,cell_id) ) + dQ_dP(WAT);
%     J( JID(WAT,cell_id), JID(WAT,cell_id) ) = J( JID(WAT,cell_id), JID(WAT,cell_id) ) + dQ_dS(WAT);
% end