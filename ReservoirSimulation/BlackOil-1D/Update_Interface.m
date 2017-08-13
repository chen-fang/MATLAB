function [DP,densf, rpermf] = Update_Interface( Nface, OIL, WAT, P, dens, rperm )

for i = 1 : Nface
    DP(i) = P(i+1) - P(i);
    densf( OIL,i ) = UpstreamNTP( dens(OIL,i), dens(OIL,i+1), DP(i) );
    densf( WAT,i ) = UpstreamNTP( dens(WAT,i), dens(WAT,i+1), DP(i) );
    
    rpermf( OIL,i ) = UpstreamNTP( rperm(OIL,i), rperm(OIL,i+1), DP(i) );
    rpermf( WAT,i ) = UpstreamNTP( rperm(WAT,i), rperm(WAT,i+1), DP(i) );
end