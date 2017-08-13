% flow term at interface in discretized domain
function [f] = ResFlow( interf_id, OIL, WAT, densf, rpermf, DP )

f( OIL, 1 ) = densf( OIL, interf_id ) * rpermf( OIL, interf_id ) * DP( interf_id );
f( WAT, 1 ) = densf( WAT, interf_id ) * rpermf( WAT, interf_id ) * DP( interf_id );
    