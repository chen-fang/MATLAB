% partial_to_left:  D[ f(i+1/2) ] / D[ x(i) ]
% partial_to_right: D[ f(i+1/2) ] / D[ x(i+1) ]

% PTN - Point To Neighbor
% Start from one point, determine the values at its neighbors with upstream
% policy
function [f] = UpstreamPTN( L, R, DP, OIL, oil_to_left, oil_to_right, WAT, wat_to_left, wat_to_right )
% 1 - left
% 2 - right
f = zeros( 2, 2 );
if DP >= 0
    f( OIL,R ) = oil_to_right;
    f( WAT,R ) = wat_to_right;
else
    f( OIL,L ) = oil_to_left;
    f( WAT,L ) = wat_to_left;
end