% NTP - Neighbor To Point
% Determine the value at one point according to the values at its neighbors with upstream
% policy

function [f] = UpstreamNTP( left, right, DP ) % DP = P(i+1) - P(i)

if DP >= 0
    f = right;
else
    f = left;
end