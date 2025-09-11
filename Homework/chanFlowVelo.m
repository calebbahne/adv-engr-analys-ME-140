function U = chanFlowVelo(n,S,B,H)
%CHANFLOWVELO: HW 3 P03.009 - 4 - calculates the channel flow velo
% see HW3_V2

U = (sqrt(S) ./ n) .* ((B .* H) ./ (B + 2 * H)).^(2/3);
end