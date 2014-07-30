function [kden,N,K] = density_und_helper(ConnMatrix, threshold)

% threshold = struct
% threshold.threshold_type = 'absolute' or 'proportional'
% threshold.thr = numeric / weight threshold
% threshold.p = numeric / proportion of weights to preserve

TM = threshold_helper(ConnMatrix,threshold);
[kden,N,K] = density_und(TM);

end