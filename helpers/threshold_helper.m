function W = threshold_helper(M, threshold)

% threshold = struct
% threshold.threshold_type = 'absolute' or 'proportional' / REQUIRED
%
% threshold.thr = numeric / weight threshold / required if threshold_type='absolute'
% threshold.threshold_p = numeric / proportion of weights to preserve / required if threshold_type='proportional'

validateattributes(M,{'numeric'},{'square','nonempty','nonzero'});
validateattributes(threshold,{'struct'},{'nonempty'});

if ~isfield(threshold,'threshold_type'), error('myApp:argChk', 'threshold_type is required.');
else threshold.threshold_type=lower(threshold.threshold_type);
end

switch threshold.threshold_type
    case 'absolute'
        if ~isfield(threshold,'thr'), error('myApp:argChk', 'Weight threshold (thr) is required.');
        elseif ~isnumeric(threshold.thr), error('myApp:argChk', 'Weight threshold (thr) must be numeric.');
        end
        W = threshold_absolute(M, threshold.thr);
        
    case 'proportional'
        if ~isfield(threshold,'threshold_p'),
            error('myApp:argChk', 'Proportion of weights to preserve (threshold_p) is required.');
        elseif ~isnumeric(threshold.threshold_p), 
            error('myApp:argChk', 'Proportion of weights to preserve (threshold_p) must be numeric.');
        elseif threshold.threshold_p>1 || threshold.threshold_p<0, 
            error('myApp:argChk', 'Proportion of weights to preserve (threshold_p) must be in the range of 0<threshold_p<1.'); 
        end
        W = threshold_proportional(M, threshold.threshold_p);
    otherwise
        error('threshold_type unknown');
end
end