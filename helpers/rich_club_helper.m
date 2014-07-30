function R = rich_club_helper(CM, type)

% type = struct
% type.name
% type.binary = logical
% type.klevel = numeric ---- can it be negative or zero?
% threshold_type optional
% threshold.thr required if threshold_type='absolute'
% threshold.threshold_p = numeric / proportion of weights to preserve / default 1


W = threshold_helper(CM, type);

if ~isfield(type,'binary'), type.binary=false;
elseif ~islogical(type.binary), error('myApp:argChk', 'Varargin.binary must be logical.'); end

if type.binary
    
    if ~isfield(type,'klevel'), type.klevel = max(sum(W));
    elseif ~isnumeric(type.klevel), error('myApp:argChk', 'Varargin.klevel must be a number.'); end

    WC = weight_conversion(W, 'binarize');
    
    [R,Nk,Ek] = rich_club_bu(WC,type.klevel);
    
else
    
    if ~isfield(type,'klevel'), type.klevel = max(sum(W~=0));
    elseif ~isnumeric(type.klevel), error('myApp:argChk', 'Varargin.klevel must be a number.'); end

    R = rich_club_wu(W,type.klevel);
end

end