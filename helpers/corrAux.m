function p_nm = corrAux(x1,x2, opts)

% p_nm = corrAux(x1, x2, opts)
%
%   This function returns the correlation value between two filtered
%   vectors. This function calls partialcorr.m.
%
%   Input:     x1,        required - filtered vector for region n
%
%              x2,        required - filtered vector for region m
%
%              opts,      required - matlab structure
%                             opts.corrType   correlation type
%                                                   'pearson' (default) 
%                                                   'kendall'
%                                                   'spearman'
%                                                   'partial'
%
%                             opts.partialZ   correlation control
%                                             (required to calculate
%                                             partial correlation type)
%                                             
%
%   Output:    p_nm,      correlation value
%
%
%   Example:   x1 = rand(1,225); 
%              x2 = rand(1,225);
%              opts.corrType = 'kendall';
%              y = corrAux(x1, x2, opts);
% 

warning off;

if ~exist('x1','var'), error('myApp:argChk', 'Parameter "x1" (vector) is required.');
elseif ~exist('x2','var'), error('myApp:argChk', 'Parameter "x2" (vector) is required.');
elseif ~isnumeric(x1) || ~isnumeric(x2), error('myApp:argChk', 'Parameter "x1" and "x2" must be a valid vector.');
elseif isempty(x1) || isempty(x2), error('myApp:argChk', 'Parameter "x1" and "x2" cannot be empty.');
elseif ~(length(x1)==length(x2)), error('myApp:argChk', 'Parameter "x1" and "x2" must have the same size.');
end

if size(x1,1) > size(x1,2), x1=x1'; end;
if size(x2,1) > size(x2,2), x2=x2'; end;

X = [x1 ; x2];

if ~exist('opts', 'var'), error('myApp:argChk','Parameter "opts" (MatLab struct) is required.');
elseif ~isstruct(opts), error('myApp:argChk', 'Parameter "opts" must be a MatLab struct.');
elseif ~isfield(opts,'corrType'), opts.corrType = 'pearson';
elseif ~ischar(opts.corrType), error('myApp:argChk', 'Parameter "opts" must be a MatLab string.');
else opts.corrType = lower(opts.corrType);
end

switch opts.corrType
    case {'kendall','pearson','spearman'}
        R = corr(X', 'type', opts.corrType);
    case 'partial'
        if ~isfield(opts,'partialZ'), error('myApp:argChk', 'Parameter "opts.partialZ" is required to calculate Partial Correlation.'); 
        elseif ~isnumeric(opts.partialZ) || isempty(opts.partialZ), error('myApp:argChk', 'Parameter "opts.partialZ" must be a valid vector.');
        elseif ~size(X,1)==size(opts.partialZ,1), error('myApp:argChk', 'X and Z must have the same number of rows.');
        end
        R = partialcorr( X', opts.partialZ');
    otherwise
        error('myApp:argChk', 'Parameter "opts.corrType" is invalid.');
end

p_nm = R(1,2);

end