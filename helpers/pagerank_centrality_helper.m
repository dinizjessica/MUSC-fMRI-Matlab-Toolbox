function r = pagerank_centrality_helper(A, type) 

% A,               adjacency matrix
% type.d,          damping factor [REQUIRED] --- Can it be negative or zero as well?
% type.falff,      initial page rank probability (non-negative)

%if nargin<2, error('myApp:argChk', 'At least 2 arguments are required (A and d).'); end

if ~isfield(type,'d'), error('myApp:argChk', 'Damping factor (varargin.d) is required.');
elseif ~isnumeric(type.d), error('myApp:argChk', 'Damping factor (varargin.d) must be a number.'); end

if ~isfield(type,'falff'), 
    N = size(A,1);
    type.falff = ones(N,1)/N;
elseif ~isvector(type.falff), error('myApp:argChk', 'Initial page rank probability (varargin.falff) must be a vector.');
elseif length(find(type.falff<0))>0, error('myApp:argChk', 'Initial page rank probability (varargin.falff) must be non-negative.');
else
    r = pagerank_centrality(A, type.d, type.falff);
    return;
end
    
r = pagerank_centrality(A, type.d);

end