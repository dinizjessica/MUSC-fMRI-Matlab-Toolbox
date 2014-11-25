function M = connMap(X, opts)

% M = connMap(X, opts)
%
% This function calculates the correlation value of the given matrix X.
% This function calls filter_helper.m and corrAux.m.
%
% Input:    X,        Required. Matrix of values
%           
%           opts,     Required. Matlab structure.
%                               opts.f_l        low pass frequency cutoff
%                                               in Hz - required
%                                                   (0 < opts.f_l < 1)
%
%                               opts.f_h        high pass frequency cutoff
%                                               in Hz - required
%                                                   (0 < opts.f_h < 1)
%
%                               opts.f_s        sampling frequency in Hz -
%                                               optional
%                                                   (opts.f_s > 0)
%
%                               opts.n          number of FIR coefficients
%                                               - optional
%                                                   (default 10)
%
%                               opts.type       filter type - required
%                                                   'fir' (default)
%                                                   'butter'
%
%                               opts.corrType   correlation type
%                                                   'pearson' (default) 
%                                                   'kendall'
%                                                   'spearman'
%                                                   'partial'
%
%
% Output:   M,     returns a matrix with the correlation value. The matrix
%                  M has dimension n (where n is equal to the number of
%                  nodes of the matrix M. In other words, the number of
%                  rows of the given matrix).
%
%
% Example: 
%         X = rand(225,225); 
%         opts.f_l = 0.01;      % required
%         opts.f_h = 0.1;       % required
%         opts.type = 'butter';
%         opts.corrType = 'kendall';
%         y = connMap(X, opts);
% 

if ~exist('X','var'), error('myApp:argChk', 'Parameter "X" (matrix) is required.');
elseif ~isnumeric(X), error('myApp:argChk', 'Parameter "X" must be a numeric matrix.');
elseif  isempty(X), error('myApp:argChk', 'Parameter "X" must be a valid matrix.');
end

if ~exist('opts', 'var'), error('myApp:argChk','Parameter "opts" (MatLab struct) is required.');
elseif ~isstruct(opts), error('myApp:argChk', 'Parameter "opts" must be a MatLab struct.');
end

M = zeros(size(X,2));

if isfield(opts,'type') && ~strcmpi( opts.type,'none')

    for y = 1:size(X,2)
        Xf(:,y) = filter_helper(X(:,y),opts);
    end
else
    Xf = X;
end;

if isfield(opts,'corrType'), 
    switch opts.corrType
        case 'partial'
            M=partialcorr_shrinkage(Xf);           % supports shrinkage
        case {'kendall','pearson','spearman'}
            M=corr(Xf, 'type', opts.corrType);
        
    end;
end;

end

%imshow(M); colormap(jet);
