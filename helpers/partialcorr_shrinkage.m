function [out,lc,lv] = partialcorr_shrinkage(x, ztrans)

% this function is made to estimate the partial correlation, calling 
% shrinkage if neccessary with shrinkage
% modified by Davy Vanderweyen

if ismatrix(x)==0
    error(message('stats:partialcorr: Inputs Must Be Matrices'));
end

if nargin ==1
    ztrans = 0;
end
[~,d] = size(x);

R = rank(cov(x));
if R >= d % do simple inversion
    disp('Simple inversion')
    ic=-inv(cov(x));
    r=(ic ./ repmat(sqrt(diag(ic)),1,d)) ./ repmat(sqrt(diag(ic))',d,1);
    r=r - diag(diag(r));
else % we need to adjust covariance matrix
    disp('Using shrinkage')
    [CV,lc,lv] = covshrinkage(x,1);
%     CV = shrinkDiag(x, .02);
    ic = -inv(CV);
    r=(ic ./ repmat(sqrt(diag(ic)),1,d)) ./ repmat(sqrt(diag(ic))',d,1);
    r=r - diag(diag(r));
end
[m,~]=size(r);
id=eye(m);
corrmat= r+id;

if ztrans
    out = 0.5 * log((1+r) ./ (1-r));
else
    out =corrmat;
end
return
