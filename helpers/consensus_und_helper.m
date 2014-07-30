function [ciu Q] = consensus_und_helper(ConnMatrix,type)

% type.buffsz = numeric / 'agreement' / optional second argument to set buffer size
% type.tau = numeric / threshold which controls the resolution of the reclustering
% type.reps = positive number / number of times that the clustering algorithm is reapplied

validateattributes(ConnMatrix,{'numeric'},{'square','nonempty','nonzero'});
validateattributes(type,{'struct'},{'nonempty'});

buffsz = 1000; % default

if isfield(type,'buffsz') && isnumeric(type.buffsz), buffsz=type.buffsz; end

if ~isfield(type,'tau'), error('myApp:argChk', 'Threshold which controls the resolution of the reclustering (type.tau) is required.'); 
elseif ~isnumeric(type.tau) || type.tau<0 || type.tau>1, error('myApp:argChk', 'Threshold which controls the resolution of the reclustering (type.tau) must be in the range of: 0 < tau <= 1.'); 
end

if ~isfield(type,'reps'), error('myApp:argChk', 'Number of times that the clustering algorithm is reapplied (type.reps) is required.');
elseif ~isnumeric(type.reps) || type.reps<1, error('myApp:argChk', 'Number of times that the clustering algorithm is reapplied (type.reps) must be a number greater that 0.');
end

Ci = modularity_louvain_und_sign(ConnMatrix);
D = agreement(Ci,buffsz);
[ciu Q] = consensus_und(D,type.tau,type.reps);

end