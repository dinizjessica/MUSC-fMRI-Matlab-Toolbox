function [Ppos Pneg]=participation_coef_sign_helper(CM,type)

% type = struct
% type.name 
% type.Ci = vector / community affiliation vector


if ~isfield(type,'Ci') 
    type.Ci = consensus_und_helper(CM,type);
elseif ~isvector(type.Ci)
    error('myApp:argChk', 'Community affiliation vector (type.Ci) must be a vector.');
end

[Ppos Pneg]=participation_coef_sign(CM,type.Ci);

end