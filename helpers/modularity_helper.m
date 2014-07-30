function [Ci Q] = modularity_helper(CM, type)
    
% CM = connection matrix
% type = struct / especifies type of modularity
%
% type.name = char ('modularity')
% type.modularity_type = char
% type.sign = logical
% type.qtype = char / modularity type
% type.Ci0 = numeric vector / initial community affiliation vector
% type.gamma = numeric / modularity resolution parameter
% type.p = numeric / probability of random node moves / positive or
% negative
%                       
%
% Ci = refined community affiliation vector
% Q = modularity (qtype dependent)


validateattributes(CM,{'numeric'},{'square','nonempty','nonzero'});
validateattributes(type,{'struct'},{'nonempty'});

qtype = 'sta';                  % default values defined by 
Ci0 = 1:length(CM);             % all modularity functions
gamma = 1;

if isfield(type,'qtype') && ischar(type.qtype), qtype=type.qtype; end

if isfield(type,'Ci0') 
    if isvector(type.Ci0) && length(type.Ci0)==length(CM), Ci0 = type.Ci0;
    else error('myApp:argChk', 'Ci0 must be a vector and must have the same length as CM.'); end
end

if isfield(type,'gamma')
    if isnumeric(type.gamma) && type.gamma>=0, gamma = type.gamma; 
    else error('myApp:argChk', 'Gamma must be a number and greater or equal to 0.'); end
end

if ~isfield(type,'modularity_type'), type.modularity_type = 'modularity_und'; 
else type.modularity_type=lower(type.modularity_type); end


switch type.modularity_type
    case 'finetune'
        if isfield(type, 'sign') && type.sign
            [Ci Q] = modularity_finetune_und_sign(CM,qtype,Ci0);
        else
            [Ci Q] = modularity_finetune_und(CM,Ci0,gamma);
        end
        
    case 'louvain'
        if isfield(type, 'sign') && type.sign            
            [Ci Q] = modularity_louvain_und_sign(CM,qtype);
        else
            [Ci Q] = modularity_louvain_und(CM,gamma);
        end
        
    case 'probtune'
        if ~isfield(type,'p') || ~isnumeric(type.p), error('myApp:argChk', 'Probability of random node moves (type.p) is required.'); end
        
        [Ci Q] = modularity_probtune_und_sign(CM,qtype,Ci0,type.p);
        
    case 'modularity_und'
        [Ci Q] = modularity_und(CM,gamma);
        
    otherwise
        error('modularity_type unknown');
        
end

end