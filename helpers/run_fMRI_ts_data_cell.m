function M = run_fMRI_ts_data_cell(X, opts, measure)

if nargin<3, error('myApp:argChk', 'At least 3 arguments are required.'); end

M.ID = X.ID;
M.nodes = X.Nodes;
resize=false;
X2=X.TS;
idx=find(sum(X2,1)==0);

if ~isempty(idx)
    resize=true;
    warning_msg = ['ATTENTION: Your TS ID ', M.ID, ' contains column(s) with just 0. This column(s) is going to be removed in order to get safisfactory results. Note that its output is going to have -1 in the column(s) removed.'];
    disp(warning_msg);
    % warning
    for k=1:length(idx)
        idx2=find(sum(X2,1)==0);
        X2(:,idx2(1))=[];
    end
end

if ~isempty(find(X2==0)), error('TS cannot contain any zero.'); end

CM = connMap(X2, opts);
null_model = false;

%M.CM = CM;

if isfield(opts, 'null_model') && opts.null_model
    null_model = true;
    [W0,R] = null_model_und_sign_helper(CM,opts);
end

if ~isfield(measure,'name'), error('myApp:argChk', 'A measure (varargin{i}.name) is required.'); end
measure.name=lower(measure.name);

switch measure.name
    case 'betweenness'
        BC = betweenness_wei_helper(CM);
        if resize, BC = vector_resize(BC, idx); end
        M.BC=BC';
        
        if null_model,
            BC_null_model = betweenness_wei_helper(W0);
            if resize, BC_null_model = vector_resize(BC_null_model, idx); end
            M.BC_null_model=BC_null_model';
        end
        
    case 'degrees'
        [deg] = degrees_und(CM);
        if resize, deg = vector_resize(deg, idx); end
        M.deg=deg;
        
        if null_model
            deg_null_model = degrees_und(W0);
            if resize, deg_null_model = vector_resize(deg_null_model, idx); end
            M.deg_null_model=deg_null_model;
        end
        
    case 'strength'
        if isfield( measure, 'sign' ) && measure.sign,
            [Spos Sneg vpos vneg] = strengths_und_sign(CM);
            if resize
                Spos = vector_resize(Spos, idx);
                Sneg = vector_resize(SNeg, idx);
            end
            M.strength_Spos = Spos;
            M.strength_Sneg = Sneg;
        else,
            str = strengths_und(CM);
            if resize, str = vector_resize(str, idx); end
            M.strength = str;
            
        end;
        
        if null_model,
            if isfield( measure, 'sign' ) && measure.sign,
                [Spos_null_model Sneg_null_model vpos vneg] = strengths_und_sign(W0);
                if resize
                    Spos_null_model = vector_resize(Spos_null_model, idx);
                    Sneg_null_model = vector_resize(Sneg_null_model, idx);
                end
                M.strength_Spos_null_model = Spos_null_model;
                M.strength_Sneg_null_model = Sneg_null_model;
            else,
                str_null_model = strengths_und(W0);
                if resize, str_null_model = vector_resize(str_null_model, idx); end
                M.str_null_model=str_null_model;
            end;
        end
        
    case 'clustering'
        C=clustering_coef_wu(CM);
        if resize, C = vector_resize(C, idx); end
        M.C = C;
        
        if null_model
            C_null_model=clustering_coef_wu(W0);
            if resize, C_null_model = vector_resize(C_null_model, idx); end
            M.C_null_model=C_null_model;
        end
        
    case 'eigenvector'
        v = eigenvector_centrality_und(CM); % binary/weighted undirected adjacency matrix.
        if resize, v = vector_resize(v, idx); end
        M.v = v;
        
        if null_model
            v_null_model = eigenvector_centrality_und(W0);
            if resize, v_null_model = vector_resize(v_null_model, idx); end
            M.v_null_model=v_null_model;
        end
        
    case 'modularity'
        [Ci Q] = modularity_helper(CM, measure);
        if resize, Ci = vector_resize(Ci, idx); end
        M.Q = Q;
        M.Ci=Ci;
        
        if null_model
            [Ci_null_model Q_null_model] = modularity_helper(W0, measure);
            if resize, Ci_null_model = vector_resize(Ci_null_model, idx); end
            M.Ci_null_model=Ci_null_model;
            M.Q_null_model=Q_null_model;
        end
        
    case 'consensus'
        [ciu Q] = consensus_und_helper(CM, measure);
        if resize, ciu = vector_resize(ciu, idx); end
        M.ciu = ciu;
        M.Q=Q;
        if null_model
            [ciu_null_model Q_null_model] = consensus_und_helper(W0, measure);
            if resize, ciu_null_model = vector_resize(ciu_null_model, idx); end
            M.ciu_null_model=ciu_null_model;
            M.Q_null_model=Q_null_model;
        end
        
    case 'rich_club'
        rich = rich_club_helper(CM,measure);
        M.rich=rich;
        
        if null_model
            rich_null_model = rich_club_helper(W0,measure);
            M.Rich_null_model=rich_null_model;
        end
        
    case 'pagerank'
        pr = pagerank_centrality_helper(CM, measure); % adjacency matrix
        if resize, pr = vector_resize(pr, idx); end
        M.pagerank=pr;
        
        if null_model
            pr_null_model = pagerank_centrality_helper(W0, measure);
            if resize, pr_null_model = vector_resize(pr_null_model, idx); end
            M.pagerank_null_model=pr_null_model;
        end
        
    case 'diversity'
        [Hpos Hneg] = diversity_coef_sign_helper(CM, measure);
        if resize
            Hpos = vector_resize(Hpos, idx);
            Hneg = vector_resize(Hneg, idx);
        end
        M.Hpos=Hpos; M.Hneg=Hneg;
        
        if null_model
            [Hpos_null_model Hneg_null_model] = diversity_coef_sign_helper(W0, measure);
            if resize
                Hpos_null_model = vector_resize(Hpos_null_model, idx);
                Hneg_null_model = vector_resize(Hneg_null_model, idx);
            end
            M.Hpos_null_model=Hpos_null_model; M.Hneg_null_model=Hneg_null_model;
        end
        
    case 'participation'
        [Ppos Pneg] = participation_coef_sign_helper(CM, measure);
        if resize
            Ppos = vector_resize(Ppos, idx);
            Pneg = vector_resize(Pneg, idx);
        end
        M.Ppos=Ppos; M.Pneg=Pneg;
        
        if null_model
            [Ppos_null_model Pneg_null_model] = participation_coef_sign_helper(W0, measure);
            if resize
                Ppos_null_model = vector_resize(Ppos_null_model, idx);
                Pneg_null_model = vector_resize(Pneg_null_model, idx);
            end
            M.Ppos_null_model=Ppos_null_model;
            M.Pneg_null_model=Pneg_null_model;
        end
        
    otherwise
        error('Measure named %s is unknown.', measure.name);
end

end

function r = vector_resize(v, idx)
for i=1:length(idx)
    if size(v,1) > size(v,2), v=v'; end;
    v2 = v;
    v = [v2(1:idx(i)-1) -1 v2(idx(i):end)];
    r = v';
end
end