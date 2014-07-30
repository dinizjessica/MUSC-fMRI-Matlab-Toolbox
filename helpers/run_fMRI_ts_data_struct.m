function M = run_fMRI_ts_data_struct(X, opts, measure)

if nargin<3, error('myApp:argChk', 'All 3 arguments are required.'); end

if isempty(X), error('myApp:argChk', 'TS must not be empty.'); end

nodes=size(X,2);
resize=false;
idx=find(sum(X,1)==0);

if ~isempty(idx)
    resize=true;
    disp('ATTENTION: Your TS contains column(s) with just 0. This column(s) is going to be removed in order to get safisfactory results');
    % warning
    for k=1:length(idx)
        idx2=find(sum(X,1)==0);
        X(:,idx2(1))=[];
    end
end

if ~isempty(find(X==0)), error('TS cannot contain any zero.'); end

if ~isfield(opts, 'filename'), opts.filename='brain_network_measures.csv'; end

CM = connMap(X, opts);
null_model = false;
aux=1;
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
        
        betweenness.name='Betweenness';
        betweenness.values=BC';
        spreadsheet{aux}=betweenness;
        aux=aux+1;
        M.BC = BC;
        if null_model,
            BC_null_model = betweenness_wei_helper(W0);
            if resize, BC_null_model = vector_resize(BC_null_model, idx); end
            
            betweenness_null_model.name='Betweenness - Null Model';
            betweenness_null_model.values=BC_null_model';
            spreadsheet{aux}=betweenness_null_model;
            aux=aux+1;
            M.BC_null_model=BC_null_model;
        end
    case 'degrees'
        [deg] = degrees_und(CM);
        if resize, deg = vector_resize(deg, idx); end
        
        degrees.name='Degrees'; degrees.values=deg;
        spreadsheet{aux}=degrees;
        aux=aux+1;
        M.deg = deg;
        if null_model
            deg_null_model = degrees_und(W0); M.deg_null_model=deg_null_model;
            if resize, deg_null_model = vector_resize(deg_null_model, idx); end
            
            degrees_null_model.name='Degress - Null Model';
            degrees_null_model.values=deg_null_model;
            spreadsheet{aux}=degrees_null_model;
            aux=aux+1;
        end
    case 'strength'
        if isfield( measure, 'sign' ) && measure.sign,
            
            [Spos Sneg vpos vneg] = strengths_und_sign(CM);
            if resize
                Spos = vector_resize(Spos, idx);
                Sneg = vector_resize(Sneg, idx);
            end
            
            strength_Spos.name='Strength Positive';
            strength_Spos.values=Spos;
            spreadsheet{aux}=strength_Spos;
            aux=aux+1;
            strength_Sneg.name='Strength Negative';
            strength_Sneg.values=Sneg;
            spreadsheet{aux}=strength_Sneg;
            aux=aux+1;
            M.strength_Spos = strength_Spos.values;
            M.strength_Sneg = strength_Sneg.values;
            
        else,
            str = strengths_und(CM);
            if resize, str = vector_resize(str, idx); end
            
            strength.name='Strength';
            strength.values = str;
            spreadsheet{aux}=strength;
            aux=aux+1;
            M.str = str;
            
        end;
        
        if null_model,
            if isfield( measure, 'sign' ) && measure.sign,
                [Spos Sneg vpos vneg] = strengths_und_sign(W0);
                if resize
                    Spos = vector_resize(Spos, idx);
                    Sneg = vector_resize(Sneg, idx);
                end
                
                strength_Spos.name='Strength Positive - Null Model';
                strength_Spos.values=Spos;
                spreadsheet{aux}=strength_Spos;
                aux=aux+1;
                strength_Sneg.name='Strength Negative - Null Model';
                strength_Sneg.values=Sneg;
                spreadsheet{aux}=strength_Sneg;
                aux=aux+1;
                M.strength_null_model_Spos = strength_Spos.values;
                M.strength_null_model_Sneg = strength_Sneg.values;
            else,
                str_null_model = strengths_und(W0);
                if resize, str_null_model = vector_resize(str_null_model, idx); end
                
                strength_null_model.values = str_null_model;
                strength_null_model.name='Strength - Null Model';
                spreadsheet{aux}=strength_null_model;
                aux=aux+1;
                M.strength_null_model=str_null_model;
            end;
            
        end
    case 'clustering'
        C=clustering_coef_wu(CM);
        if resize, C = vector_resize(C, idx); end
        
        M.C = C;
        clustering.name='Clustering'; clustering.values=C;
        spreadsheet{aux}=clustering;
        aux=aux+1;
        if null_model
            C_null_model=clustering_coef_wu(W0);
            if resize, C_null_model = vector_resize(C_null_model, idx); end
            
            M.C_null_model=C_null_model;
            clustering_null_model.name='Clustering - Null Model';
            clustering_null_model.values=C_null_model;
            spreadsheet{aux}=clustering_null_model;
            aux=aux+1;
        end
    case 'eigenvector'
        v = eigenvector_centrality_und(CM); % binary/weighted undirected adjacency matrix.
        if resize, v = vector_resize(v, idx); end
        
        M.v = v;
        eigenvector.name='Eigenvector'; eigenvector.values=v;
        spreadsheet{aux}=eigenvector;
        aux=aux+1;
        if null_model
            v_null_model = eigenvector_centrality_und(W0);
            if resize, v_null_model = vector_resize(v_null_model, idx); end
            eigenvector_null_model.name='Eigenvector - Null Model'; eigenvector_null_model.values=v_null_model;
            spreadsheet{aux}=eigenvector_null_model;
            aux=aux+1;
            M.v_null_model=v_null_model;
        end
    case 'modularity'
        [Ci Q] = modularity_helper(CM, measure);
        if resize, Ci = vector_resize(Ci, idx); end
        M.Ci = Ci; M.modularity_Q = Q;
        modularity.name='Modularity'; modularity.values=Ci;
        spreadsheet{aux}=modularity;
        aux=aux+1;

        z=zeros(1,nodes-1);
        Q= [Q z];
        modularity_q.name='Q - Modularity'; modularity_q.values=Q;
        spreadsheet{aux}=modularity_q;
        aux=aux+1;
        
        if null_model
            [Ci_null_model Q_null_model] = modularity_helper(W0, measure);
            if resize
                Ci_null_model = vector_resize(Ci_null_model, idx);
            end
            modularity_null_model.name='Modularity - Null Model';
            modularity_null_model.values=Ci_null_model;
            spreadsheet{aux}=modularity_null_model;
            aux=aux+1;
            
            Q_null_model = [Q_null_model z];
            modularity_q_null_model.name='Q - Modularity - Null Model'; 
            modularity_q_null_model.values=Q_null_model;
            spreadsheet{aux}=modularity_q_null_model;
            aux=aux+1;
            M.Ci_null_model=Ci_null_model; M.modularity_Q_null_model=modularity_q_null_model;
        end
    case 'consensus'
        [ciu Q] = consensus_und_helper(CM, measure);
        if resize, ciu = vector_resize(ciu, idx); end
        
        M.ciu = ciu;
        consensus.name='Consensus'; consensus.values=ciu;
        spreadsheet{aux}=consensus;
        aux=aux+1;
        
        z=zeros(1,nodes-1);
        Q= [Q z];
        consensus_q.name='Q - Consensus'; consensus_q.values=Q;
        spreadsheet{aux}=consensus_q;
        aux=aux+1;
        M.consensus_Q = Q;
        
        if null_model
            [ciu_null_model Q_null_model]= consensus_und_helper(W0, measure);
            if resize, ciu_null_model = vector_resize(ciu_null_model, idx); end
            
            M.ciu_null_model=ciu_null_model;
            consensus_null_model.name='Consensus - Null Model';
            consensus_null_model.values=ciu_null_model;
            spreadsheet{aux}=consensus_null_model;
            aux=aux+1;
            
            Q_null_model = [Q_null_model z];
            consensus_q_null_model.name='Q - Consensus - Null Model'; 
            consensus_q_null_model.values=Q_null_model;
            spreadsheet{aux}=consensus_q_null_model;
            aux=aux+1;
        end
    case 'rich_club'
        r = rich_club_helper(CM,measure);
        if length(r)<nodes
           l=nodes-length(r); 
           if resize, l=l+length(idx); end
           z=zeros(1,l);
           r= [r z];
        end
        rich.values = r;
        rich.name = 'Rich Club';
        spreadsheet{aux}=rich;
        aux=aux+1;
        M.rich=r;
        
        if null_model
            r_null_model = rich_club_helper(W0,measure);
            if length(r_null_model)<nodes
                l=nodes-length(r_null_model);
                if resize, l=l+length(idx); end
                z=zeros(1,l);
                r_null_model=[r_null_model z];
            end
            rich_null_model.values = r_null_model;
            rich_null_model.name = 'Rich Club - Null Model';
            spreadsheet{aux}=rich_null_model;
            aux=aux+1;
            M.Rich_null_model=r_null_model;
        end
    case 'pagerank'
        pr = pagerank_centrality_helper(CM, measure); % adjacency matrix
        if resize, pr = vector_resize(pr, idx); end
        
        pagerank.values = pr;
        pagerank.name='Pagerank Centrality';
        spreadsheet{aux}=pagerank;
        aux=aux+1;
        M.pagerank=pr;
        
        if null_model
            pr_null_model = pagerank_centrality_helper(W0, measure);
            if resize, pr_null_model = vector_resize(pr_null_model, idx); end
            
            pagerank_null_model.values = pr_null_model;
            pagerank_null_model.name = 'Pagerank Centrality - Null Model';
            spreadsheet{aux}=pagerank_null_model;
            aux=aux+1;
            M.pagerank_null_model=pr_null_model;
        end
    case 'diversity'
        [Hpos Hneg] = diversity_coef_sign_helper(CM, measure);
        if resize
            Hpos = vector_resize(Hpos, idx);
            Hneg = vector_resize(Hneg, idx);
        end
        
        diversity_Hpos.values = Hpos;
        diversity_Hpos.name = 'Diversity Positive';
        spreadsheet{aux}=diversity_Hpos;
        aux=aux+1;
        diversity_Hneg.values = Hneg;
        diversity_Hneg.name = 'Diversity Negative';
        spreadsheet{aux}=diversity_Hneg;
        aux=aux+1;
        M.Hpos=Hpos; M.Hneg=Hneg;
        
        if null_model
            [Hpos_null_model Hneg_null_model] = diversity_coef_sign_helper(W0, measure);
            if resize
                Hpos_null_model = vector_resize(Hpos_null_model, idx);
                Hneg_null_model = vector_resize(Hneg_null_model, idx);
            end
            
            diversity_Hpos_null_model.values = Hpos_null_model;
            diversity_Hpos_null_model.name = 'Diversity Positive - Null Model';
            spreadsheet{aux}=diversity_Hpos_null_model;
            aux=aux+1;
            diversity_Hneg_null_model.values = Hneg_null_model;
            diversity_Hneg_null_model.name = 'Diversity Negative - Null Model';
            spreadsheet{aux}=diversity_Hneg_null_model;
            aux=aux+1;
            M.Hpos_null_model=Hpos_null_model; M.Hneg_null_model=Hneg_null_model;
        end
        
    case 'participation'
        [Ppos Pneg] = participation_coef_sign_helper(CM, measure);
        if resize
            Ppos = vector_resize(Ppos, idx);
            Pneg = vector_resize(Pneg, idx);
        end
        
        participation_Ppos.values = Ppos;
        participation_Ppos.name = 'Participation Positive';
        spreadsheet{aux}=participation_Ppos;
        aux=aux+1;
        participation_Pneg.values = Pneg;
        participation_Pneg.name = 'Participation Negative';
        spreadsheet{aux}=participation_Pneg;
        aux=aux+1;
        M.Ppos=Ppos; M.Pneg=Pneg;
        
        if null_model
            [Ppos_null_model Pneg_null_model] = participation_coef_sign_helper(W0, measure);
            if resize
                Ppos_null_model = vector_resize(Ppos_null_model, idx);
                Pneg_null_model = vector_resize(Pneg_null_model, idx);
            end
            
            participation_Ppos_null_model.values = Ppos_null_model;
            participation_Ppos_null_model.name = 'Participation Positive - Null Model';
            spreadsheet{aux}=participation_Ppos_null_model;
            aux=aux+1;
            participation_Pneg_null_model.values = Pneg_null_model;
            participation_Pneg_null_model.name = 'Participation Negative - Null Model';
            spreadsheet{aux}=participation_Pneg_null_model;
            aux=aux+1;
            M.Ppos_null_model=Ppos_null_model; M.Pneg_null_model=Pneg_null_model;
        end
    otherwise
        error('Measure named %s is unknown.', measure.name);
end


M.spreadsheet=spreadsheet;

end


function r = vector_resize(v, idx)
for i=1:length(idx)
    if size(v,1) > size(v,2), v=v'; end;
    v2 = v;
    v = [v2(1:idx(i)-1) -1 v2(idx(i):end)];
    r = v';
end
end