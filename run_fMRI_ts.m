function M = run_fMRI_ts(X, opts, varargin)


%   M = run_fMRI_ts(X, opts, varargin);
%   
% 
%   This function generates a spreadsheet containing all results from the
%   selected measures for each node specified by the timeseries X.
%   This function can calculate node betweenness centrality, degree,
%   clustering coefficient, eigenector centrality, modularity, consensus
%   clustering, strength, rich club coefficient, pageRank centrality,
%   Shannon-entropy based diversity coefficient, and participation
%   coefficient. It is also possible to obtain the null model outcome to 
%   all measures. Some measures require additional parameters.
%
%   Note: Only rich_club outcome is not node related.
%
% Inputs:    X,             Required.
%                           X can either be a Matlab matrix or a cell.
%                           
%                           To run a single subject, X should be the
%                           corresponding timeserie matrix. To run
%                           multiple subjects, X must be a Matlab cell
%                           containing a Matlab structure to each subject.
%                           In this case, each structure of the cell
%                           must have ID (unique id that identifies the
%                           subject), TS (timeserie matrix), and Nodes
%                           (number of timeseries nodes) as a field.
%                            
%
%            opts,          Required.
%                           Matlab structure containing all possible
%                           options to filter the timeserie data and to 
%                           define if a null model outcome to all selected
%                           measures. It is also an option not to filt the
%                           data (this option should be chosen in case of
%                           the timeserie be already filtered) and/or not
%                           to get null model outcomes.
%
%                               opts.type       filter type (if ommitted -
%                                               no filtering is performed)
%                                                   'none' (default - no filtering performed)
%                                                   'fir' 
%                                                   'butter'
%
%                               opts.f_l        low pass frequency cutoff
%                                               in Hz - required to filter
%                                               the data
%                                                   (0 < opts.f_l < 1)
%
%                               opts.f_h        high pass frequency cutoff
%                                               in Hz - required to filter
%                                               the data
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
%                               opts.corrType   correlation type - optional
%                                                   'pearson' (default) 
%                                                   'kendall'
%                                                   'spearman'
%                                                   'partial'
%
%                               opts.null_model logical - optional
%                                                   (default false)
%
%                               opts.bin_swaps  Average number of swaps of each edge in binary randomization. 
%                                               Used to calculate null model - optional
%                                                    bin_swaps=5 is the default (each edge rewired 5 times)
%                                                    bin_swaps=0 implies no binary randomization
%
%                               opts.wei_freq   Frequency of weight sorting in weighted randomization.
%                                               Used to calculate null model - optional
%                                                    (0 < wei_freq <= 1)
%                                                    wei_freq=1 implies that weights are resorted at each step (default)
%                                                    wei_freq=0.1 implies that weights are sorted at each 10th step
%
%
%            varargin,      Matlab structure responsable for determining
%                           which measure must be calculate?and their
%                           respective parameters.
%                               name            selected measure. Required
%                                                   'betweenness'
%                                                   'strength'
%                                                   'degrees'
%                                                   'clustering'
%                                                   'eigenvector'
%                                                   'modularity'
%                                                   'consensus'
%                                                   'rich_club'
%                                                   'pagerank'
%                                                   'diversity'
%                                                   'participation'
%
%                               modularity_type     specifies modularity
%                                                   type - optional
%                                                       'modularity_und' (default)
%                                                           varargin{i}.gamma (optional)
%
%                                                       'louvain'
%                                                           varargin{i}.sign (optional)
%                                                                true
%                                                                   varargin{i}.qtype (optional)
%                                                                false 
%                                                                   varargin{i}.gamma (optional)
%
%                                                       'probtune'
%                                                           varargin{i}.p is required
%                                                           varargin{i}.Ci0 (optional)
%                                                           varargin{i}.qtype (optional)
%
%                                                       'finetune'
%                                                           varargin{i}.sign (optional)
%                                                                true
%                                                                   varargin{i}.qtype (optional)
%                                                                   varargin{i}.Ci0 (optional)
%                                                                false 
%                                                                   varargin{i}.gamma (optional)
%                                                                   varargin{i}.Ci0 (optional)
%
%                               sign                specifies if modularity_type is sign.
%                                                   Logical - optional
%                                                        (Default false)
%
%                               qtype               modularity type. (see Rubinov and Sporns, 2011)
%                                                       'sta',  Q_* (default)
%                                                       'pos',  Q_+
%                                                       'smp',  Q_simple
%                                                       'gja',  Q_GJA
%                                                       'neg',  Q_-
%
%                               Ci0                 initial community affiliation vector - optional
%
%                               gamma               modularity resolution parameter - optional
%                                                       gamma>1     detects smaller modules
%                                                       0<=gamma<1  detects larger modules
%                                                       gamma=1     (default) leads to the 'classic' modularity function
%
%                               p                   probability of random
%                                                   node moves. Required to
%                                                   calculate probtune.
%
%                               buffsz              optional second argument to
%                                                   set buffer size. Used
%                                                   by agreement.m before
%                                                   calculate consensus -
%                                                   optional
%                                                       (default 1000)
%
%                               tau                 threshold which controls the 
%                                                   resolution of the reclustering.
%                                                   Required to calculate consensus.
%
%                               reps                number of times that the clustering
%                                                   algorithm is reapplied.
%                                                   Required to calculate consensus.
%
%                               Ci                  community affiliation vector. Used to
%                                                   calculate diversity and participation 
%                                                   - optional
% 
%                               threshold_type      specifies threshold type. Used to threshold
%                                                   the connectivity matrix before calculate 
%                                                   rich_club - optional
%                                                       'absolute'
%                                                       'proportional' (default) 
%
%                               thr                 weight threshold. Used to threshold
%                                                   the connectivity matrix before calculate 
%                                                   rich_club - Required if
%                                                   threshold_type='absolute' 
%
%                               threshold_p         proportion of weights to preserve. Used to
%                                                   threshold the connectivity matrix before 
%                                                   calculate rich_club - optional
%                                                       (default 1)
%                                                               
%
% Output:    M,             Matlab cell containing all results from the 
%                           selected measures. If the input X is a Matlab
%                           matrix, then the output M is a cell containing
%                           a Matlab structure to each measure selected.
%                           However, if the input X is a Matlab cell, then
%                           the output M is a cell containing another cell
%                           to each measure selected. This measure cell 
%                           contains a Matlab structure with the 
%                           corresponding outcome to each subject of the
%                           data.
%                           
%
%            .csv file,     spreadsheet containing outcomes from all chosen
%                           measures. If the input X is a Matlab matrix,
%                           then the spreadsheet is going to have the
%                           results of all chosen measures in a single 
%                           file. In this case, the name of the file is 
%                           "fMRI_spreadsheet_mm-dd-yyyy_HH-MM-SS.csv"
%                           (mm-dd-yyyy is the current date, and HH-MM-SS
%                           the current time). If the input X is a Matlab
%                           cell, a new file is going to be created to each
%                           measure selected. In this case, the name of the file 
%                           is "measure_spreadsheet_mm-dd-yyyy_HH-MM-SS.csvs"
%                           (measure is the chosen measure, mm-dd-yyyy is
%                           the current date, and HH-MM-SS the current time)
%
%
% NULL MODEL
%
%   This function randomizes an undirected network with positive and
%   negative weights, while preserving the degree and strength
%   distributions. If null model outcome is expected, the chosen measure is
%   going to be calculate using a null model output as its input. 
%   This function calls null_model_und_sign.m (BCT toolbox).
%
%   Input:      object.name,        [required]
%               opts.null_model     [required]
%               opts.bin_swap       [optional]
%               opts.wei_freq       [optional]
%
%   Output:     null model output of the chosen measure.
%
%   Example:    opts.type = 'none';
%               y.name = 'betweenness';
%               M = run_fMRI_ts(X, opts, y); 
%               
%
% BETWEENNESS CENTRALITY
%
%   Node betweenness centrality is the fraction of all shortest paths in 
%   the network that contain a given node. Nodes with high values of 
%   betweenness centrality participate in a large number of shortest paths.
%   This function calls weight_conversion.m and betweenness_wei.m (BCT toolbox).
%
%   Input:      object.name,        [required]
%
%   Output:     M.BC,         node betweenness centrality vector.
%
%   Example:    opts.f_l = 0.1;
%               opts.f_h = 0.2;
%               opts.type = 'fir';
%               a.name = 'betweenness';
%               M = run_fMRI_ts(X, opts, a); 
% 
%
% NODE DEGREES 
%
%   Node degree is the number of links connected to the node. This
%   function calls degrees_und.m (BCT toolbox).
%   
%   Input:      object.name,        [required]
%
%   Output:     M.deg,       node degree.
%
%   Example:    opts.f_l = 0.1;
%               opts.f_h = 0.2;
%               b.name = 'degrees';
%               M = run_fMRI_ts(X, opts, b);
%
%
% CLUSTERING COEFFICIENT
%
%   The weighted clustering coefficient is the average "intensity" of 
%   triangles around a node. This function calls clustering_coef_wu.m
%   (BCT toolbox).
%
%   Input:      object.name,        [required]
%
%   Output:     M.C,        clustering coefficient vector.
%
%   Example:    opts.f_l = 0.1;
%               opts.f_h = 0.2;
%               opts.type='fir'
%               c.name = 'clustering';
%               M = run_fMRI_ts(X, opts, c);
%
%
% EIGENVECTOR CENTRALITY        Spectral measure of centrality
%
%   Eigenector centrality is a self-referential measure of centrality:
%   nodes have high eigenvector centrality if they connect to other nodes
%   that have high eigenvector centrality. The eigenvector centrality of
%   node i is equivalent to the ith element in the eigenvector 
%   corresponding to the largest eigenvalue of the adjacency matrix.
%   eigenvector. This function calls eigenvector_centrality_und.m (BCT
%   toolbox)
%
%   Input:      object.name,        [required]
%
%   Output:     M.v,        eigenvector associated with the largest
%                           eigenvalue of the adjacency matrix CIJ.
%
%   Example:    opts.f_l = 0.1;
%               opts.f_h = 0.2;
%               opts.type = 'fir';
%               d.name = 'eigenvector';
%               M = run_fMRI_ts(X, opts, d);
%
%
% CONSENSUS CLUSTERING
%
%   Consensus clustering seeks a consensus partition of the 
%   agreement matrix D. The algorithm used here is almost identical to the
%   one introduced in Lancichinetti & Fortunato (2012): The agreement
%   matrix D is thresholded at a level TAU to remove an weak elements. The
%   resulting matrix is then partitions REPS number of times using the
%   Louvain algorithm (in principle, any clustering algorithm that can
%   handle weighted matrixes is a suitable alternative to the Louvain
%   algorithm and can be substituted in its place). This clustering
%   produces a set of partitions from which a new agreement is built. If
%   the partitions have not converged to a single representative partition,
%   the above process repeats itself, starting with the newly built
%   agreement matrix. This function calls modularity_louvain_und_sign.m,
%   agreement.m, and consensus_und.m (BCT toolbox).
%   Note: It is not necessary run modularity before consensus. This script
%   is already running a modularity to calculate consensus. Even if
%   "modularity" and "consensus" are chosen as measures to run run_fMRI_ts,
%   the modularity output is not going to be consensus input, a new
%   modularity is calculated).
%
%   Input:      object.name,         [required]
%               object.buffsz,       [optional]
%               object.tau,          [required]
%               object.reps,         [required]
%
%   Output:     M.ciu,      consensus partition
%               M.Q,        maximized modularity
%
%   Example:    opts.null_model = false;
%               e.name = 'consensus';
%               e.buffsz = 1000;
%               e.tau = 0.5;    %required
%               e.reps = 5;     %required
%               M = run_fMRI_ts(X, opts, e);
%
%
% MODULARITY
%
%   The optimal community structure is a subdivision of the network into
%   nonoverlapping groups of nodes in a way that maximizes the number of
%   within-group edges, and minimizes the number of between-group edges.
%   The modularity is a statistic that quantifies the degree to which the
%   network may be subdivided into such clearly delineated groups. This
%   function calls modularity_und.m, modularity_finetune_und_sign.m,
%   modularity_finetune_und.m, modularity_louvain_und_sign.m,
%   modularity_louvain_und.m, or modularity_probtune_und_sign.m (BCT 
%   toolbox) according to the chosen modularity type.
%
%   Input:      object.name,                [required]
%               object.modularity_type,     [optional]
%               object.sign                 [optional]
%               object.qtype                [optional]
%               object.Ci0                  [optional]
%               object.gamma                [optional]
%               object.p                    [required to calculate probtune]
%
%   Output:     M.Ci,       optimal community structure
%               M.Q,        maximized modularity
%
%   Example:    opts.null_model = true;
%               f.name = 'modularity';
%               f.modularity_type = 'finetune';
%               f.sign = false;
%               f.gamma = 1;
%               M = run_fMRI_ts(X, opts, f);
%
%
% STRENGTHS
%
%   Node strength is the sum of weights of links connected to the node.
%   This function calls strengths_und_sign.m and strengths_und.m (BCT
%   toolbox).
%
%   Input:      object.name,        [required]
%               object.sign,        [optional]
%
%   Output:     M.Spos/M.Sneg,     nodal strength of positive/negative
%                                  weights (output to sigh strengths)
%               M.str,             node strength (output to no sigh strengths)
%
%   Example:    opts.f_l = 0.1;
%               opts.f_h = 0.2;
%               opts.type = 'fir';
%               g.name = 'strength';
%               g.sign = true;
%               M = run_fMRI_ts(X, opts, g);
%
%
% RICH CLUB
%
%   The rich club coefficient, R, at level k is the fraction of edges that
%   connect nodes of degree k or higher out of the maximum number of edges
%   that such nodes might share. A threshold is performed before the rich_club
%   This function calls rich_club_bu.m,
%   rich_club_wu.m, threshold_proportional.m, threshold_absolute.m, and
%   weight_conversion.m (BCT toolbox).
%   Note1: rich club outcome is not node related.
%   Note2: if the outcome has length less than the number of nodes of the
%   given timeserie and X is a matrix, zeros is put in the output just to
%   complete the size. However, remember that the output is not node based.
%
%   Input:      object.name,            [required]
%               object.k_level,         [optional]
%               object.binary,          [optional]
%               object.threshold_type   [required]
%               object.thr              [required if threshold_type is 'absolute']
%               object.threshold_p      [optional - used by threshold_type = 'proportional']
%
%   Output:     M.Rw,               rich-club curve (binary-false outcome)
%               M.R,                vector of rich-club coefficients for
%                                   levels 1 to k_level. (binary-true
%                                   outcome)
%
%   Example:    opts.f_l = 0.1;
%               opts.f_h = 0.2;
%               opts.type = 'fir';
%               h.name = 'rich_club';
%               h.binary = false;
%               h.threshold_type = 'proportional';
%               h.threshold_p = 0.5;
%               M = run_fMRI_ts(X, opts, h);
%
%
% PAGERANK CENTRALITY
%
%   The PageRank centrality is a variant of eigenvector centrality. This
%   function computes the PageRank centrality of each vertex in a graph.
%   Formally, PageRank is defined as the stationary distribution achieved
%   by instantiating a Markov chain on a graph. The PageRank centrality of
%   a given vertex, then, is proportional to the number of steps (or amount
%   of time) spent at that vertex as a result of such a process. 
%   The PageRank index gets modified by the addition of a damping factor,
%   d. In terms of a Markov chain, the damping factor specifies the
%   fraction of the time that a random walker will transition to one of its
%   current state's neighbors. The remaining fraction of the time the
%   walker is restarted at a random vertex. A common value for the damping
%   factor is d = 0.85.
%   This function calls pagerank_centrality.m (BCT toolbox).
%
%   Input:      object.name,          [required]
%               object.d,             [required]
%               object.falff,         [optional]
%
%   Output:     M.pagerank,         vectors of page rankings
%
%   Example:    opts.null_model = true;
%               i.name = 'pagerank';
%               i.d = 0.85              %required
%               M = run_fMRI_ts(X, opts, i);
%
%
% DIVERSITY_COEF_SIGN     Shannon-entropy based diversity coefficient
%
%   The Shannon-entropy based diversity coefficient measures the diversity
%   of intermodular connections of individual nodes and ranges from 0 to 1.
%   If Ci (community affiliation vector) is not provided, consensus is
%   going to be executed (note that consensus has additional required
%   parameters which must be provide) in order to generate Ci. This
%   function calls diversity_coef_sign.m and consensus_und_helper.m
%   (BCT toolbox).
%
%   Input:      object.name,         [required]
%               object.Ci            [optional]
%               object.tau,          [required if Ci is not provided]
%               object.reps,         [required if Ci is not provided]
%               object.buffsz,       [used if Ci is not provided - optional]
%
%   Output:     M.Hpos/M.Hneg,        diversity coefficient based on
%                                     positive/negative connections.
%
%   Example:    opts.type = 'none';
%               j.name = 'diversity';
%               j.tau = 0.5;
%               j.reps = 50;
%               M = run_fMRI_ts(X, opts, j);
%
%
% PARTICIPATION_COEF_SIGN     Participation coefficient
%
%   Participation coefficient is a measure of diversity of intermodular
%   connections of individual nodes. If Ci (community affiliation vector)
%   is not provided, consensus is going to be executed (note that consensus
%   has additional required parameters which must be provide) in order to
%   generate Ci. This function calls participation_coef_sign.m and
%   consensus_und_helper.m (BCT toolbox).
%
%   Input:      object.name,         [required]
%               object.Ci,           [optional]
%               object.tau,          [required if Ci is not provided]
%               object.reps,         [required if Ci is not provided]
%               object.buffsz,       [used if Ci is not provided - optional]
%
%   Output:     Ppos/Pneg,           participation coefficient from
%                                    positive/negative weights
%
%   Example:    opts.type = 'none';
%               l.name = 'participation';
%               l.tau = 0.5;
%               l.reps = 50;
%               M = run_fMRI_ts(X, opts, l);
%
%

if nargin<3, error('myApp:argChk', 'All 3 arguments are required.'); end

addpath(genpath('2014_04_05 BCT'));
addpath(genpath('helpers'));

if iscell(X)
    for i=1:length(varargin)
        
        for k=1:length(X)
            aux{k}=run_fMRI_ts_data_cell(X{k}, opts, varargin{i});
        end
        M{i} = aux;
        
        filename  = strcat(lower(varargin{i}.name), '_',opts.corrType, datestr(clock,'_mm-dd-yyyy_HH-MM-SS'), '.csv');
        filename_msg = ['The output is going to be in a file named ',filename];
        disp(filename_msg);
        
        creating_spreadsheet_data_cell(filename,varargin{i}.name, aux);
    end
    
elseif ismatrix(X)
    for i=1:length(varargin)
        M{i} = run_fMRI_ts_data_struct(X, opts, varargin{i});
    end
    
    filename  = strcat('fMRI_spreadsheet', datestr(clock,'_mm-dd-yyyy_HH-MM-SS'), '.csv');
    filename_msg = ['The output is going to be in a file named ',filename];
    disp(filename_msg);
    
    creating_spreadsheet_data_struct(filename,M);
    
else
    error('myApp:argChk', 'X  must be either a matrix (a single TS) or a cell (multiple TS).');
end

end


