function run_fMRI_viz(M, opts, varargin)

% run_fMRI_viz(M, opts, varargin)
% 
% Inputs:   M,          Required
%                       Output from run_fMRI_ts.m.
%
%           opts,       Required
%                       Matlab structure containing options to generate
%                       brain vizualization. It can be the same opts used
%                       on run_fMRI_ts.m. However, it must contain
%                       additional informations.
%
%                           opts.atlas         Atlas to be used - required.
% 
%                           opts.colormap      Colormap - optional
%                                                   Default 'cool'
%
%                           opts.mapping       Mapping - optional
%
%           varargin,   Matlab structure responsable for determining
%                       which measure must be plotted. It must be the same
%                       parameter used on run_fMRI_ts.m.
%
%
% Note: For now, opts.mapping is not used.


if nargin<3, error('myApp:argChk', 'All 3 arguments are required.'); end
if ~iscell(M), error('myApp:argChk', 'Parameter M must be a Matlab cell.'); end
if ~isstruct(opts), error('myApp:argChk', 'Parameter opts must be a Matlab structure.'); end
if ~isfield(opts, 'atlas'), error('myApp:argChk', 'An atlas (opts.atlas) is required.'); end
if ~isfield(opts, 'colormap'), opts.colormap='cool'; end    

addpath(genpath('surfstat'));

for i=1:length(M)
    for j=1:length(varargin)
        
        if ~isstruct(varargin{j}), error('myApp:argChk', 'Every varargin parameter must be a Matlab structure.'); end
        if ~isfield(varargin{j}, 'name'), error('myApp:argChk', 'Every varargin parameter must have a field name (varargin.name).'); end
        
        measure_name = M{i}.spreadsheet{1}.name;
        if strcmpi(measure_name, 'Strength Positive'), measure_name='strength';
        elseif strcmpi(measure_name, 'Rich Club'), measure_name='rich_club';
        elseif strcmpi(measure_name, 'Diversity Positive'), measure_name='diversity';
        elseif strcmpi(measure_name, 'Participation Positive'), measure_name='participation';
        elseif strcmpi(measure_name, 'Pagerank Centrality'), measure_name='pagerank'; end
            
        if strcmpi(varargin{j}.name, measure_name) % checking the name
            for y=1:length(M{i}.spreadsheet)
               
                if (strcmpi(measure_name, 'modularity') || strcmpi(measure_name, 'consensus')) && (mod(y, 2)==0), continue; end % making sure that Q values will not be plotted.
                
                sizeM = size(M{i}.spreadsheet{y}.values);
                if sizeM(1)>sizeM(2), x = M{i}.spreadsheet{y}.values; 
                else x = M{i}.spreadsheet{y}.values'; end

                M{i}.spreadsheet{y}.name
                
                plotBrainView(x, opts);
                
            end
       end
    end
end

end


function plotBrainView(x, opts)

            figure;
            
            x(find(x==-1))=0; % taking the -1 values out before checking min(x)
            
            minX = min(x);
            maxX = max(x);
            
            x=(x-minX)./maxX;
            
            Value=zeros(size(opts.atlas));

            aux = [1:length(x)];
            x = [aux' x];

            for r=1:length(opts.atlas) 
                for s=1:max(opts.atlas) % regions
                    if opts.atlas(r)==x(s) 
                        Value(r)=x(s,2); 
                    end
                end
            end

            AVG_SURF = SurfStatReadSurf({'AVGmid_LEFT.obj','AVGmid_RIGHT.obj'});
            
            SurfStatView(Value, AVG_SURF, 'Cortical Surface Area');

            SurfStatColormap(opts.colormap);

            SurfStatColLim([0 1]);
end