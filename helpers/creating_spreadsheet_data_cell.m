function creating_spreadsheet_data_cell(filename,measure_name, data)

% data = cell
% data{i}.ID = string
% data{i}.TS = vector
% data{i}.nodes = number

if nargin<2, error('myApp:argChk', 'All 2 arguments are required.'); end

delim = ',';
existing_file = false;

if exist(filename, 'file'), existing_file=true; end

fid = fopen(filename,'a');
nodes=1:data{1}.nodes;

% Print the heading
if existing_file, fprintf(fid, '\n\n'); end
fprintf(fid, '%s', upper(measure_name));

if strcmpi(measure_name, 'modularity') || strcmpi(measure_name, 'consensus') 
    modularity = false;
    if strcmpi(measure_name, 'modularity'), modularity=true; end
    
    if modularity, fprintf(fid, '\nQ and Ci\n'); %print output name
    else fprintf(fid, '\nQ and ciu\n'); end
    fprintf(fid, 'ID/Node%sQ%s', delim,delim);
    arrayfun(@(w) fprintf(fid, '%d%s', w, delim), nodes);
    fprintf(fid, '\n');
    
    % printing data
    for i=1:length(data)
        id=data{i}.ID;
        fprintf(fid, '%s%s', id, delim);
        q=data{i}.Q;
        if modularity, C=data{i}.Ci;
        else C=data{i}.ciu; end
        k = [q C'];
        arrayfun(@(y) fprintf(fid, '%d%s', y, delim), k);
        fprintf(fid, '\n');
    end
    
    if isfield(data{1},'Ci_null_model') || isfield(data{1},'ciu_null_model')
        if modularity, fprintf(fid, '\nQ_null_model and Ci_null_model\n'); %print output name
        else fprintf(fid, '\nQ_null_model and ciu_null_model\n'); end
        fprintf(fid, 'ID/Node%sQ_null_model%s', delim,delim);
        arrayfun(@(w) fprintf(fid, '%d%s', w, delim), nodes);
        fprintf(fid, '\n');
        
        % printing data
        for i=1:length(data)
            id=data{i}.ID;
            fprintf(fid, '%s%s', id, delim);
            Q_null_model=data{i}.Q_null_model;
            if modularity, C_null_model=data{i}.Ci_null_model;
            else C_null_model=data{i}.ciu_null_model; end
            k = [Q_null_model C_null_model'];
            arrayfun(@(y) fprintf(fid, '%d%s', y, delim), k);
            fprintf(fid, '\n');
        end
    end
    
else
    
    fieldname = fieldnames(data{1});
    
    for y=3:length(fieldname) % 3 because data{i} also has ID and nodes as field
        
        %print output name
        fprintf(fid, '\n%s\n', fieldname{y});
        
        % print the nodes (as a row)
        fprintf(fid, 'ID/Node%s', delim);
        
        %nodes=1:156;
        arrayfun(@(w) fprintf(fid, '%d%s', w, delim), nodes);
        fprintf(fid, '\n');
        
        % printing data
        for i=1:length(data)
            id=data{i}.ID;
            fprintf(fid, '%s%s', id, delim);
            k=getfield(data{i}, fieldname{y});
            arrayfun(@(y) fprintf(fid, '%d%s', y, delim), k);
            fprintf(fid, '\n');
        end
    end
end

fclose(fid);

end

