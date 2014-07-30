function creating_spreadsheet_data_struct(filename, data)

% data = cell
% data{i}.name
% data{i}.values = vector

if nargin<2, error('myApp:argChk', 'All 2 arguments are required.'); end


delim = ',';

fid = fopen(filename,'w');

print_node = false;
for k=1:length(data)
    % Print the heading cells
    for i=1:length(data{k}.spreadsheet),
        
        if i==1 && ~print_node,
            
            fprintf(fid, 'Node%s%s%s', delim, upper(data{k}.spreadsheet{i}.name), delim);
            print_node = true;
            
        elseif k==length(data) && i==length(data{k}.spreadsheet),
            
            fprintf(fid, '%s', upper(data{k}.spreadsheet{i}.name));
            
        else
            
            fprintf(fid, '%s%s', upper(data{k}.spreadsheet{i}.name), delim);
            
        end;
    end
end

fprintf(fid,'\n');

%V=zeros(num_nodes,length(data));

idx = 1;
for y=1:length(data)
    for j=1:length(data{y}.spreadsheet),
        
        V(:,idx)=data{y}.spreadsheet{j}.values;
        idx=idx+1;
    end
end

V=[(1:size(V,1))' V];

for i=1:size(V,1),
    
    d=V(i,:);
    
    arrayfun(@(y) fprintf(fid,'%d%s',y,delim),d);
    
    fprintf(fid,'\n');
    
end;


fclose(fid);

end

