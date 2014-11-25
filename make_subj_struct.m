function[X]=make_subj_struct(folder)

% Author: Davy Vanderweyen
% this function makes a subject structure X in the right format for the 
% graph theory fmriTB_v1 pipeline by analyzing a specific folder containing 
% one .txt file per subject and appending each to a cell of X.
% the text file files must have the brain regions as columns and the time
% series as rows. The header must not be taller than 1 row.
% you only need to specify the folder location in the function argument!

% WARNING: there cannot be other .txt files; they will mess up the resulting
% X file. Do not mess up the X files.

cd(folder);                                            % go into putative folder
files=dir('*.txt');                                    % put all .txt files in a cell array
l=length(files);                                       % evaluate the number of subjects
X=cell(1,l);                                           % create a cell array of the same length as there are subjects
for i=1:l
    s=files(i).name;                                   % set up a variable equal to a string
    [token,~]=strtok(s,'.');                           % use it to separate the subject ID from the file extension
    X{i}.ID=token;                                     % append subject name to appropriate cell
    impo=importdata(files(i).name);                    % import data from file
    X{i}.TS=impo.data(:,1:264);                        % append data to subj structure
    X{i}.Nodes=length(X{i}.TS);                        % append the number of nodes to the subj structure
end

end
