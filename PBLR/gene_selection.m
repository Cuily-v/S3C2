%% gene selection
% if flag equals to 0, genes are selected by Fano factor index; if flag equals to 1, genes are selected by Gini index; if flag equals to 2,
% genes are selected by both Fano factor index and Gini index. 
% by default, flag = 0. If flag ~=0, please check the path for Rscript
function id = gene_selection(M0,flag)
if ~exist('flag','var') || isempty(flag)
    flag = 0; % select gene based on Fano factor and the mean expression
end
if flag == 0
    id = HVGs(M0);
elseif flag == 1
    % Replace the following line by the appropriate path for Rscript
    Rscript = '"C:\Program Files\R\R-3.5.0\bin\Rscript"'; % for 64-bit windows
    % Rscript = '"/usr/local/bin/Rscript"'; % for Mac OS
    filefolder = pwd;
    T = array2table(M0,'RowNames',strcat('gene',cellstr(num2str([1:size(M0,1)]'))));
    writetable(T,'raw_temporal.txt','Delimiter','\t','WriteRowNames',1);
    % Calling R's GiniIndex
    eval([' system([', '''', Rscript, ' GiniIndex.R ', '''', ' filefolder]);']);
    id = importdata('Gini_ID.txt');
else
    id1 = HVGs(M0);
    % Replace the following line by the appropriate path for Rscript
    Rscript = '"C:\Program Files\R\R-3.5.0\bin\Rscript"'; % for 64-bit windows
    % Rscript = '"/usr/local/bin/Rscript"'; % for Mac OS
    filefolder = pwd;
    T = array2table(M0,'RowNames',strcat('gene',cellstr(num2str([1:size(M0,1)]'))));
    writetable(T,'raw_temporal.txt','Delimiter','\t','WriteRowNames',1);
    % Calling R's GiniIndex
    eval([' system([', '''', Rscript, ' GiniIndex.R ', '''', ' filefolder]);']);
    id2 = importdata('Gini_ID.txt');
    id = union(id1,id2);
end


