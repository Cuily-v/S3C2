% SpectrClustEvaluation.m
%   Run spectral embedding and repeat k-means for T times
%
% Spectral clustering evaluation
%
%  2013-02-8
%  By Li Chun-Guang
% Revised: Sept. 19, 2013


function [me_idx,acc,acc_std, Theta] =SpectrClustEvaluation(L,idx0,k,T,epsilon, Theta_ini)
if nargin<3
	k=1;
end
if nargin<4
    T =20;  
end
if nargin<5
    epsilon =1e-40;
end
if nargin < 6
    Theta_ini =[];
end
if isempty(Theta_ini) || rank(Theta_ini) < 2  
    warmstart =0;
else
    warmstart =1;
    previous_idx = find_idx_from_Theta(Theta_ini, k);
end

acc_tmp =zeros(1,T);
%idx_arr =zeros(T,size(L,2));
acc =0;
acc_std =0;

%% Check full zero L
if norm(L,'fro') > epsilon
    
    L = abs(L) + abs(L');
    %% check the full-zero columns
    sumL =sum(L,1);
    tmp_idx =find(sumL < epsilon, 1);
    if (~isempty(tmp_idx))
        L(:,tmp_idx) =1e-4*abs(randn(size(L,1),length(tmp_idx)));
        L = (L + L') / 2;
    end
    
    %% Building graph Laplacian and do Spectral Embedding
    %L = BuildAdjacency(thrC(L,1));
    D = diag(sum(L,2).^(-1./2)); 
    L = eye(size(L,1)) - D*L*D;
    [U,~,~] = svd(L);
    V = U(:,end - k + 1 : end);
    
    %% Repeat k-means T times
    for ii =1 :T
        %idx = clu_ncut(L,k);     
        %idx = kmeans(V, k,'emptyaction','singleton','replicates',20,'display','off');% 10
        %acc_tmp(ii) = compacc((idx)' , idx0);  
        %% via kmeans with  'warmstart'        
        %idx = kmeans(V, k,'emptyaction','singleton','replicates',20,'display','off');% 10
        
        if (~warmstart) 
            idx  =kmeans(V, k,'emptyaction','singleton','replicates',20,'display','off');
        else
            %idx= kmeans(V, k, 'Start', 'warmstart', 'initiallabel', previous_idx, 'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
            idx= kmeans(V, k,'emptyaction','singleton','replicates',20,'display','off');     
            %previous_idx = idx;     
        end
        %acc_tmp(ii) = compacc((idx)' , idx0);
        %acc_tmp(ii)  =1- Misclassification(idx , idx0');  
        acc_tmp(ii) = evalAccuracyHungarian((idx)' , idx0);  %%%%计算ACC值的函数
        
    end
    
    Theta =form_structure_matrix(idx);
    %idx_assemble =Clustering_assemble(idx_arr,size(L,2),nbcluster);
    %acc_assemble = compacc(idx_assemble,idx0)
    acc = acc + mean(acc_tmp);
    me_idx=idx;
    acc_std =acc_std + std(acc_tmp);
    disp(num2str(mean(acc_tmp)));
    
else
    Theta =zeros(size(L));
    acc = acc + 0;
    acc_std =acc_std + 1;
    disp('failed in the bad condition of L...');
    
end

function M = form_structure_matrix(idx,n)
if nargin<2
    n =size(idx,2);
end
M =zeros(n);
id =unique(idx);
for i =1:length(id)
    idx_i =find(idx == id(i));
    M(idx_i,idx_i)=ones(size(idx_i,2));
end

function [idx]= find_idx_from_Theta(M,k,thresh)
if nargin<3
    thresh=0.5;
end
logic_M =(M >thresh);
n =size(M,2);
id =1:max(k,n);
% Find the k connected components in the graph
markers =ones(1,n);
%idx_cell =cell(1,k);
idx =zeros(1,n);
c =0;
%% Detect the biggest groups and assign........
for i=1:n
    if (markers(i)) % find the index which are still turn-on.
        c =c+1;
        j_idx =find(logic_M(i,:)==1);
        idx(j_idx) =id(c);
        %idx_cell{c} =j_idx;
        markers(j_idx) =0; % turn-off all markers which are connecting with i
    end
end