function [group,coph] = clusteing_NMFs(M,n,K)
% perform unsupervised clustering of single cell data based on NMF
% Inputs:
% M: the data matrix for clustering
% n: repeat times of NMF
% K: the number of clusters
% Outputs:
%    group: the cluster label of each cell
%    coph; clustering stability index. This index can be used to determine the cluster number.

% distance part do PCA transformation
t=3; Ds = cell(1,t); Ds{1,1} = pdist2(M',M','correlation'); Ds{1,2} = pdist2(M',M','Spearman'); Ds{1,3} = pdist2(M',M','cosin');
% compute the similarity matrix
As = cell(1,t);
for i = 1:t
    D = Ds{1,i}; As{1,i} = exp(-D./max(D(:)));
end
% symNMF on these three data matrix with repeat times equaling to 20
Result = cell(1,t+1);
for i = 1:t
    A = As{1,i};
    Re = cell(1,n);
    for j = 1:n
        H = symnmf_anls(A, K);
        Re{1,j}.H = H;
    end
    Result{1,i} = Re;
end
% compute the INMF result
A = zeros(size(M)); index = find(M ~= 0); A(index) = 1;
Re = cell(1,n);
tol = 10^(-6); maxiter = 500;
for i = 1:n
    i
    [~,H] = INMF(M,A,K,tol,maxiter);
    Re{1,i}.H = H';
end
Result{1,t+1} = Re;
% compute the consensus matrix
N = size(M,2);
Cs = cell(1,t+2);
for flag = 1:t+1
    Re = Result{1,flag};
    label_record = cell(1,n);
    for i = 1:n
        V = Re{1,i}.H';
        label = zeros(1,N);
        for j = 1:N
            label(j) = find(V(:,j) == max(V(:,j)));
        end
        label_record{1,i} = label;
    end
    % compute consensus matrix
    consensus = zeros(N);
    for i = 1:n
        label = label_record{1,i}; C = zeros(N);
        for j = 1:N
            for k = j:N
                if label(j) == label(k)
                    C(j,k) = C(j,k)+1;
                    C(k,j) = C(j,k);
                end
            end
        end
        consensus = consensus + C;
    end
    Cs{1,flag} = consensus/n;
end
% the last one is the sum of the first three
C = zeros(N);
for i = 1:t
    C = C + Cs{1,i};
end
Cs{1,t+2} = C/t;
% compute the coph of consensus matrix
clusts = cell(1,t+2); cophs = zeros(1,t+2);
for i = 1:t+2
    [~,clusts{1,i},~,cophs(i)] = nmforderconsensus0(Cs{1,i},K);
end
% when the first three consensus matrix is small, select the last one
cophid = cophs(t+1:t+2);
if abs(max(cophid)-min(cophid)) > 0.1
    id = find(cophid == max(cophid));
    group = clusts{1,t+id}; coph = cophs(t+id);
else
    C = (Cs{1,t+1}+Cs{1,t+2})/2;
    [~,group,~,coph] = nmforderconsensus0(C,K);
end

