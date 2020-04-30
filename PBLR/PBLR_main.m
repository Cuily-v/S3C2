function Xi = PBLR_main(M,id,group,boundary)
%% run PBLR on each group
M1 = M(id,:); [m,n] = size(M); nid = setdiff(1:m,id); M2 = M(nid,:);
group2 = ones(1,n);
Xi1 = PBLR(M1,group,boundary);
disp('remaining submatrix')
Xi2 = PBLR(M2,group2,boundary);
% integrate
Xi = zeros(m,n);
Xi(id,:) = Xi1; Xi(nid,:) = Xi2;
