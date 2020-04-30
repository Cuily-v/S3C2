clear
clc
addpath('PBLR_v1.0')
addpath('Data')
% load('Test_1_Biase.mat')
% load('Test_2_Grun.mat')
% load('Test_3_Loh.mat')
% load('Test_4_Ting.mat')
load('Test_5_pollen.mat')
% load('Test_6_Zeisel.mat')
M0=in_X';
flag=0;
id = gene_selection(M0,flag)
in_X=(in_X(:,id))';
X=in_X;
A=zeros(size(X,2),size(X,2));
for i=(-1):0.01:0
Omega=logical(X);
idx=true_labs;
lambda=10^(i);
Gamma=10^(1);
gamma0 =0.5;
t_max =10;
K=12;
[me_idx,acc_i, Z, X_fill] =S3C2(X, Omega, idx,lambda, Gamma, gamma0,t_max)
A=A+Z;
end
Z=A/101;

