clear all
clc

%% Import Data
addpath('PBLR')
addpath('Data')
load('Test_5_pollen.mat')

%% Data preprocessing
M0=in_X';
flag=0;
id = gene_selection(M0,flag)
in_X=(in_X(:,id))';


%% The main body of the method (lambda and Gamma require parameter sensitivity analysis, here is just an example)
X=in_X;
Omega=logical(X);
idx=true_labs;
lambda=10^(-1);
Gamma=10^(0);
gamma0 =0.5;
t_max =10;
K=12;
[me_idx,acc_i, Z, X_fill] =S3C2(X, Omega, idx,lambda, Gamma, gamma0,t_max)


