% This is a demo showing how to running PBLR using the synthetic dataset 2
clear;clc
addpath(genpath('./'))
%% load data
M = readtable('raw.txt','Delimiter','\t','ReadRowNames',true,'ReadVariableNames',true);% M is raw data in table form, rows are genes and columns are cells.
M0 = table2array(M); % data matrix
%% select informative genes for clustering
id = gene_selection(M0); % id is the index of the selected informative genes
%% clustering
M0 = log10(M0+1); % log10 transformated data
M00 = M0(id,:); % submatrix with selected genes
K = 3; % the number of desired clusters. One can also infer the number of clusters based on clustering stability index coph, see function clusteing_NMFs.m and our paper 
n = 20; % run symNMF and INMF 20 times, respectively
[group,coph] = clusteing_NMFs(M00,n,K);
%% boundary selection by visually checking the estimated boundary. by_default = sophisticated
boundary_selection(M0);
%% run PBLR with selected boundary function
% 1: exponetial function; 2: simple piecewise function; 3: sophisticated piecewise function
boundary = 3; % by default 
X = PBLR_main(M0,id,group,boundary); % return the imputated data matrix
Xi = max(10.^X-1,0);
dlmwrite('PBLR_impute.txt',Xi)


