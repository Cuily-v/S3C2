%% select high variable genes (HVGs) based on its average expression and Fano factor
function id = HVGs(M,low_mu,high_mu,low_F)
if ~exist('low_mu','var') || isempty(low_mu)
    low_mu = 0.01; % select gene above this average expression level
end
if ~exist('high_mu','var') || isempty(high_mu)
    high_mu = 3.5; % select gene below this average expression level
end
if ~exist('low_F','var') || isempty(low_F)
    low_F = 0.05; % select gene above this Fano factor
end
% M is the raw count data matrix
% id is the set of select gene position
%% remove outlier genes and cells, keep all genes expressed in  >=3  cells and all cells with at least 200 detected genes
fM = zeros(size(M)); fM(M ~= 0) = 1;
sR = sum(fM,2); sC = sum(fM);
ir = find(sR >= 3); ic = find(sC >= 200);
M_new = M(ir,ic);
%% normalization:we employ a global-scaling normalization method LogNormalize that normalizes the gene expression measurements for each cell by the total expression,
%multiplies this by a scale factor (10,000 by default), and log-transforms the result:
sM = sum(M_new);
sM=sM+0;
data = log(M_new*pinv(diag(sM))*10000+1);
%% select HVGs: calculates the average expression and Fano factor for each gene, places these genes into bins,
% and then calculates a z-score for Fano factor within each bin.
mu = log(1+mean(exp(data)-1,2));
F = log(var(exp(data)-1,0,2)./mean(exp(data)-1,2));%fano因子的定义
mu(isnan(mu)) = 0;F(isnan(F)) = 0;
[Y,E] = discretize(mu,20);%把平均表达值分为20组
%[Y,E] = discretize(X,N) 将 X 中的数据分成宽度一致的 N 个 bin，还会返回 bin 边界 E。
idx = setdiff(1:20,unique(Y));%C = setdiff(A,B) 返回 A 中存在但 B 中不存在的数据，不包含重复项。C 是有序的。
if ~isempty(idx)
    E(idx+1) = [];
    Y = discretize(mu,E);
end
mean_y = grpstats(F,Y,'mean');%求每一组的fano因子的平均值
sd_y = grpstats(F,Y,'std');%求每一组的fano因子的方差
F_scaled = (F - mean_y(Y))./sd_y(Y);%对于fano因子的标准化
F_scaled(isnan(F_scaled)) = 0;
idx_pass = (mu > low_mu) & (mu < high_mu) & F_scaled > low_F;
id = ir(idx_pass);

% subplot(2,2,3)
plot(mu,F_scaled,'k.')
% title(['\fontsize{15} Biase'])
xlabel('Average expression');
ylabel('Fano factor')