function  boundary_selection(M)
%% select optimal boundary function visually
% we provided three ways (i.e., three boundary functions) for estimating boundary of gene expression in single cell data.
% First, the boundary of each gene is defined as the upper one-sided 95% confidence bound by fitting the log-transformed mean gene expression value and the ratio of zeros using exp(-lambda*x^2) ;
% Second, boundary is estimated by a simple piecewise function (see Methods); Third, boundary is estimated by a sophisticated piecewise functions.
[m,n]=size(M);
M_flag = zeros(m,n); M_flag(M == 0) = 1; s_M = sum(M_flag,2);
r = s_M/n; r = max(r,0.3);
R = zeros(m,n);
for i = 1:m
    R(i,:) = binornd(1,r(i),[1,n]);
end
%generate new
M_new = M; M_new(R == 1) = 0;
% fit with various boundary function
% compute the dropout
zero_ratio = zeros(m,1);
avgs = zeros(m,1); % real values in zero space
avgl = zeros(m,1); % observed values in nonzero space
for j = 1:m
    idz = find(M_new(j,:) == 0);
    if ~ isempty(idz)
        zero_ratio(j) = length(idz)/n;
        avgs(j) = sum(M(j,idz))/length(idz);
    else
        zero_ratio(j) = 0;
    end
    if length(idz) ~= n
        avgl(j) = sum(M_new(j,:))/(n-length(idz));
    end
end
ix = intersect(find(zero_ratio > 0),find(zero_ratio < 1));% do not consider all zeros or all non-zeros
zrs = zero_ratio(ix);avgl = avgl(ix); avgs = avgs(ix);zrs0 = zrs; avgl0 = avgl; avgs0 = avgs;
% order avgl
[avgl,I] = sort(avgl); zrs = zrs(I);
mavg = mean(avgl); mzr = mean(zrs);
a0 = -log(mzr+1)/(mavg^2+eps);
ft = fittype('exp(-a*x^2)');% exponential function
fitobject = fit(avgl,zrs,ft,'StartPoint',a0);
p11 = predint(fitobject,avgl,0.95);% 95% confidece intervel

hFig1 = figure('position', [300, 50, 800, 200]);
subplot(1,3,1)
plot(fitobject,avgl,zrs),
hold on; plot(avgl,p11,'k--')
hold on;
scatter(avgs0,zrs0,4,'r')
hold on;
% xlim([0,2.5])
set(gca,'Ytick',[0 0.5 1])
box on;
ax=gca;ax.LineWidth=1.2;
title('exponential function')
xlabel('Mean expression')
ylabel('Ratio of zeros')
legend off
ylim([0 1])
% fit with piecewise function
bound_avg = zeros(length(ix),2);
for k = 1:length(ix)
    % between 0.05
    id = intersect(find(zrs0> zrs0(k)-0.05), find(zrs0 < zrs0(k)+0.05));
    % corresponding avg values
    va_id = avgl0(id);
    % if the rate samll than 0.8 take 0 as lower bound, min as
    % upper bound
    % small than 0.8: min and max
    if zrs0(k) >= 0.8
        bound_avg(k,1) = min(va_id);
    else
        bound_avg(k,1) = max(va_id);
    end
    if zrs0(k) >= 0.8
        bound_avg(k,2) = min(va_id);
    elseif zrs0(k) >= 0.6 && zrs0(k) < 0.8
        bound_avg(k,2) = median(va_id(find(va_id<median(va_id))));%1/4
    elseif zrs0(k)>0.4 && zrs0(k)< 0.6
        bound_avg(k,2) = median(va_id(find(va_id>median(va_id))));% 3/4
    else
        bound_avg(k,2) = max(va_id);
    end
end
% plot
subplot(1,3,2)
scatter(avgl0,zrs0,4,'b')
hold on;
scatter(avgs0,zrs0,4,'r')
hold on;
scatter(bound_avg(:,1),zrs0,4,'k')
% xlim([0,2.5])
box on;
ax=gca;ax.LineWidth=1.2;
set(gca,'Ytick',[0 0.8 1])
title('simple piecewise function')
xlabel('Mean expression')

subplot(1,3,3)
scatter(avgl0,zrs0,4,'b')
hold on;
scatter(avgs0,zrs0,4,'r')
hold on;
scatter(bound_avg(:,2),zrs0,4,'k')
% xlim([0,2.5])
box on;
ax=gca;ax.LineWidth=1.2;
set(gca,'Ytick',[0 0.4 0.6 0.8 1])
title('sophisticated piecewise function')
xlabel('Mean expression')
