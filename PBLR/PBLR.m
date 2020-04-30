function Xi = PBLR(M,group,boundary)
% <INPUT>
% M: the raw data
% group: the label information is know in advance or estimated by clustering
% boundary: the estimated boundary type. 1: exponetial function; 2: simple piecewise function; 3: sophisticated piecewise function
% <OUTPUT>
% Xi: the imputed data matrix
K = max(group);
Xs_record = cell(1,K); id_record = cell(1,K);
for i = 1:K
    id = find(group == i);
    x = M(:,id);
    id_record{1,i} = id; Xs_record{1,i} = x;
end
% compute the dropout on each group
m = size(M,1); U_record = cell(1,K);
for i = 1:K
    x = Xs_record{1,i};
    zero_ratio = zeros(m,1); avgs = zeros(m,1);
    for j = 1:m
        idz = find(x(j,:) ~= 0);
        if ~ isempty(idz)
            avgs(j) = sum(x(j,:))/length(idz); zero_ratio(j) = 1-length(idz)/size(x,2);
        else
            zero_ratio(j) = 1;
        end
    end
    ix = intersect(find(zero_ratio > 0),find(zero_ratio < 1));% do not consider all zeros or all non-zeros
    avgsx = avgs(ix); zrs = zero_ratio(ix);
    % three situations
    if boundary == 1 %exponetial function
        % order avgl
        [avgso,I] = sort(avgsx); zrso = zrs(I);
        mavg = mean(avgso); mzr = mean(zrso);
        a0 = -log(mzr+1)/(mavg^2+eps);
        ft = fittype('exp(-a*x^2)');% exponetial function
        fitobject = fit(avgso,zrso,ft,'StartPoint',a0);
        p = predint(fitobject,avgso,0.95);% 95% confidece intervel
        bound_values = exp(-p(2)*zrs.^2);
    elseif boundary == 2 % simple piecewise function
        bound_values = zeros(length(ix),1);
        for k = 1:length(ix)
            % between 0.05
            id = intersect(find(zrs> zrs(k)-0.05), find(zrs < zrs(k)+0.05));
            % corresponding avg values
            va_id = avgsx(id);
            % if the rate high than 0.8 take min as upper bound
            if zrs(k) >= 0.8
                bound_values(k) = min(va_id);
            else
                bound_values(k) = max(va_id);
            end
        end
        
    elseif boundary == 3 % sophisticated piecewise function
        bound_values = zeros(length(ix),1);
        for k = 1:length(ix)
            % between 0.05
            id = intersect(find(zrs> zrs(k)-0.05), find(zrs < zrs(k)+0.05));
            % corresponding avg values
            va_id = avgsx(id);
            if zrs(k) >= 0.8
                bound_values(k) = min(va_id);
            elseif zrs(k) >= 0.6 && zrs(k) < 0.8
                bound_values(k) = median(va_id(find(va_id<median(va_id))));%1/4
            elseif zrs(k)>0.4 && zrs(k)< 0.6
                bound_values(k) = median(va_id(find(va_id>median(va_id))));% 3/4
            else
                bound_values(k) = max(va_id);
            end
        end
    end
    U = avgs;  U(ix) = bound_values; U_record{1,i} = U;
end
%% impute by PBLR
Xi_record = cell(1,K); delta = 0;
for i = 1:K
    i
    x = Xs_record{1,i};  U1 = U_record{1,i}; U = repmat(U1,1,size(x,2));
    xi = MC_ADMM(x,delta,U);
    Xi_record{1,i} = xi;
end
% integrate to one matrix
Xi = M;
for i = 1:K
    id = id_record{1,i};
    Xi(:,id) = Xi_record{1,i};
end
