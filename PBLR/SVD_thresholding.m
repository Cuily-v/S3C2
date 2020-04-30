function Z = SVD_thresholding(L,tau)
% SVD thresholding on L
[U,S,V] = svd(L);
St = zeros(size(S));
for i = 1:min(size(S))
    s = S(i,i);
    St(i,i) = thresh(s,tau);
end
Z = U*St*V';
    