function X = MC_ADMM(M,delta,U,varargin)
%imputation matrix with ADMM method
% agrmin ||X||_*, s.t. ||X-Y||_F^2 <= delta
% <Inputs>:
% M:     The observed data matrix
% delta: Parameter for relaxation model
% U: the upper boundary

%        (Below are optional arguments: can be set by providing name-value pairs)
%        X_INIT:     The initial value of decision matrix X       
%        Z_INIT:     The initial value of Lagrange multiplier
%        ITER_MAX:   The maximum iterations
%        TOL:        The predefined presision


% <Outputs>:
% X :     Record the decision matrix of each iteration
%% initialization
[m,n] = size(M); X_init = zeros(m,n); Z_init = zeros(m,n); iter_max = 500; tol = 10^(-6); 
if (rem(length(varargin),2) == 1)
    error('Error:Optional parameters should always go by pairs');
else
    for i = 1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'X_INIT',              X_init = varargin{i+1};
            case 'Z_INIT',              Z_init = varargin{i+1};
            case 'ITER_MAX',            iter_max = varargin{i+1};
            case 'TOL',                 tol = varargin{i+1};
            otherwise
                error(['Unrecognized option: ',varargin{i}]);
        end
    end
end
X = X_init;  Z = Z_init; 
r = 1.6; beta = 2.5/sqrt(m*n);
X_record = cell(1,iter_max); Y_record = cell(1,iter_max); Z_record = cell(1,iter_max); stop_value = zeros(1,iter_max);
for iter = 1:iter_max
    % update Y
    O = find(M == 0);    
    B = X-1/beta*Z; 
    % project to the range
    L = zeros(size(B));
    d1 = B-L; index1 = find(d1 < 0); B(index1) = L(index1);
    d2 = U-B; index2 = find(d2 < 0); B(index2) = U(index2);   
    b_m = B-M; b_m(O) = 0; para_y = min(delta/norm(b_m,'fro'),1);
    Y = (para_y - 1)*b_m + B; Y_record{1,iter} = Y;
    % update X
    A = Y + 1/beta*Z;
    X = SVD_thresholding(A,1/beta);
    X_record{1,iter} = X;
    % update Z
    Z = Z-r*beta*(X-Y); Z_record{1,iter} = Z;
    % stop control
    if iter > 1
        stop_value(iter) = norm(X-X_record{iter-1},'fro')/norm(X_record{iter-1},'fro');
        if stop_value(iter) < tol
            break,
        end
    end
end
X = X_record{1,iter};
