% SMC_StrSR.m
function [Z, C, E, N] =SMC_StrSR(X, Omega,Theta, lambda, Gamma, gamma0, tau, opt, X_init, Z_init,  relax, affine)
if nargin< 6
    gamma0 =0;
end 
if nargin< 7
    tau =Inf; % by default, X is noise-free.
end
if nargin < 8
    opt.tol =1e-4;
    opt.maxIter =1e6;  
    opt.rho =1.1;
    opt.mu_max =1e4;
    opt.norm_sr ='1';
    opt.norm_mc ='1';    
end   
if nargin < 9
    X_init =zeros(size(X));
end
if nargin < 10    
    Z_init =zeros(size(X,2));
end
if nargin < 11    
    relax =1;
end
if nargin < 12    
    affine =0;
end

[d, n] = size(X);

%% Initializing optimization variables
C = X_init; 
C_old = C;        
E = sparse(d,n);  
N = E;
In = eye(n);
Z = Z_init; 
Z_old = Z;
Y1 = zeros(d,n);
Y2 = Y1;
Zeros =zeros(d,n);
e_zeros =zeros(1,n);

P_Omega_star_P_Omega =logical(Omega);
P_Omega = find(Omega);  


X_norm_inf = norm( X(:), inf);
%% Find idx from Theta
k =rank(double(1-Theta)); 
thresh =0.8;
if  (relax ==2) || (sum(double(Theta(:))) <0.5)  
    idx = ones(1,n);
    k = max(idx);
else
    idx = find_idx_from_Theta(1-Theta, k, thresh);
end

%% Start main loop
iter = 0;
mu = opt.mu;      
rho =opt.rho;%Updating strategy for the penalty paramter \mu
mu_max = opt.mu_max;
tol = opt.tol; %tol =0.05; epsilon =1e-3;
maxIter = opt.maxIter;
while iter < maxIter
    iter = iter + 1;

    %% Solving C and N by linearized ADM: Structured Matrix Completion
    %% 1. updating C
    Z_norm_inf =norm(Z(:), inf);
    if (sum(Omega(:)) == d * n)        
        C = X; % If there is no missing entries, i.e. |Omega| = n*d, it is smart to let C = X;
    else        
        I_Z = In - Z;
        I_Z_norm_two =(normest(I_Z,0.1))^2;
        eta =1.01 + I_Z_norm_two;
		
        temp =C + N - X - Y1/mu; 
        Zeros(P_Omega) = temp(P_Omega_star_P_Omega);
        temp = C - (Zeros + Y2*( In - Z')/mu +(C - C*Z -E)*(In - Z'))/eta;

        for ii =1:k
            idx_ii =(idx ==ii);
            tmp_ii =temp(:,idx_ii);

            [U,sigma,V] = svd(tmp_ii,'econ');    

            sigma = diag(sigma);
            svp = length(find(sigma >= 1/(mu*eta)));
            if svp >= 1
                sigma = sigma(1:svp) - 1/(mu*eta);
            else
                svp = 1;
                sigma = 0;
            end        
            C(:,idx_ii) = U(:,1:svp)*diag(sigma)*V(:,1:svp)';
        end
        
    end
    
   %% 2. updating N (However N may not be activated forever: in noise free case tau is Inf and it noise case we prefer to 
   %                  take a relaxion ||P_\Omega(C-X)|| rather than introducing N)
    eta =1;
    switch opt.norm_mc
        case '21'
            temp =(N + C - X - Y1/mu)/eta;
            Zeros(P_Omega) = temp(P_Omega_star_P_Omega);    
            temp = N-Zeros;
            N = solve_l1l2(temp, tau/(mu*eta));

        case 'fro' % Need a check...
            temp =(N + C - X - Y1/mu)/eta;
            Zeros(P_Omega) = temp(P_Omega_star_P_Omega);
            N = (mu/(2*tau + mu*eta)) *(eta*N - Zeros);  

        case '1'
                temp =(N + C - X - Y1/mu)/eta;
                Zeros(P_Omega) = temp(P_Omega_star_P_Omega);
                temp = N - Zeros;
                N = max(0,temp - tau/(mu*eta)) + min(0,temp + tau/(mu*eta));  
    end
    
    %% Solving Z and E by linearized ADM: Structured Sparse Representation
    
    C_norm_inf =norm(C(:), inf);
    
    C_norm_two =normest(C, 0.1);
    xi =C_norm_two^2 +1.01;    

    %% 3. updating Z: simultaneously and separately 
    Z_t = Z;
    [u, s, v] = svd(Z_t);
    temp = Z - ( C' * ( C*Z + E - C - Y2/mu ) ) /xi- (Gamma*(u*v'))/(mu*xi);
    Z = max(0, temp - ( Gamma * Theta + gamma0) / (mu * xi) ) + min(0, temp +  ( Gamma * Theta + gamma0) /(mu* xi) ); 
    Z = Z - diag(diag(Z)); % to force Zii = 0
   
	
	
    %% 4. updating E
    temp = C - C*Z + Y2/mu;
    switch opt.norm_sr
        case '21'
            E = solve_l1l2(temp,lambda/mu); 

        case 'fro'
            E = (mu /(2*lambda + mu))*temp;

        case '1'
            E = max(0,temp - lambda/mu) + min(0,temp + lambda/mu);
    end
    if (affine)
        E(1,:) = e_zeros; % used for the problem with affinity constraint when augmenting C with a full-one row
    end
    
    %% Checking convergence conditions:   
    temp = X -C -N;
    leq1 = temp(P_Omega);        
    stopC1 =max(abs(leq1(:)))/X_norm_inf;
    leq2 = C - C*Z - E;
    stopC2 =max(abs(leq2(:)))/C_norm_inf; 
   
    %% Checking whole looping convergence conditions:
    leq3 = Z_old - Z;
    leq4 = C_old - C;
    stopC3 =max(abs(leq3(:)))/Z_norm_inf;  
    stopC4 =max(abs(leq4(:)))/C_norm_inf; 
    
    
    %% Checking convergence conditions:
    if (stopC1 < tol) && (stopC2<tol) ||(stopC3 < tol) && (stopC4<tol) || mu > mu_max %1e4%mu >1e6  % protect the over-looping|| mu >1e3;
        break;
    else
        
   %% Updating Y1, Y2, and mu
        temp =X - C - N;
        Zeros(P_Omega) = temp(P_Omega_star_P_Omega);
        Y1 = Y1 + mu * Zeros;
        Y2 = Y2 + mu * (C - C*Z -E);
               
        %% Updating strategy for the penalty paramter \mu
        mu = mu * rho;
        
        %% Adaptive updating strategy for the penalty paramter \mu
        %if ( mu * max( eta^0.5 * max(abs(C - C_old)), xi^0.5 * max(abs(Z - Z_old)))/X_norm_inf <= epsilon)
        %    mu = mu * rho;
        %end
        C_old =C;
        Z_old =Z;
        
    end    
end

function [idx]= find_idx_from_Theta(M,k,thresh)
if nargin<3
    thresh=0.5;
end
logic_M =(M >thresh);
n =size(M,2);
id =1:max(k,n);
% Find the k connected components in the graph
markers =ones(1,n);
%idx_cell =cell(1,k);
idx =zeros(1,n);
c =0;
%% Detect the biggest groups and assign........
for i=1:n
    if (markers(i)) % find the index which are still turn-on.
        c =c+1;
        j_idx =find(logic_M(i,:)==1);
        idx(j_idx) =id(c);
        %idx_cell{c} =j_idx;
        markers(j_idx) =0; % turn-off all markers which are connecting with i
    end
end



function [E] = solve_l1l2(W,lambda)
n = size(W,2);
E = W;
for i=1:n
    E(:,i) = solve_l2(W(:,i),lambda);
end

function [x] = solve_l2(w,lambda)
% min lambda |x|_2 + |x-w|_2^2
nw = norm(w);
if nw>lambda
    x = (nw-lambda)*w/nw;
else
    x = zeros(length(w),1);   
end
