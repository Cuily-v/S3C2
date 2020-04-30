% S3C2.m
function [me_idx,acc_i, Z, X_fill] =S3C2(X, Omega, idx,lambda, Gamma, gamma0, t_max, relax, affine,  opt, T, nu1, nu2, tau, lambda_max, Gamma_max)
if nargin < 8      
    relax =1;
end
if nargin <9        
    affine =0;
end
if nargin < 10       
    opt.tol =1e-4;
    opt.maxIter =1e6;
    opt.rho =1.1;
    opt.mu_max =1e4;
    opt.norm_sr ='1';
    opt.norm_mc ='1';
end
if nargin < 11     
    T =1; % T=1 for spectral clustering and T>2 for spectral assumble clustering
end
if nargin < 12
    nu1 =0.8; % or nu1 =1.0-2.0,  e.g. 1.1, 1.2, 1.5, 2, which is used to increase the lambda and Gamma
end
if nargin < 13
    nu2 =1; % or nu2 =1
end
if nargin < 14
    tau =Inf;
end
if nargin < 15
    lambda_max =Inf; 
    Gamma_max =Inf;
end
[d,n]=size(X);
SC_assemble =0;

%% Initialization
nbcluster =max(idx);  % the number of clusters is given, or it need to be estimated.
Theta = zeros(size(X,2));
Z = zeros(size(X,2));
X = zeros(size(X)) + X;
X_fill = X; 
E = zeros(size(X));

X_norm_two =norm(X,2);
opt.mu =1.25/X_norm_two;
X_fill_old = X_fill;
Theta_old =ones(size(X,2));
lambda_t =lambda;
Gamma_t =Gamma;  
%% main loop of S3C2
t=0;
while (t < t_max)     
    
    t = t+1;
    % parameters updating
    lambda_t =min(lambda_t * nu1, lambda_max);
    Gamma_t = min(Gamma_t * nu2, Gamma_max);
    [Z, X_fill, E] =SMC_StrSR(X, Omega,Theta, lambda_t, Gamma_t, gamma0, tau, opt, X_fill, Z,  relax, affine);
    
	
    %%%Network_Diffusion
 	A=Z; 
 	A = A-diag(diag(A));
 	P = (dominateset(double(abs(A)),min(30,length(A)-1))).*sign(A);
 	DD = sum(abs(P'));
 	P = P + (eye(length(P))+diag(sum(abs(P'))));
 	P = (TransitionFields(P));
 	[U,D] = eig(P);
 	d = real((diag(D))+eps);
 	alpha = 0.8;
 	beta = 2;
 	d = (1-alpha)*d./(1-alpha*d.^beta);
 	D = diag(real(d));
 	W = U*D*U';
 	W = (W.*(1-eye(length(W))))./repmat(1-diag(W),1,length(W));
 	D=sparse(1:length(DD),1:length(DD),DD);
 	W = (W+W')/2;     
 	Z=W;    
   
    
    
    %Step II: Structure Estimation via Spectral Assemble Clustering
    if (SC_assemble)
        rho =1;CKSym = BuildAdjacency(thrC(Z,rho));
        [acc_i, Theta, ~, ~] = SpectrAssembleClustEvaluation(CKSym, idx, nbcluster, T);
    else
        [me_idx,acc_i, ~, Theta] =SpectrClustEvaluation(Z, idx, nbcluster, T, 1e-40, 1-Theta);
    end

    Theta =1-Theta;    
    %% Stop criterion checking    
    deltaC =norm(X_fill - X_fill_old,'fro')/norm(X_fill,'fro');
	deltaTheta=norm(Theta-Theta_old,'inf')
	
    if (~affine)
        tol_D =1e-5;
    else
        tol_D =5e-3;
    end
    if ((deltaTheta < 1) || (deltaC < tol_D))            
           disp(['* * * * * * STOP @ t = ',num2str(t),'. acc = ',num2str(acc_i)]);
            break;
    end   
    X_fill_old =X_fill;
    Theta_old = Theta;
end
