function [X, err,iter]=MDT_HTNN_DCT(X0, tau, Pomega, miss_num)

%% default paremeters setting
dim = size(X0);
stau = tau(1);
ttau = tau(2);
wtau = tau(3);
hankel_dim=[stau, dim(1)-stau+1, ttau, dim(2)-ttau+1, wtau, dim(3)-wtau+1];
Pomegac=1-Pomega;
tol        = 1e-4; 
max_iter   = 500;
max_mu     = 1e10;
rho        = 1.2;
mu         = 1e-4;
%% variables initialization
X  = Pomega.*X0;
G  = zeros(hankel_dim); 
M  = zeros(hankel_dim); 
[v0, S] = DE_tensor_all(ones(dim),tau,[]);
D  = DE_tensor_all_adjoint(v0,S);
%% main loop
iter = 0;
while iter<max_iter
    iter = iter + 1;  
   
    %% Update G
    X1= DE_tensor_all( X, tau, S);
    G = prox_htnn_C(X1-M/mu,1/mu); 

    %% Updata X -- proximal operator of TNN
    Z=DE_tensor_all_adjoint(G+M/mu,S)./D;
    X=Pomegac.*Z+Pomega.*X0;
   
    %% Stop criterion
   MAE=(1/miss_num)*sum(abs(Pomegac.*X0-Pomegac.*X),'all');
    RMSE=1/sqrt(miss_num)*sqrt(sum((Pomegac.*X0-Pomegac.*X).^2,'all'));

    dY   =  G-DE_tensor_all(X,tau,S);    
    chg  = max(abs(dY(:)));
    if chg < tol
         break;
    end 
      
    %% Update detail display
        if iter == 1 || mod(iter, 1) == 0
            err = norm(dY(:),'fro');
            disp(['iter= ' num2str(iter) ', mu=' num2str(mu) ...
                  ', chg=' num2str(chg) ', err=' num2str(err)...
                   ', MAE=' num2str(MAE) ',RMSE=' num2str( RMSE)...
                     ]); 
        end
     
    %% Update mulipliers
    % M = M+mu*(G-mat2hankel_3D(X,tau));
    M = M+mu*dY;
    mu = min(rho*mu,max_mu);
end