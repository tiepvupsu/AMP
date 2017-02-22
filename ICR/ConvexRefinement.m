function [x_j] = ConvexRefinement(x_j_old,Avg_x,y,A,AtA,Rho,Lambda,MaxInterIter,Unconstrained)

if Unconstrained ==1
    Scale = 200;
    rho_ADMM =0.1;
    F = Scale*diag(Rho(:,1)./ abs(Avg_x) );
    x_k_old = zeros(size(A,2),1);
    z_k_old = zeros(size(A,2),1);
    u_k_old = zeros(size(A,2),1);
    
    for Internal_iter = 1:MaxInterIter
        x_k = (AtA+rho_ADMM*(F'*F)+Lambda*eye(size(AtA))) \  ( A'*y + rho_ADMM * F'*(z_k_old - u_k_old) );
        z_k = wthresh(F*x_k+u_k_old,'s',1/(2*Scale*rho_ADMM) );
        u_k = u_k_old + F*x_k - z_k;
        
        x_k_old = x_k;
        z_k_old = z_k;
        u_k_old = u_k;
    end
    x_j = F\z_k;
    
else
    
    options = optimset('Algorithm','interior-point-convex','display','off' );
    f=Rho' ./ Avg_x'  - 2*y'*A;
    H=1*( AtA + Lambda* eye(size(AtA)) );
    [x_j,fval1P,exitflag,~,lagrange]= quadprog(H,f,[],[],[],[],zeros(1,size(A,2)),[],[],options);
    
end
