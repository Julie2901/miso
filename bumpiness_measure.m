function hn=bumpiness_measure(x,Data,target,  lambda, gamma,PHI,P,phi0,pdim,rbf_flag)

% compute the bumpiness of the surrogate model for a potential sample point
% x
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%Input:
%x - variable vector for which we want to compute mu
%Data - structure array with all problem information 
%lambda, gamma - RBF parameters
%PHI - matrix with RBF values of distances
%P - sample site matrix equipped with a vector of 1's
%phi0 - rbf value of distance 0
%pdim - dimension + 1 
%rbf_flag - indicator which rbf type should be used
%
%Output:
%hn - bumpiness corresponding to x
%--------------------------------------------------------------------------

R_y=sqrt(sum((Data.S(1:Data.m,:)-repmat(x,Data.m,1)).^2,2)); %compute distance between x and all already sampled points

if any(R_y < Data.tol) %point x is too close to already sampled points
    hn=0; %give the bumpiness a bad bumpiness function value -> avoid sampling at this point
else
    new_phi = rbfvalue(R_y,rbf_flag); %compute rbf value of the new point x
    n_old = Data.m; %old number of function evaluations (without x)
    kk=1;
    nwr=Data.m+1;
    PHI(nwr,1:nwr-1) = new_phi; %update rbf matrices with new potential sample point        
    PHI(1:n_old+kk-1,n_old+kk) = new_phi';
    PHI(n_old+kk,n_old+kk) = phi0;
    P(n_old+kk,1:Data.dim) = x;
    %set up matrices for solving the linear system
    A_aug = [PHI(1:Data.m+1,1:Data.m+1),P(1:Data.m+1,:);P(1:Data.m+1,:)',zeros(pdim)];
    rhs=[zeros(Data.m,1); 1; zeros(Data.dim+1,1)]; %(m+1)th variable corresponds to mu
    eta = sqrt((10^-16) * norm(A_aug, 1) * norm(A_aug, inf));
    coeff = (A_aug + eta * eye(Data.m+1 + pdim)) \[rhs];
    mu=coeff(Data.m+1);
    
    if abs(mu) <1e-6
        mu=0; %mu is too inexact, give a bad value
    elseif (mu<0) 
        mu = 0;
    end
    
    if mu== 0
        hn=100;
    else
        m0=1;
        yhat = rbf_prediction(x,Data,lambda, gamma, rbf_flag); %predict RBF value of x
        gn=(-1)^(m0+1)*mu*(yhat-target)^2; %bumpiness measure
        hn=-1/gn; %minimize -1/gn to avoid numerical difficulties at already sampled points
    end
end

end%function