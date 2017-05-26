function mu=inf_step(x, Data, PHI, P, phi0, pdim, rbf_flag)

% compute the value of mu in the inf step of the target value sampling
% strategy
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%Input:
%x - variable vector for which we want to compute mu
%Data - structure array with all problem information 
% PHI - matrix with RBF values of distances
%P - sample site matrix equipped with a vector of 1's
%phi0 - rbf value of distance 0
%pdim - dimension + 1 
%rbf_flag - indicator which rbf type should be used
%
%Output:
%mu - parameter value if x was a new sample point
%--------------------------------------------------------------------------

R_y=sqrt(sum((Data.S(1:Data.m,:)-repmat(x,Data.m,1)).^2,2)); %compute distance between x and all already sampled points

if any(R_y < Data.tol) %point x is too close to already sampled points
    mu=99999; %assign bad value
else
    new_phi = rbfvalue(R_y, rbf_flag); %compute rbf value of the distances of x to Data.S
    n_old = Data.m; %the old number of function evaluations
    kk=1; %equip RBF matrices with new data
    PHI(n_old+kk,1:n_old+kk-1) = new_phi;
    PHI(1:n_old+kk-1,n_old+kk) = new_phi';
    PHI(n_old+kk,n_old+kk) = phi0;
    P(n_old+kk,1:Data.dim) = x;
    %set up matrixes for solving linear system
    A_aug = [PHI(1:Data.m+1,1:Data.m+1),P(1:Data.m+1,:);P(1:Data.m+1,:)',zeros(pdim)];
    rhs=[zeros(Data.m,1); 1; zeros(Data.dim+1,1)]; %(m+1)th variable corresponds to mu

    eta = sqrt((10^-16) * norm(A_aug, 1) * norm(A_aug, inf));
    coeff = (A_aug + eta * eye(Data.m+1 + pdim)) \[rhs];
    mu=coeff(Data.m+1);
    if abs(mu) <1e-6
        mu=0;
    elseif (mu<0) %mu is too imprecise, assign bad value 
        mu=99999;
    end
end

end%function