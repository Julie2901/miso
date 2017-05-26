function yhat = rbf_prediction(CandPoint,Data,lambda, gamma, rbf_flag )


% compute the predicted objuective function value
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%Input:
%CandPoint - matrix with potential sample points 
%Data - structure array with all problem information 
%lambda, gamma - rbf parameters
%rbf_flag - indicator which rbf type should be used
%
%Output:
%yhat - predicted objective function values for all points in CandPoint
%--------------------------------------------------------------------------

[mX,nX]=size(CandPoint); %dimensions of the points where function value should be predicted
[mS,nS]=size(Data.S); %dimensions of sample sites taken so far 
R = pdist2(CandPoint,Data.S); %compute pairwise dstances between points in X and S. pdist2 is MATLAB built-in function

%compute RBF value
if strcmp(rbf_flag,'cub') %cubic RBF
    Phi=R.^3;
elseif strcmp(rbf_flag,'lin') %linear RBF
    Phi=R;
elseif strcmp(rbf_flag,'tps') %thin-plate spline RBF
    Phi=R.^2.*log(R+realmin);
    Phi(logical(eye(size(R)))) = 0;
end
   
p1 = Phi*lambda; %first part of response surface - weighted sum of RBF values of distances
p2 = [CandPoint,ones(mX,1)]*gamma; % polynomial tail of response surface
yhat=p1+p2; %predicted function value

end%function