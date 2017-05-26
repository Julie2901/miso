function [lambda, gamma] = rbf_params(Data, rbf_flag)

% computing the radial basis function parameters
% This function is called by all sampling strategies whenever a new data
% point (x, f(x)) has bee obtained.
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%Input:
%Data - structure array with all problem information 
%rbf_flag - indicator of the type of RBF model wanted
%
%Output:
%lambda - weight parameters in the sum of RBF values (first term)
%gamma - parameters of polynomial tail
%--------------------------------------------------------------------------


distances=pdist2(Data.S,Data.S); %compute pairwise dstances between points in S, pdist2 is MATLAB built-in function
if strcmp(rbf_flag,'cub') %cubic RBF
    PairwiseDistance=distances.^3;
elseif strcmp(rbf_flag,'lin') %linear RBF
    PairwiseDistance=distances;
elseif strcmp(rbf_flag,'tps') %thin-plate spline RBF
    PairwiseDistance=distances.^2.*log(distances+realmin);
    PairwiseDistance(logical(eye(size(distances)))) = 0;
end


PHI(1:Data.m,1:Data.m)=PairwiseDistance; %matrix with radial bsais function values
P = [Data.S,ones(Data.m,1)];% sample site matrix equipped with vector of ones
[m,n]=size(Data.S); %determine how many points are in S and what is the dimension


A=[PHI,P;P', zeros(n+1,n+1)]; %set up matrix for solving for parameters
RHS=[Data.Y;zeros(n+1,1)]; %right hand side of linear system
params=A\RHS; %compute parameters
lambda=params(1:m); %weights 
gamma=params(m+1:end); %parameters of polynomial tail


end





