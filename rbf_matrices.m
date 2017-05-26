function [PHI,P,phi0, pdim]=  rbf_matrices(Data,rbf_flag)

% Set up the matrices needed for computing the RBF parameters
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%Input:
%Data - structure array with all problem information 
%rbf_flag - indicator which rbf type should be used
%
%Output:
%PHI - matrix with RBF values of distances
%P - sample site matrix equipped with a vector of 1's
%phi0 - rbf value of distance 0
%pdim - dimension + 1 
%--------------------------------------------------------------------------

distances=pdist2(Data.S,Data.S); %compute pairwise dstances between points in S, pdist2 is MATLAB built-in function
%compute RBF values with the selected RBF type
if strcmp(rbf_flag,'cub') %cubic RBF
    PairwiseDistance=distances.^3;
elseif strcmp(rbf_flag,'lin') %linear RBF
    PairwiseDistance=distances;
elseif strcmp(rbf_flag,'tps') %thin-plate spline RBF
    PairwiseDistance=distances.^2.*log(distances+realmin);
    PairwiseDistance(logical(eye(size(distances)))) = 0;
end

PHI(1:Data.m,1:Data.m)=PairwiseDistance; %matrix with RBF values of distances
phi0 = rbfvalue(0, rbf_flag); %phi-value where distance of 2 points =0 (diagonal entries)
pdim = Data.dim + 1; 
P = [Data.S,ones(Data.m,1)];% sample site matrix equipped with a vector of 1's

end %function


