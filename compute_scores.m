function xnew = compute_scores(Data,CandPoint, w_r, lambda, gamma, rbf_flag )

% computing the scores of the candidate points from which we select the
% next sample point
% This function is called by the sampling strategies that use a candidate
% point sampling approach (cp.m, rs.m, cptv.m, cptvl.m)
% We compute 2 scores: one derived from the objective function value
% predictions by the RBF, one derived from the distance of the candidates
% to the set of already sampled points
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%Input:
%Data - structure array with all problem information 
%CandPoint - matrix with candidate sample points
%w_r weight for the surrogate model criterion (weight for distance
%criterion is 1-w_r)
%lambda, gamma - RBF parameters
%rbf_flag - indicator of the type of RBF model wanted
%
%Output:
%xnew - new sample point at which the next function evaluation will be done
%--------------------------------------------------------------------------

[mX,nX]=size(CandPoint); %dimensions of the points where function value should be predicted
[mS,nS]=size(Data.S); %dimensions of sample site matrix (already evaluated points)  
R = pdist2(CandPoint,Data.S); %compute pairwise dstances between points in X and S. pdist2 is MATLAB built-in function

%compute RBF matrix values
if strcmp(rbf_flag,'cub') %cobic RBF
    Phi=R.^3;
elseif strcmp(rbf_flag,'lin') %linear RBF
    Phi=R;
elseif strcmp(rbf_flag,'tps') %thin-plate spline RBF
    Phi=R.^2.*log(R+realmin);
    Phi(logical(eye(size(R)))) = 0;
end
   
p1 = Phi*lambda; %first part of response surface - weighted sum of distances
p2 = [CandPoint,ones(mX,1)]*gamma; % polynomial tail of response surface
yhat=p1+p2; %predicted function value

%% surrogate model score
%scale predicted objective function values to [0,1]
min_yhat = min(yhat); %find min of predicted objective function value
max_yhat = max(yhat); %find maximum of predicted objective function value
if min_yhat == max_yhat  %compute scaled objective function value scores
    scaled_yhat=ones(length(yhat),1);
else
    scaled_yhat = (yhat-min_yhat)/(max_yhat-min_yhat);
end

%% distance score
dist=pdist2(CandPoint,Data.S(1:Data.m,:))'; %compute the pairwise distances between the candidate points and the already sampled points
%scale distances to already evaluated points to [0,1]
min_dist = (min(dist,[],1))'; % distance of every candidate point to the set of already sampled points  
max_min_dist = max(min_dist); %maximum of distances
min_min_dist = min(min_dist); %minimum of distances
if  max_min_dist == min_min_dist  %scale distance values to [0,1]
        scaled_dist =ones(length(min_dist),1);
else
        scaled_dist = (max_min_dist-min_dist)/(max_min_dist-min_min_dist);
end

%% compute weighted score for all candidates
score = w_r*scaled_yhat + (1 -w_r)*scaled_dist;

%assign bad scores to candidate points that are too close to already sampled
%points
score(min_dist < Data.tol) = Inf; 
%find candidate with best score (lowest score) -> becomes new sample point
[~,minID] = min(score);
xnew = CandPoint(minID,:);  %new sample point


end %function