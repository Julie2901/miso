function [new_p, pred_f] = mlsl(Data,lambda, g, rbf_flag)

% multi-level single linkage algorithm for finding the various minima of
% the BRF surrogate surface
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%Input:
%Data - structure array with all problem information 
%lambda, g - RBF parameters
%rbf_flag - indicator of the type of RBF to use
%
%Output:
%new_p - matrix with suggested new sample points
%pred_f - predicted objective function values corresponding to new_p
%--------------------------------------------------------------------------

%% algorithm parameters
min_func =@(x)rbf_prediction(x,Data,lambda, g, rbf_flag); %function to be minimized (here: rbf)
xrange = Data.xup- Data.xlow;     % variable ranges 
minxrange = min(xrange);  % minimum variable range
pop_size = 20;   % population size for GA
maxevals = 1000*Data.dim; %max number of RBF evaluations allowed
remevals = maxevals;     % number of remaining function evaluations
maxiter = 10;     % max number of iterations
n_per_iter = 50*Data.dim;  % number of evaluations per iteration
fract_new = 0.1;      % parameter for percentage of points from which local optimization starts
sigma = 4;        % default value for computing crit distance
maxsample = maxiter*n_per_iter; %maximum number of allowed objective function evaluations (RBF surface)

ctr_fevals = 0;   % number of calls to the function to be minimized
ctr_locminsearch = 0;   % number of local minimization runs
ctr_local = 0;    % number of distinct local minima found

dist_tol = 0.005*minxrange*sqrt(Data.dim); %min distance that new points have to be away from eah other
StartPoints = zeros(maxevals, Data.dim);   % initialize start point matrix 
new_p = zeros(maxevals,Data.dim); % initialize matrix for points to output (local minima of RBF)
pred_f = zeros(maxevals,1); % initialize vector with RBF values corresponding to local minima

start_sample = rand(maxsample, Data.dim); %randomly generated set of points
SamplePoints = repmat(Data.xlow, maxsample,1) +  repmat(Data.xup-Data.xlow, maxsample,1).*start_sample; %scale to true ranges
spsample_dist = pdist(start_sample); %pairwise distance between points
pairw_dist = spsample_dist; %initialize pairwise distances
critical_dist = ((gamma(1+(Data.dim/2))*prod(Data.xup - Data.xlow)*sigma)^(1/Data.dim))/sqrt(pi); %critical distance when 2 points are equal

iter_ctr = 0; %  iteration counter
stop = false;   % indicator of termination of MLSL using the Bayesian criterion
SampleVals = zeros(maxsample, 1);   % function values of the sample points
startp_ID = zeros(maxsample, 1);  % indicates whether a sample point has been used as a starting point      

%% optimization loop
while ~stop  && (iter_ctr < maxiter) && (ctr_fevals < maxevals) 
    iter_ctr = iter_ctr+1; %increment iteration counter
    i_start = (iter_ctr-1)*n_per_iter+1;
    i_end = min(iter_ctr*n_per_iter, (iter_ctr-1)*n_per_iter + remevals);
        
    % RBF evaluations
    for ii = i_start:i_end
        SampleVals(ii) = feval(min_func, SamplePoints(ii,:));
        ctr_fevals=ctr_fevals+1;
    end

    % update number of remaining function evaluations
    remevals = maxevals - ctr_fevals;

    % if computational budget is exhausted, then return
    if ctr_fevals >= maxevals
        new_p = new_p(1:ctr_local,:);
        pred_f = pred_f(1:ctr_local);
        return;
    end

    [~, ids] = sort(SampleVals(1:i_end));% sort the sample points according to their function values
    criticaldist = critical_dist*((log(iter_ctr*n_per_iter))/(iter_ctr*n_per_iter))^(1/Data.dim);  %  critical distance
    ctr_loc_start = floor(fract_new*iter_ctr*n_per_iter); % consider only the best points to start local minimization 
    id_matrix = zeros(ctr_loc_start, 1);
    %check if any local search start points are in critical distance
    n_selected = 0;
    for ii = 1:ctr_loc_start
        if startp_ID(ids(ii)) == 0
            select = true;
            for jj = 1:(ii-1)
                if ids(ii) < ids(jj)
                    distance = pairw_dist( (ids(ii)-1)*maxsample + (ids(jj) - (ids(ii)*(ids(ii)+1))/2) );
                else
                    distance = pairw_dist( (ids(jj)-1)*maxsample + (ids(ii) - (ids(jj)*(ids(jj)+1))/2) );
                end
                if  (distance <= criticaldist) && (SampleVals(ids(jj)) < SampleVals(ids(ii))) 
                    select = false;
                    break;
                end
            end
            if select
                n_selected = n_selected + 1;
                id_matrix(n_selected) = ids(ii);
                startp_ID(ids(ii)) = 1;
            end
        end
    end

    % uase startpoints for local min as partial population in GA 
    for ii = 1:n_selected        
        % assign starting point and value
        ctr_locminsearch = ctr_locminsearch + 1;
        StartPoints(ctr_locminsearch,:) = SamplePoints(id_matrix(ii),:);
        n_generations=max(2,round(remevals/pop_size)); %number of generaltions
        part_pop = StartPoints(ctr_locminsearch,:); %partial initial population 
        ga_options = gaoptimset('Display','off','Generations',n_generations,'PopulationSize',pop_size, 'InitialPopulation', part_pop);
        [min_point,min_value,~,out]= ga(min_func,Data.dim,[],[],[],[],Data.xlow, Data.xup,[],Data.integer,ga_options);
        ctr_fevals = ctr_fevals +out.funccount;
        remevals = maxevals - ctr_fevals; %remaining function evaluations
        
        % check if minimum point is too close to previous points
        if newp(min_point, new_p(1:ctr_local,:), dist_tol) == 1
            ctr_local = ctr_local + 1;
            new_p(ctr_local,:) = min_point;
            pred_f(ctr_local) = min_value;
        end
        
        if ctr_fevals >= maxevals %check if number of allowed evaluations has been reached
            break;
        end
    end
     
    % check termination condition
    if ((ctr_fevals >= maxevals) )
        break;
    end
    
    % check stopping condition
    e_nlocmin = ctr_local*(ctr_loc_start - 1)/(ctr_loc_start - ctr_local - 2);
    e_domain = (ctr_loc_start - ctr_local - 1)*(ctr_loc_start + ctr_local)/(ctr_loc_start*(ctr_loc_start - 1));
    if (e_nlocmin - ctr_local < 0.5) && (e_domain >= 0.995)
        stop = true;
    end
           
end

new_p = new_p(1:ctr_local,:);
pred_f = pred_f(1:ctr_local);


end %function