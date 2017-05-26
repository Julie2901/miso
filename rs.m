function  Data = rs(Data)

% sampling procedure that creates random candidate points by perturbing ALL
% variable values of the best point found so far. restarts from scratch if
% it thinks a local min was found
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%Input:
%Data - structure array with all problem information
%
%Output:
%Data - updated structure array with all problem information
%--------------------------------------------------------------------------

%% initializations
rs.ncand= min(100*Data.dim,5000); % number of candidate points
rs.xrange = Data.xup - Data.xlow; %variable range in every dimension
rs.minxrange = min(rs.xrange); %smallest variable range
rs.sigma_stdev_default = 0.2*rs.minxrange; %default variable perturbation range    
rs.sigma_stdev = rs.sigma_stdev_default;  % current perturbation range
rs.sigma_min=0.2*(1/2)^6*rs.minxrange; %minimum perturbation range
rs.maxshrinkparam = 5;% maximal number of perturbation range reductions 
rs.failtolerance = max(5,Data.dim); %threshold for consecutively failed improvement trials 
rs.succtolerance =3; %threshold for consecutively successful improvement trials
rs.iterctr = 0; % number of iterations
rs.shrinkctr = 0; % number of times sigma_stdev was halved
rs.failctr = 0; % number of consecutive unsuccessful iterations
rs.localminflag = 0;  % indicates whether or not xbest is a local minimum
rs.succctr=0; % number of consecutive successful iterations
restart = false; %restart becomes true after local min encountered and enough evaluations left to restart from scratch
rs.S = []; %initialize sample site matrix, will be merged with Data.S after local min is found
rs.Y = []; %initialize vector with objective function values, will be merged with Data.Y after local min is found
rs.T = []; %initialize vector with function evaluation times, will be merged with Data.T after local min is found
rs.xbest = []; %best point found in current restart
rs.fbest = []; %best value found in current restart
rs.m = 0; %number of evaluations done in current restart
rs.w_r = 0.95; %weight for surrogate model criterion, always emphasizes surrogate model score
start_nr = 1; %number of restarts


%% optimization steps
% do until max number of f-evals reached or local min found
while rs.m+Data.m < Data.maxeval
    if restart %after we found a local minimum, we restart from scratch 
        restart = false;
        %save data obtained so far, create new desgin
        start_nr =start_nr +1; %update counter for restarts
        rs.S = [rs.S; Data.S]; %sample site matrix
        rs.Y = [rs.Y; Data.Y]; %vector with objective function values
        rs.T = [rs.T; Data.T]; %vector with function evaluation times
        rs.xbest = [rs.xbest; Data.xbest]; %best point found
        rs.fbest = [rs.fbest; Data.fbest]; %best function value found
        rs.m = rs.m + Data.m; %number of function evaluations
        %reset parameters
        rs.sigma_stdev = rs.sigma_stdev_default; %reset perturbation range
        rs.shrinkctr = 0; % reset number of times perturbation range was decreased
        rs.failctr = 0; % reset number of consecutive unsuccessful iterations
        rs.localminflag = false;  % indicates whether or not xbest is at a local minimum
        rs.succctr=0; % reset number of consecutive successful iterations
        rs.iterctr = 0; %reset iteration counter
        if strcmp(Data.init_design, 'own')
           Data.init_design = 'slhd'; 
        end
        if Data.number_startpoints < Data.dim + 1 %may be because user previously supplied an initial design
            Data.number_startpoints = Data.dim + 1;
        end
        %compute number of start points (either the user-given number, or
        %the number of evaluations left in the budget)
        Data.number_startpoints = min(Data.number_startpoints, Data.maxeval - rs.m);
        %generate new starting design
        if strcmp(Data.init_design, 'slhd') %symmetric latin hypercube
            InitialPoints = slhd(Data); 
        elseif strcmp(Data.init_design, 'lhs') %matlab's lhsdesign
            InitialPoints = lhsdesign(Data.number_startpoints,Data.dim);
        end

        %scale S to true ranges
        Data.S=repmat(Data.xlow, Data.number_startpoints,1) + repmat(Data.xup-Data.xlow,Data.number_startpoints,1).*InitialPoints; 
        %round integer variables to integers
        Data.S(:,Data.integer)=round(Data.S(:,Data.integer));
        %do expensive evaluations --sequentially
        Data.Y=zeros(Data.number_startpoints,1); %initialize vector with function evaluations
        Data.T = zeros(Data.number_startpoints,1); %initialize vector with evaluation times
        for ii = 1: Data.number_startpoints %go through all startpoints
            fevalt = tic;  %start timer for function evaluation
            Data.Y(ii)=feval(Data.objfunction, Data.S(ii,:)); %expensive function evaluation
            Data.T(ii) = toc(fevalt); %stop timer for function evaluation
        end

        Data.m=size(Data.S,1); %number of function evaluations done
        [Data.fbest,fbest_ID]=min(Data.Y); %best function value so far
        Data.xbest=Data.S(fbest_ID,:); %best point so far    
        if rs.m+Data.m >= Data.maxeval %check if total number of evaluations (from all restarts), exceeds allowed number of evaluations
            break
        end
    end 
        
    rs.iterctr = rs.iterctr + 1; %increment iteration counter
    fprintf('\n Restart number: %d \n',start_nr);
    fprintf('\n Iteration: %d \n',rs.iterctr); %print some user information
    fprintf('\n Number of function evaluations: %d \n', rs.m +Data.m);
    fprintf('\n Best function value so far: %d \n', Data.fbest)
    
    if ~isempty(rs.fbest) %rs.fbest is a vector with the best values found in all restarts
        bp = min(rs.fbest); 
    else
        bp = Data.fbest;
    end
    fprintf('\n Overall best value: %f \n', bp);
    
    %compute RBF parameters
    if strcmp(Data.surrogate, 'rbf_c') %cubic RBF
        rbf_flag = 'cub';
        [lambda, gamma] = rbf_params(Data, rbf_flag);
    elseif strcmp(Data.surrogate, 'rbf_l') %linear RBF
        rbf_flag = 'lin';
        [lambda, gamma] = rbf_params(Data, rbf_flag);
    elseif strcmp(Data.surrogate, 'rbf_t') %thin plate spline RBF
        rbf_flag = 'tps';
        [lambda, gamma] = rbf_params(Data, rbf_flag);
    else
        error('rbf type not defined')
    end
              
    %generate candidate points by perturbing all variable values of the
    %best point found so far
    CandPoint = repmat(Data.xbest,rs.ncand,1); 
    for ii =1:rs.ncand    
        for jj =1:Data.dim            
            if ismember(jj,Data.integer) %integer perturbation has to be at least 1 unit 
                rr=randn(1);
                s_std=sign(rr)*max(1,abs(rs.sigma_stdev*rr));
            else %perturbation of continuous variables can be smaller than 1
                s_std= rs.sigma_stdev*randn(1);
            end
            CandPoint(ii,jj) = CandPoint(ii,jj) +s_std;
            %if candidate point is outside lower and upper variable bounds,
            %do reflection over boundary 
            if CandPoint(ii,jj) < Data.xlow(jj)
                CandPoint(ii,jj) = Data.xlow(jj)+ (Data.xlow(jj)-CandPoint(ii,jj));
                if CandPoint(ii,jj) >Data.xup(jj)
                    CandPoint(ii,jj) = Data.xlow(jj);
                end
            elseif CandPoint(ii,jj) > Data.xup(jj)
                CandPoint(ii,jj) = Data.xup(jj)- (CandPoint(ii,jj)-Data.xup(jj));
                if CandPoint(ii,jj) <Data.xlow(jj)
                    CandPoint(ii,jj) = Data.xup(jj);
                end
            end            
        end
    end
    CandPoint(:,Data.integer)=round(CandPoint(:,Data.integer)); %round integer variables
    
    xnew = compute_scores(Data,CandPoint,rs.w_r, lambda, gamma, rbf_flag); %select the best candidate
     
    % clear unnecessary variables
    clear CandPoint;
    
    fevalt = tic;  %start timer for function evaluation
    fnew = feval(Data.objfunction,xnew); %expensive evaluation
    timer = toc(fevalt); %stop timer for function evaluation
    Data.m=Data.m+1; %update the number of function evaluations
    Data.T(Data.m,1) = timer; %update the vector with evaluation times
    Data.S(Data.m,:)=xnew; %update the matrix with sample sites
    Data.Y(Data.m,1)=fnew; %update the vector with function values

    if fnew < Data.fbest %update best point found so far if necessary
        if (Data.fbest - fnew) > (1e-3)*abs(Data.fbest)
            % "significant" improvement
            rs.failctr = 0;
            rs.succctr=rs.succctr+1; %update success counter
        else
            %no "significant" improvement
            rs.failctr = rs.failctr + 1; %update fail counter
            rs.succctr=0;
        end  
        Data.xbest = xnew; %best point found so far
        Data.fbest = fnew; %best objective function value found so far
    else
        rs.failctr = rs.failctr + 1; %update fail counter
        rs.succctr=0;
    end

    
    % check if algorithm is in a local minimum
    rs.shrinkflag = 1;      
    if rs.failctr >= rs.failtolerance %check how many consecutive failed improvement trials
        if rs.shrinkctr >= rs.maxshrinkparam %threshold of allowable perturbation range reductions has been reached
            rs.shrinkflag = 0;
            disp('Stopped reducing perturbation range because the maximum reduction has been reached.');
        end
        rs.failctr = 0;
        
        if rs.shrinkflag == 1 %reduce perturbation range 
            rs.shrinkctr = rs.shrinkctr + 1;
            rs.sigma_stdev = max(rs.sigma_min,rs.sigma_stdev/2);
            disp('Reducing perturbation range to half!');
        else  
            rs.localminflag = true;
            disp('Local minimum, restarting the algorithm from scratch.');
        end
    end

    if rs.succctr>=rs.succtolerance %check if number of consecutive improvements is large enough
        rs.sigma_stdev=2*rs.sigma_stdev;%increase perturbation range
        rs.succctr=0;
    end
    if rs.localminflag %local min found, restart from scratch
        restart = true;
    end
    
end %while
  
%collect only the data for actual function evaluations that have been done
%(Data.m may be lower than maxeval, in which case the algorithm restarts)
if ~isempty(rs.S) %in this case, there was a restart and we have to put all evaluated points into Data.S, etc
    Data.S=[rs.S;Data.S];
    Data.Y = [rs.Y; Data.Y];
    Data.T = [rs.T;Data.T];
    Data.NumberFevals=Data.m+rs.m;
else %algorithm did not restart
    Data.S=Data.S(1:Data.m,:);
    Data.Y=Data.Y(1:Data.m,:);
    Data.T=Data.T(1:Data.m,:);
    Data.NumberFevals=Data.m;
end
end %function