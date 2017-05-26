function Data = tv(Data)

% Target value sampling strategy
% Define a target value and selet the point at which it is most likely that
% the target value will be achieved
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
%Data - updated structure array
%--------------------------------------------------------------------------


%% parameters and initialization for target value search
tv.iterctr = 0; %initialize iteration counter
tv.w_j=[0, 0.50,  2, 100]; %weight factors used in global grid search,
tv.w_counter=0;
tv.n=10; %cycle length
tv.sample_sequence=[0,(1:tv.n),tv.n+1]; % 0 = inf step, (1:tv.n) = cycle step local, tv.n+1 = cycle step global
gaoptions=gaoptimset('PopulationSize',10,'Display','off'); %set options for genetic algorithm

%% optimization steps
while Data.m < Data.maxeval %do until budget of function evaluations exhausted 
    %compute Rbf parameters
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

    tv.iterctr=tv.iterctr+1; %update iteration counter
    fprintf('\n Iteration: %d \n',tv.iterctr); %print some user information
    fprintf('\n Number of function evaluations: %d \n', Data.m);
    fprintf('\n Best function value so far: %d \n', Data.fbest);
    
       
    %find the weight for computing target value - depends on iteration
    %counter
    if mod(tv.iterctr,length(tv.sample_sequence))==0
        sample_stage=tv.sample_sequence(end);
    else
        sample_stage=tv.sample_sequence(mod(tv.iterctr,length(tv.sample_sequence)));
    end

    %see Holmstrom 2008 "An adaptive radial basis algorithm (ARBF) for expensive black-box global optimization", JOGO
    if sample_stage == 0 %InfStep - minimize Mu_n
        [PHI,P,phi0, pdim]=  rbf_matrices(Data,rbf_flag); %compute matrices for RBF model
         minimize_mu2 = @(x)inf_step(x, Data,PHI,P,phi0,pdim, rbf_flag);
        [x_ga, ~]=ga(minimize_mu2,Data.dim,[],[],[],[],Data.xlow,Data.xup,[],Data.integer, gaoptions);%use GA for optimization
        xselected = x_ga; %new sample point
    elseif 1 <= sample_stage && sample_stage <= tv.n %cycle step global search
        %find min of surrogate model
        minimize_RBF=@(x)rbf_prediction(x,Data,lambda, gamma, rbf_flag);
        [~, f_rbf]=ga(minimize_RBF,Data.dim,[],[],[],[],Data.xlow,Data.xup,[],Data.integer, gaoptions);
        wk = (1-sample_stage/tv.n)^2; % select weight for computing target value
        f_target = f_rbf -wk*(max(Data.Y) - f_rbf); %target for objective function value
        minimize_bump2 = @(x)bumpiness_measure(x,Data,f_target,  lambda, gamma,PHI,P,phi0,pdim, rbf_flag);
        %use GA to minimize bumpiness measure
        [x_bump, ~]=ga(minimize_bump2,Data.dim,[],[],[],[],Data.xlow,Data.xup,[],Data.integer, gaoptions);
        xselected=x_bump; %new point is the one that minimizes the bumpiness measure
    else %cycle step local search
        %find the minimum of RBF surface
        minimize_RBF=@(x)rbf_prediction(x,Data,lambda, gamma, rbf_flag);
        [x_rbf, f_rbf]=ga(minimize_RBF,Data.dim,[],[],[],[],Data.xlow,Data.xup,[],Data.integer, gaoptions);
        if f_rbf < min(Data.Y)-10^(-6)*abs(min(Data.Y))
            xselected = x_rbf; %select minimum point as new sample point if sufficient improvements
        else %otherwise, do target value strategy
            f_target = min(Data.Y)-10^(-2)*abs(min(Data.Y)); %target value
            minimize_bump2 = @(x)bumpiness_measure(x,Data,f_target, lambda, gamma,PHI,P,phi0,pdim,rbf_flag);
            [x_bump,~]=ga(minimize_bump2,Data.dim,[],[],[],[],Data.xlow,Data.xup,[],Data.integer, gaoptions);
            xselected=x_bump;
        end
    end
    [~,dist]=knnsearch(Data.S,xselected); %find closest already sampled point to xselected
    while dist < Data.tol %the selected point is too close to already evaluated point
        %randomly select point from variable domain
        xselected = Data.xlow +(Data.xup-Data.xlow).*rand(1,Data.dim);
        xselected(Data.integer) = round(xselected(Data.integer));
        [~,dist]=knnsearch(Data.S,xselected);
    end
    fevalt = tic;     %start timer for function evaluation
    fnew = feval(Data.objfunction,xselected); %new obj function value
    timer = toc(fevalt); %stop timer
    
    Data.m=Data.m+1; %update the number of function evaluations
    Data.S(Data.m,:)=xselected; %update sample site matrix
    Data.Y(Data.m,1)=fnew; %update vector with function values
    Data.T(Data.m,1) = timer; %record objective function evaluation time
    if fnew < Data.fbest %update best point found so far if necessary 
        Data.xbest = xselected; %best point found so far
        Data.fbest = fnew; %best objective function value found so far
    end
   
end %function



