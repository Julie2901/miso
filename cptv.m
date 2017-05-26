function Data =  cptv(Data)


% sampling procedure that switches between coordinate perturbation 
% approach (cp.m) and target value approach (tv.m)
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

%% parameters and initialization for coordinate perturbation
cp.ncand= min(100*Data.dim,5000); % number of candidate points
cp.xrange = Data.xup - Data.xlow; %variable range in every dimension
cp.minxrange = min(cp.xrange); %smallest variable range
cp.sigma_stdev_default = 0.2*cp.minxrange;  %default (initital) perturbation range
cp.sigma_stdev = cp.sigma_stdev_default;   % current perturbation rnage
cp.sigma_min=0.2*(1/2)^6*cp.minxrange; %minimum perturbation range
cp.maxshrinkparam = 5;% maximal number of perturbation range reductions 
cp.failtolerance = max(5,Data.dim); %threshold for the number of consecutively failed improvement trials
cp.succtolerance =3; %threshold for consecutively successful improvement trials
cp.iterctr = 0; % number of iterations
cp.shrinkctr = 0; % number of times perturbation range was reduced
cp.failctr = 0; % number of consecutive unsuccessful iterations
cp.localminflag = 0;  % indicates whether or not xbest is at a local minimum
cp.succctr=0; % number of consecutive successful iterations
cp.weightpattern=[0.3,0.5,0.8,0.95]; %weight pattern for computing scores of candidate points

%% parameters and initialization for target value search
tv.iterctr = 0; %initialize target value sampling iteration counter
tv.failctr = 0; % number of consecutive unsuccessful iterations
tv.succctr=0; %number of consecutive successful iterations
tv.failtolerance = max(5,Data.dim); %threshold for the number of consecutively failed improvement trials
tv.w_j=[0, 0.50,  2, 100]; %weight factors for computing target values
tv.w_counter=0;
tv.n=10; %number of cycles in global search
tv.sample_sequence=[0,(1:tv.n),tv.n+1]; %sample sequence: inf-step, cycle step gobal, cycle step local
tv.failtolerance = length(tv.sample_sequence); %threshold for the number of failed improvement trials
gaoptions=gaoptimset('PopulationSize',10,'Display','off'); %options for genetic algorithm

cp_search=true; %start search with coordinate perturbation 

%% optimization steps
while Data.m < Data.maxeval %do until budget of function evaluations exhausted 
    %compute parameters of the surrogate model
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
 
    fprintf('\n Number of function evaluations: %d \n', Data.m);
    fprintf('\n Best function value so far: %d \n', Data.fbest);
    if cp_search %search with coordinate perturbation
        cp.iterctr = cp.iterctr + 1; %increment iteration counter
        %select weight for score computation
        mw=mod(cp.iterctr,length(cp.weightpattern));
        if mw==0
            w_r=cp.weightpattern(end);
        else
            w_r=cp.weightpattern(mw);
        end
        
        % Perturbation probability
        pert_p = min(20/Data.dim,1)*(1-(log(Data.m-2*(Data.dim+1)+1)/log(Data.maxeval-2*(Data.dim+1))));
        %create candidate points
        CandPoint = repmat(Data.xbest,cp.ncand,1); 
        for ii =1:cp.ncand
            r=rand(1,Data.dim);
            ar = r<pert_p;
            if ~(any(ar))
                r = randperm(Data.dim);
                ar(r(1))=1;
            end
            for jj =1:Data.dim
                if ar(jj)==1
                    if ismember(jj,Data.integer) %integer perturbation has to be at least 1 unit 
                        rr=randn(1);
                        s_std=sign(rr)*max(1,abs(cp.sigma_stdev*rr));
                    else
                        s_std= cp.sigma_stdev*randn(1);
                    end
                    CandPoint(ii,jj) = CandPoint(ii,jj) +s_std;
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
        end
        CandPoint(:,Data.integer)=round(CandPoint(:,Data.integer)); %round integer values of candidate points
        xnew = compute_scores(Data,CandPoint, w_r, lambda, gamma, rbf_flag); %score candidates and select the best candidate
        clear CandPoint;
        fevalt = tic; %start timer for function evaluation
        fnew = feval(Data.objfunction,xnew); %expensive function value
        timer = toc(fevalt); %stop timer for function evaluation
        Data.m=Data.m+1; %update the number of function evaluations
        Data.S(Data.m,:)=xnew; %update the sample site matrix with newly sampled point
        Data.Y(Data.m,1)=fnew; %update vector with function values
        Data.T(Data.m,1)= timer; %update vector with function evaluation time
        
        if fnew < Data.fbest %update best point found so far if necessary
            if (Data.fbest - fnew) > (1e-3)*abs(Data.fbest)
                % "significant" improvement
                cp.failctr = 0;
                cp.succctr=cp.succctr+1; %update counter for successful trials
            else
                %no "significant" improvement
                cp.failctr = cp.failctr + 1; %update counter for unsuccessful trials
                cp.succctr=0;
            end  
            Data.xbest = xnew; %best point found so far
            Data.fbest = fnew; %best objective function value found so far
        else
            cp.failctr = cp.failctr + 1; %update counter for unsuccessful trials
            cp.succctr=0;
        end
        % check if algorithm is in a local minimum
        cp.shrinkflag = 1;      
        if cp.failctr >= cp.failtolerance %number of consecutive failed improvement trials is at threshold
            if cp.shrinkctr >= cp.maxshrinkparam %number of allowed perturbation range reductions is exceeded
                cp.shrinkflag = 0;
                disp('Stopped reducing perturbation range because the maximum number of reductions has been reached.');
            end
            cp.failctr = 0;

            if cp.shrinkflag == 1 % reduce perturbation range
                cp.shrinkctr = cp.shrinkctr + 1;
                cp.sigma_stdev =  max(cp.sigma_min,cp.sigma_stdev/2);
                disp('Reducing perturbation range to half!');
            else  
                cp.localminflag = 1; %local min encountered
                cp_search=false; %end coordinate perturbation search
                cp.shrinkctr=0; %reset counters
                cp.failctr = 0; 
                cp.succctr=0;
                cp.sigma_stdev =cp.sigma_stdev_default; %reset perturbation range 
            end
        end
    else %search with target value strategy
        tv.iterctr=tv.iterctr+1; %increase target value iteration counter
        %check which sample stage we are in
        if mod(tv.iterctr,length(tv.sample_sequence))==0
            sample_stage=tv.sample_sequence(end);
        else
            sample_stage=tv.sample_sequence(mod(tv.iterctr,length(tv.sample_sequence)));
        end

        %see Holmstrom 2008 "An adaptive radial basis algorithm (ARBF) for expensive black-box global optimization", JOGO
        if sample_stage == 0 %InfStep -minimize Mu_n
            [PHI,P,phi0, pdim]=  rbf_matrices(Data,rbf_flag); %set up matrices for RBF model
             minimize_mu2 = @(x)inf_step(x, Data,PHI,P,phi0,pdim, rbf_flag); %objective function
            [x_ga, ~]=ga(minimize_mu2,Data.dim,[],[],[],[],Data.xlow,Data.xup,[],Data.integer, gaoptions); %use GA to minimize      
            xselected = x_ga; %new sample point
        elseif 1 <= sample_stage && sample_stage <= tv.n %cycle step global search
            %find min of surrogate model
            minimize_RBF=@(x)rbf_prediction(x,Data,lambda, gamma, rbf_flag);
            [~, f_rbf]=ga(minimize_RBF,Data.dim,[],[],[],[],Data.xlow,Data.xup,[],Data.integer, gaoptions);
            wk = (1-sample_stage/tv.n)^2; %compute weight for target value computation
            f_target = f_rbf -wk*(max(Data.Y) - f_rbf); %target for objective function value  
            minimize_bump2 = @(x)bumpiness_measure(x,Data,f_target,  lambda, gamma,PHI,P,phi0,pdim, rbf_flag);
            [x_bump,~]=ga(minimize_bump2,Data.dim,[],[],[],[],Data.xlow,Data.xup,[],Data.integer, gaoptions); %use GA to minimize the bumpiness measure
            xselected=x_bump; %new point is the one that minimizes the bumpiness measure
        else
            minimize_RBF=@(x)rbf_prediction(x,Data,lambda, gamma, rbf_flag); %minimize the objective function
            [x_rbf, f_rbf]=ga(minimize_RBF,Data.dim,[],[],[],[],Data.xlow,Data.xup,[],Data.integer, gaoptions); %use GA to find min of surface
            if f_rbf < min(Data.Y)-10^(-6)*abs(min(Data.Y))  %min of surface is a significant improvement
                xselected = x_rbf; %new point is min of RBF surface
            else
                f_target = min(Data.Y)-10^(-2)*abs(min(Data.Y)); %define target value
                minimize_bump2 = @(x)bumpiness_measure(x,Data,f_target, lambda, gamma,PHI,P,phi0,pdim,rbf_flag);
                [x_bump,~]=ga(minimize_bump2,Data.dim,[],[],[],[],Data.xlow,Data.xup,[],Data.integer, gaoptions);%use GA to minimize the bumpiness measure
                xselected=x_bump; %new point is at min of bumpiness measure
            end
        end
        [~,dist]=knnsearch(Data.S,xselected); %find closest already sampled point to xselected
        while dist < Data.tol %the selected point is too close to already evaluated point
            %randomly select point from the variable domain
            xselected = Data.xlow +(Data.xup-Data.xlow).*rand(1,Data.dim);
            xselected(Data.integer) = round(xselected(Data.integer));
            [~,dist]=knnsearch(Data.S,xselected);
        end
        fevalt = tic; %start timer for objective function evaluation
        fnew = feval(Data.objfunction,xselected); %new obj function value
        timer = toc(fevalt); %stop timer for function evaluations
        Data.m=Data.m+1; %update the number of function evaluations
        Data.S(Data.m,:)=xselected; %update the sample size matrix
        Data.Y(Data.m,1)=fnew; %update vector with function values
        Data.T(Data.m,1) = timer; %update vector with evaluation times
        
        if fnew < Data.fbest %update best point found so far if necessary
            if (Data.fbest - fnew) > (1e-3)*abs(Data.fbest)
                % "significant" improvement
                tv.failctr = 0;
                tv.succctr=tv.succctr+1;%increment success counter
            else
                %no "significant" improvement
                tv.failctr = tv.failctr + 1; %increment fail counter
                tv.succctr=0;
            end  
            Data.xbest = xselected; %best point found so far
            Data.fbest = fnew; %best objective function value found so far
        else
            tv.failctr = tv.failctr + 1; % increment fail counter
            tv.succctr=0;
        end
        % check if algorithm is in a local minimum
        if tv.failctr >= tv.failtolerance 
            tv.localminflag = 1;
            cp_search=true; %back to coordinate perturbation search
            tv.failctr=0; %reset fail counter
            tv.succctr=0; %reset success counter
        end
    end    
end


end%function