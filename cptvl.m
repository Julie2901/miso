function Data =  cptvl(Data)

% sampling procedure that switches between coordinate perturbation 
% approach (cp.m), target value approach (tv.m), and a local minimization
% on the continuous variables using fmincon.m
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

global sampledata %global variable

%% parameters and initialization for coordinate perturbation
cp.ncand= min(100*Data.dim,5000); % number of candidate points
cp.xrange = Data.xup - Data.xlow; %variable range in every dimension
cp.minxrange = min(cp.xrange); %smallest variable range
cp.sigma_stdev_default = 0.2*cp.minxrange; %default perturbation range    
cp.sigma_stdev = cp.sigma_stdev_default;  % current perturbation range
cp.sigma_min=0.2*(1/2)^6*cp.minxrange; %minimum perturbation range
cp.maxshrinkparam = 5;% maximum number of times we can reduce the perturbation range 
cp.failtolerance = max(5,Data.dim); %threshold for the number of consecutively failed improvement trials
cp.succtolerance =3; %threshold for the number of consecutive successful improvement trials 
cp.iterctr = 0; % counter of iterations in coordinate perturbation search
cp.shrinkctr = 0; % counter of perturbation range reductions
cp.failctr = 0; % countert of consecutive unsuccessful iterations
cp.localminflag = 0;  % indicates whether or not xbest is at a local minimum
cp.succctr=0; % counter of consecutive successful iterations
cp.weightpattern=[0.3,0.5,0.8,0.95]; %weight pattern for computing candidate point scores

%% parameters and initialization for target value search
tv.iterctr = 0; %counter for target value iterations
tv.failctr = 0; % counter of consecutive unsuccessful iterations
tv.succctr=0; %counter for consecutive successful trials
tv.failtolerance = max(5,Data.dim); %threshold for the number of consecutively failed improvement trials
tv.w_j=[0, 0.50,  2, 100]; %weight factors used for computing target values
tv.w_counter=0; %counter for going through the weight factors for determining which weight to use 
tv.n=10; %length of cycle step global search
tv.sample_sequence=[0,(1:tv.n),tv.n+1]; %inf step, cycle step global, cycle step local
tv.failtolerance = length(tv.sample_sequence); %threshold for the number of consecutively failed improvement trials
gaoptions=gaoptimset('PopulationSize',10,'Display','off'); %options for genetic algorithm

cp_search=true; %start search by candidate point perturbation
tv_search=false;

cp_sw=[]; %initialize arrays for recording switch between search modes
tv_sw=[];

while Data.m < Data.maxeval %do until budget of function evaluations exhausted 
    %compute the parameters of the RBF model
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
        CandPoint(:,Data.integer)=round(CandPoint(:,Data.integer)); %round integer variables of the candidate points
        xnew = compute_scores(Data,CandPoint, w_r, lambda, gamma, rbf_flag); %select the best candidate based on scores
        clear CandPoint;
        fevalt = tic; %start timer for function evaluation
        fnew = feval(Data.objfunction,xnew); %new function evaluation
        timer = toc(fevalt); %stop timer for function evaluation
        Data.m=Data.m+1; %update the number of function evaluations
        Data.S(Data.m,:)=xnew; %update sample site matrix with newly sampled point
        Data.Y(Data.m,1)=fnew; %update function value vector
        Data.T(Data.m,1) = timer; %update evaluation time vector
        
        if fnew < Data.fbest %update best point found so far if necessary
            if (Data.fbest - fnew) > (1e-3)*abs(Data.fbest)
                % "significant" improvement
                cp.failctr = 0;
                cp.succctr=cp.succctr+1; %update success counter
            else
                %no "significant" improvement
                cp.failctr = cp.failctr + 1; %update fail counter
                cp.succctr=0;
            end  
            Data.xbest = xnew; %best point found so far
            Data.fbest = fnew; %best objective function value found so far
        else
            cp.failctr = cp.failctr + 1; %update fail counter
            cp.succctr=0;
        end
        % check if algorithm is in a local minimum
        cp.shrinkflag = 1;      
        if cp.failctr >= cp.failtolerance %number of consecutively failed improvement trials is at threshold
            if cp.shrinkctr >= cp.maxshrinkparam %check if perturbation range can be reduced
                cp.shrinkflag = 0;
                disp('Stopped reducing perturbation range because the maximum number of allowed reductions has been reached.');
            end
            cp.failctr = 0;

            if cp.shrinkflag == 1 %reduce perturbation range by half
                cp.shrinkctr = cp.shrinkctr + 1; %update counter for perturbation range reductions
                cp.sigma_stdev =  max(cp.sigma_min,cp.sigma_stdev/2);
                disp('Reducing perturbation range to half!');
            else  
                cp.localminflag = 1;
                cp_search=false; %local min encountered, switch to target value strategy
                cp.shrinkctr=0; %reset counters
                cp.failctr = 0;
                cp.succctr=0;
                cp.sigma_stdev =cp.sigma_stdev_default;
                cp_sw=[cp_sw;Data.m, Data.fbest]; %record the number of evaluations when sampling methods switched and best point so far
                if (size(cp_sw,1)==1) || (size(cp_sw,1)>1 &&  cp_sw(end,2) < cp_sw(end-1,2))
                    tv_search=true;%go to target value search
                else
                    tv_search=false;
                end
            end
        end
    elseif tv_search %search with target value strategy
        tv.iterctr=tv.iterctr+1; %update iteration counter for target value strategy
        %find the weight for computing target value
        if mod(tv.iterctr,length(tv.sample_sequence))==0
            sample_stage=tv.sample_sequence(end);
        else
            sample_stage=tv.sample_sequence(mod(tv.iterctr,length(tv.sample_sequence)));
        end

        %see Holmstrom 2008 "An adaptive radial basis algorithm (ARBF) for expensive black-box global optimization", JOGO
        if sample_stage == 0 %InfStep -minimize Mu_n
            [PHI,P,phi0, pdim]=  rbf_matrices(Data,rbf_flag); %set up matrices for computing RBF parameters
             minimize_mu2 = @(x)inf_step(x, Data,PHI,P,phi0,pdim, rbf_flag);
            [x_ga, ~]=ga(minimize_mu2,Data.dim,[],[],[],[],Data.xlow,Data.xup,[],Data.integer, gaoptions);     
            xselected = x_ga; %new sample point
        elseif 1 <= sample_stage && sample_stage <= tv.n %cycle step global search
            %find min of surrogate model
            minimize_RBF=@(x)rbf_prediction(x,Data,lambda, gamma, rbf_flag);
            [~, f_rbf]=ga(minimize_RBF,Data.dim,[],[],[],[],Data.xlow,Data.xup,[],Data.integer, gaoptions);
            wk = (1-sample_stage/tv.n)^2; %compute weight for determining target value
            f_target = f_rbf -wk*(max(Data.Y) - f_rbf); %target for objective function value          
            minimize_bump2 = @(x)bumpiness_measure(x,Data,f_target,  lambda, gamma,PHI,P,phi0,pdim, rbf_flag);
            [x_bump,~]=ga(minimize_bump2,Data.dim,[],[],[],[],Data.xlow,Data.xup,[],Data.integer, gaoptions);
            xselected=x_bump; %new point is the one that minimizes the bumpiness measure
        else %cycle step local search
            minimize_RBF=@(x)rbf_prediction(x,Data,lambda, gamma, rbf_flag); %find min of RBF surface
            [x_rbf, f_rbf]=ga(minimize_RBF,Data.dim,[],[],[],[],Data.xlow,Data.xup,[],Data.integer, gaoptions);
            if f_rbf < min(Data.Y)-10^(-6)*abs(min(Data.Y)) %min of RBF surface is sufficient improvement
                xselected = x_rbf; %new point is min point of RBF surface
            else
                f_target = min(Data.Y)-10^(-2)*abs(min(Data.Y)); %define target value and minimize bumpiness
                minimize_bump2 = @(x)bumpiness_measure(x,Data,f_target, lambda, gamma,PHI,P,phi0,pdim,rbf_flag);
                [x_bump,~]=ga(minimize_bump2,Data.dim,[],[],[],[],Data.xlow,Data.xup,[],Data.integer, gaoptions);
                xselected=x_bump; %new point is the one that minimizes bumpiness
            end
        end
        [~,dist]=knnsearch(Data.S,xselected); %find closest already sampled point to xselected
        while dist < Data.tol %the newly selected point is too close to already evaluated points
            %randomly select point from variable domain
            xselected = Data.xlow +(Data.xup-Data.xlow).*rand(1,Data.dim);
            xselected(Data.integer) = round(xselected(Data.integer));
            [~,dist]=knnsearch(Data.S,xselected);
        end
        fevalt = tic; %start timer for function evaluation
        fnew = feval(Data.objfunction,xselected); %new obj function value
        timer = toc(fevalt); %stop timer for function evaluation
        Data.m=Data.m+1; %update the number of function evaluations
        Data.S(Data.m,:)=xselected; %update sample site matrix
        Data.Y(Data.m,1)=fnew; %update vector with function values
        Data.T(Data.m,1) = timer; %update vector with evaluation times
        if fnew < Data.fbest %update best point found so far if necessary
            if (Data.fbest - fnew) > (1e-3)*abs(Data.fbest)
                % "significant" improvement
                tv.failctr = 0;
                tv.succctr=tv.succctr+1; %update success counter
            else
                %no "significant" improvement
                tv.failctr = tv.failctr + 1; %update fail counter
                tv.succctr=0;
            end  
            Data.xbest = xselected; %best point found so far
            Data.fbest = fnew; %best objective function value found so far
        else
            tv.failctr = tv.failctr + 1; %update fail counter
            tv.succctr=0;
        end
        % check if algorithm is in a local minimum
        if tv.failctr >= tv.failtolerance %number of consecutively failed improvement trials is at threshold
            tv.localminflag = 1; %local min found
            cp_search=true; %back to coordinate perturbation search
            tv.failctr=0; %reset counters
            tv.succctr=0;           
            tv_sw=[tv_sw;Data.m, Data.fbest]; %record number of evals until search strategy switch and best function value 
        end
    else %local search
        A=zeros(length(Data.integer), Data.dim); %initialize matrix A and vector b for linear constraints that force integer variables to not change
        b=zeros(length(Data.integer),1);
        for ii =1:length(Data.integer)
            A(ii,ii) = 1;
            b(ii) = Data.xbest(ii);
        end
        fmc_options = optimset('Display','off', 'MaxFunEvals',Data.maxeval-length(Data.Y)); %set fmincon options
        %use fmincon to minimize the true objective function keeping the
        %integer variables fixed (Ax=b constraints)
        [x_val,f_val] = fmincon(Data.objfunction,Data.xbest,[],[],A,b,Data.xlow,Data.xup,[],fmc_options); 
        
        if f_val <Data.fbest %update best point found
            Data.fbest=f_val;
            Data.xbest=x_val;
        end
        %save all sampledata elements as points in Data.S (columns
        %1:dimension); Data.Y is column dimension +1; evaluation times is
        %column dimension+2
        Data.S=sampledata(1:min(Data.maxeval,size(sampledata,1)),1:Data.dim);
        Data.Y=sampledata(1:min(Data.maxeval,size(sampledata,1)),Data.dim+1);
        Data.T=sampledata(1:min(Data.maxeval,size(sampledata,1)),Data.dim+2);
        Data.m=length(Data.Y);
        cp_search =true; %go back to candidate point search
    end
end

end%function