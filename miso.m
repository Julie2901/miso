function [xbest, fbest] = miso(datafile,maxeval, surrogate, n_start, init_design, sampling, ...
    own_design)


% miso.m is a surrogate model based optimization algorithm that aims at
% finding a (near) optimal solution of computationally expensive, black-box, 
% global optimization problems with mixed-integer constraints.
% When using MISO to solve your optimization problem, please cite the paper
% J. Mueller: 'MISO: mixed-integer surrogate optimization framework', 2015,
% to appear in Optimization and Engineering, DOI: 10.1007/s11081-015-9281-2
% This MATLAB code comes with a documentation on how to use it. 

% MISO input files:
% datafile -- string, user defined problem
% maxeval -- integer, max number of allowed function evaluations
% surrogate -- string with surrogate model type
% n_start -- integer, number of points for initial design
% init_design -- string, initial design strategy
% sampling -- string, sampling strategy
% own_design -- matrix (m x dimension) with initial points

% MISO output files:
% xbest - best point found during optimization
% fbest - best function value found during optimization

%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------




global sampledata %global variable to collect sample data
sampledata = []; %initialize sampledata to empty matrix 
    
%% check user inputs
if nargin <1 %user forgot to give file name with problem definition
    error('You must provide a data file with the problem definition')
end
%check how many input arguments are given and assign [ ] to values that
%were not supplied
if nargin ==1
    maxeval =[];
    surrogate=[];
    n_start=[];
    init_design=[];
    sampling=[];
    own_design=[];
elseif nargin == 2
    surrogate=[];
    n_start=[];
    init_design=[];
    sampling=[];
    own_design=[];
elseif nargin == 3
    n_start=[];
    init_design=[];
    sampling=[];
    own_design=[];
elseif nargin == 4
    init_design=[];
    sampling=[];
    own_design=[];
elseif nargin == 5
    sampling=[];
    own_design=[];
elseif nargin == 6
    own_design=[];
end
    
    
%load optimization problem data from user-defined file
Data=feval(datafile); % load problem data

%check if max number of function evaluations is defined by user 
if isempty(maxeval) %use default max number of allowed function evaluations
    disp('Using default number of max. evals (50*dimension)!')
    Data.maxeval = 50*Data.dim;
elseif maxeval <Data.dim + 1 %number of allowable evaluations is too low to even fit model
    error('The user given maximum number of allowed function evaluations is too low to fit an RBF surrogate')
else
    Data.maxeval = maxeval;
end

%check if surrogate model defined
if isempty(surrogate) %user did not define whcih RBF model type is wanted, use default
    disp('Using default surrogate model (cubic RBF)!')
    Data.surrogate = 'rbf_c';
else
    Data.surrogate = surrogate;
end

%check if the number of initial design points is defined
if isempty(n_start) && isempty(own_design) %number of start points not defined, no user-given design 
    disp('Using default number of initial design points (2*(dimension+1))!')
    n_start = 2*(Data.dim+1); %default initial design size
else
    if ~isempty(n_start) && n_start<Data.dim+1 && isempty(own_design)
        disp('The number of deisgn points must be at least (dimension + 1). Using (dimension +1) design points!')
        n_start = Data.dim + 1;
    end
end

%check if the type of initial experimental design is defined
fixrank = false;
if isempty(init_design) %no initial design strategy defined
    disp('Using the default initial experimental design strategy (symmetric Latin hypercube)!')
    Data.init_design = 'slhd';
elseif strcmp(init_design, 'own') %user has to supply a matrix with selected sample points to be included in initial design
    Data.init_design = 'own';
    if isempty(own_design)
        error('If you use init_design = "own", you must provide a matrix (own_design) with initial points as input!')
    else
        Data.own_design = own_design;
    end
    rod = own_design; %check if user-given design is integer feasible
    rod(:,Data.integer) = round(own_design(:,Data.integer));
    if any(sum((own_design-rod).^2,2)>10e-6) %initial design did not satisfy integer constraints
        error('The user given initial design does not satisfy the integer constraints!')
    end
    if size(own_design,2) < Data.dim %check if given points are of the correct dimension
        error('The matrix with initial design points has the wrong dimension (number of columns must equal problem dimension)!')
    elseif (size(own_design,1) < Data.dim+1 && isempty(n_start)) || (~isempty(n_start) && n_start +size(own_design,1) < Data.dim+1)
        disp('User given initial design does not have sufficiently many points to fit RBF. Filling in remaining points by slhd!')
        n_start = Data.dim + 1 - size(own_design,1);
    else %user given design contains sufficiently many points
        %check rank of matrix
        if rank([own_design, ones(size(own_design,1),1)]) < Data.dim + 1
            disp('Rank of user given initial design matrix is too low. Generating more points to satisfy rank condition')
            n_start = Data.dim-rank(own_design);
        else
            n_start = 0;
        end
    end
elseif strcmp(init_design, 'slhd')
    Data.init_design = 'slhd';
elseif strcmp(init_design, 'lhs')
    Data.init_design = 'lhs';    
end

%check what sampling strategy is required
if isempty(sampling)
    disp('Using the default sampling strategy (cptvl)!')
    Data.sampling = 'cptvl';
else
    Data.sampling = sampling;
end

%set random number seed according to taskID
%rand('state',1); %set random number seed 
%randn('state',1);

TimeStart=tic;%record total computation time
Data.number_startpoints=n_start;
Data.tol = 0.001*min(Data.xup-Data.xlow); %tolerance when 2 points are considered the same

%% initial experimental design
if Data.number_startpoints > 0  
    if strcmp(Data.init_design, 'slhd') %symmetric latin hypercube
        InitialPoints = slhd(Data); 
    elseif strcmp(Data.init_design, 'lhs') %matlab's lhsdesign
        InitialPoints = lhsdesign(n_start,Data.dim);
    elseif strcmp(Data.init_design,'own')%use slhd in case own design does not have sufficiently many points
        InitialPoints = slhd(Data); 
    else
        error('Undefined initial experimental design strategy')
    end
    %scale S to true ranges
    Data.S=repmat(Data.xlow, Data.number_startpoints,1) + repmat(Data.xup-Data.xlow,Data.number_startpoints,1).*InitialPoints; 
else
    Data.S = [];
end
%check if user gave (partial) initial design
if ~isempty(own_design)
    Data.S = [own_design; Data.S];
    Data.number_startpoints =size(Data.S,1); 
end
%round integer variables to integers
Data.S(:,Data.integer)=round(Data.S(:,Data.integer));
% check if rank of the initial sample size matrix is large enough
if rank([Data.S,ones(size(Data.S,1),1)]) < Data.dim + 1
    fixrank = true;
end 
while fixrank %rank is too small to fit RBF, add random points to initial design to satisfy rank condition
    n_new = Data.dim+1-rank([Data.S,ones(size(Data.S,1),1)]); %minimum number of new points needed
    randpoint = repmat(Data.xlow,n_new,1) + repmat((Data.xup-Data.xlow), n_new,1).*rand(n_new,Data.dim);
    randpoint(:,Data.integer) = round(randpoint(:,Data.integer));
    if rank([[Data.S;randpoint], ones(size(Data.S,1)+n_new,1)]) == Data.dim + 1
        Data.S = [Data.S; randpoint];
        fixrank = false;
    end 
end
    
%% do expensive evaluations 
Data.m=size(Data.S,1);
Data.Y=zeros(Data.m,1);
Data.T = zeros(Data.m,1);
for ii = 1: Data.m
    fevalt = tic;
    Data.Y(ii)=feval(Data.objfunction, Data.S(ii,:));
    Data.T(ii) = toc(fevalt);
end
[Data.fbest,fbest_ID]=min(Data.Y); %best function value so far
Data.xbest=Data.S(fbest_ID,:); %best point so far    
    
  
%% sampling
if strcmp(Data.sampling, 'cp') %candidates by perturbation of randomly selected variables
    sol =  cp(Data);
elseif strcmp(Data.sampling, 'tv') %target value strategy
    sol = tv(Data);
elseif strcmp(Data.sampling,'ms') %surface minimum
    sol = ms(Data);
elseif strcmp(Data.sampling, 'rs') %random candidates by perturbing all variables
    sol = rs(Data);
elseif strcmp(Data.sampling, 'cptv') % cp+tv
    sol = cptv(Data);
elseif strcmp(Data.sampling,'cptvl') % cp + tv + fmincon on continuous variables
    sol = cptvl(Data);
else
   error('Sampling procedure not included'); 
end
xbest = sol.xbest;
fbest = sol.fbest;
sol.total_T = toc(TimeStart); %stop timer for optimization

%% saving
save results.mat sol %save results

%% plot progress curve
%yvals contains the best function value found so far
yvals = zeros(sol.m,1);
yvals(1) = sol.Y(1);
for ii = 2: sol.m
   if yvals(ii-1) < sol.Y(ii)
       yvals(ii) = yvals(ii-1);
   else
        yvals(ii) = sol.Y(ii);
   end
end

figure
plot((1:sol.m),yvals);
xlabel('Number of function evaluations')
ylabel('Objective function value')
title('Progress plot')

end
    