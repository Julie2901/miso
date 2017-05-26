function Data = ms(Data)

% use the minimum point of the surrogate model as new sample point
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

ms.iterctr = 0; %initialize iteration counter

while Data.m < Data.maxeval %do until budget of function evaluations exhausted 
    %compute RBF parameters for chosen RBF type
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

    ms.iterctr=ms.iterctr+1;
    fprintf('\n Iteration: %d \n',ms.iterctr); %print some user information
    fprintf('\n Number of function evaluations: %d \n', Data.m);
    fprintf('\n Best function value so far: %d \n', Data.fbest);
    
    %use MLSL algorithm to find the local minima of the RBF surface
    [xnew, ~] = mlsl(Data,lambda, gamma, rbf_flag);
    % go through the suggested new points and discard all points that are
    % too close to previously sampled points
    kk = 1;
    while kk <= size(xnew,1)
        [~, Distan]=knnsearch(Data.S, xnew(kk,:));
        if Distan < Data.tol
            xnew(kk,:)=[];
        else
            kk=kk+1;
        end
    end
    if isempty(xnew) %if all suggested points are too close to already evaluated points, pick a random point from the variable domain
        xnew=Data.xlow +(Data.xup-Data.xlow).*rand(1,Data.dim);
        [~, Distan]=knnsearch(Data.S, xnew);
        while Distan <Data.tol
            xnew=Data.xlow +(Data.xup-Data.xlow).*rand(1,Data.dim);
            [~, Distan]=knnsearch(Data.S, xnew);
        end
    end
    
    for kk = 1:size(xnew,1) %go through the list of newly selected sample points and do expensive evaluations
        fevalt = tic; %start timer for function evaluation
        fnew = feval(Data.objfunction,xnew(kk,:)); %new function value
        timer = toc(fevalt); %stop timer
        Data.m=Data.m+1; %update the number of function evaluations
        Data.S(Data.m,:)=xnew(kk,:); %update sample site matrix with newly sampled point
        Data.Y(Data.m,1)=fnew; %update function value vector
        Data.T(Data.m,1) = timer; %update vector with evaluation times
        if fnew < Data.fbest %update best point found so far if necessary
            Data.xbest = xnew(kk,:); %best point found so far
            Data.fbest = fnew; %best objective function value found so far
        end
    
    end
    
end %function