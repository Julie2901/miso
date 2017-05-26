function Data=datainput_hartman3
 
% example optimization problem (computationally cheap)
% 3-dimensional Hartman function
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%Input: none
%
%Output:
%Data - structure array with all problem information
%--------------------------------------------------------------------------

Data.xlow =zeros(1,3);  % variable lower bounds
Data.xup=ones(1,3);     % variable upper bounds
Data.objfunction=@(x)hartman3(x); %handle to objective function
Data.dim = 3;   %problem dimesnion
Data.integer = [1];   %indices of integer variables
Data.continuous = [2,3]; %indices of continuous variables

end %function

function y=hartman3(x) %objective function
global sampledata; %global variable that collects sample points, function values and evaluation times
c=[1,1.2,3,3.2]'; %c,A,b are data vectors
A=[3, 10, 30; 0.1,10,35; 3, 10, 30;0.1,10,35];
P=[0.3689    0.1170    0.2673
    0.4699    0.4387    0.7470
    0.1091    0.8732    0.5547
    0.0382    0.5743    0.8828];
x=x(:)'; % make sure vector is row vector
fevalt = tic; %start timer for function evaluation
y=-sum(c.*exp(-sum(A.*(repmat(x,4,1)-P).^2,2))); %compute objective function value
t = toc(fevalt); %stop timer for function evaluation
sampledata = [sampledata; x(:)',y, t]; %collect sample data (point x, value y, time t)
end %hartman3