function Data= datainput_rastrigin12
% example optimization problem (computationally cheap)
% 12-dimensional Rastrigin function
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
Data.xlow=-ones(1,12); %variable lower bounds
Data.xup=3*ones(1,12); %variable upper bounds
Data.dim=12; %problem dimension
Data.integer=(1:5); %indices of integer variables
Data.continuous=(6:12); %indices of continuous variables
Data.objfunction=@(x)rastrigin(x); %objective function handle
end %function

function y=rastrigin(x) %objective function
global sampledata; %global variable that collects sample points, function values and evaluation times

x=x(:)'; % make sure vector is row vector
fevalt = tic; %start timer for function evaluation
y=sum((x.^2) - cos(2*pi*x),2);  %compute objective function value
t = toc(fevalt); %stop timer for function evaluation
sampledata = [sampledata; x(:)',y, t]; %collect sample data (point x, value y, time t)
end %rastrigin