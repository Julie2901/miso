function Data = datainput_branin
% example optimization problem (computationally cheap)
% 2-dimensional Branin function
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

Data.xlow=[-5,0]; %variable lower bounds
Data.xup=[10,15]; %variable upper bounds
%objective function
Data.objfunction= @(x)branin(x); 
Data.integer=1; %indices of integer variables
Data.continuous=2; %indices of continuous variables
Data.dim = 2; %problem dimension

end %function

function y=branin(x) %objective function
global sampledata; %global variable that collects sample points, function values and evaluation times

x=x(:)'; % make sure vector is row vector
fevalt = tic; %start timer for function evaluation
y=(x(:,2)-5.1*x(:,1).^2./(4*pi^2)+5*x(:,1)./pi-6).^2 + 10*(1-1/(8*pi))*cos(x(:,1))+10; %compute objective function value
t = toc(fevalt); %stop timer for function evaluation
sampledata = [sampledata; x(:)',y, t]; %collect sample data (point x, value y, time t)
end %branin