function Data= datainput_convex
 
% example optimization problem (computationally cheap)
% convex 8-dimensional test function
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

Data.xlow=-10*ones(1,8); %variable lower bounds
Data.xup=10*ones(1,8);  %variable upper bounds
Data.dim = 8; %problem dimension
Data.integer=(1:3); %indices of integer variables
Data.continuous=[4:8];%(5:8); %indices of continuous variables
Data.objfunction=@(x)convex_testfunction(x);
end %function

function y = convex_testfunction(x)
global sampledata;
fevalt = tic; 
y = 3.1*x(:,1).^2 + 7.6* x(:,2).^2 +6.9*x(:,3).^2 +0.004*x(:,4).^2 +...
    +19*x(:,5).^2 +3*x(:,6).^2 +x(:,7).^2  +4*x(:,8).^2 ; %objective function handle
t = toc(fevalt);
sampledata = [sampledata; x(:)',y, t];
end