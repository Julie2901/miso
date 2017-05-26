function testdriver

% Run this code to test if your MATLAB knows the path to the MISO codes 
% and if toolboxes are installed.
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%Input: none
%
%Output: none
%--------------------------------------------------------------------------

disp('I am using MISO to solve a test problem. If you have all toolboxes installed, no errors should occur.')

[xbest, fbest] = miso('datainput_hartman3',100, 'rbf_t', [], 'lhs', 'cptv');
%[xbest, fbest] = miso('datainput_branin',100, 'rbf_t', [], 'slhd', 'cptv');
%[xbest, fbest] = miso('datainput_rastrigin12',100, 'rbf_t', [], 'slhd', 'cptv');
%[xbest, fbest] = miso('datainput_convex',100, 'rbf_t', [], 'slhd', 'cptv');
xbest
fbest




disp('Now type ''load results.mat'' into the command prompt to load the results.')

end %function