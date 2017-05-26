function new = newp(x,M,tol)

% takes a point x and a matrix A, and computes their distance to check
% whether x is too close to any point in A
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%Input:
%x - the point to be checked
%M - matrix with already selected points
%tol - distance below which x is considered too close to A
%
%Output:
%output - 0 if x is too close to A; 1 otherwise
%pred_f - predicted objective function values corresponding to new_p
%--------------------------------------------------------------------------

[~,b] = knnsearch(x,M);
if any(b<tol)
    new = 0;
else
    new = 1;
end


end %function