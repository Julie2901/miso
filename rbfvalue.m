function p = rbfvalue(d, rbf_flag)

% compute the RBF value of a given distance d
%--------------------------------------------------------------------------
%Author information
%Juliane Mueller
%juliane.mueller2901@gmail.com
%--------------------------------------------------------------------------
%
%Input:
%d - distance between 2 points
%rbf_flag - indicator which rbf type should be used
%
%Output:
%p - radial basis function value phi(d) 
%--------------------------------------------------------------------------

%compute RBF value based on the type of RBF required
if strcmp(rbf_flag,'cub') %cubic RBF
    p = d.^3;
elseif strcmp(rbf_flag, 'lin') %inear RBF
    p = d;
elseif strcmp(rbf_flag,'tps')
    if d > 0 %thin-plate spline RBF
        p = d.^2.*log(d+realmin);
    else
        p = zeros(size(d));
    end
end

end %function