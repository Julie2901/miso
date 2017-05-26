function InitialPoints = slhd(Data)

% Create a symmetric Latin hypercube design 
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
%InitialPoints - points in the starting design
%--------------------------------------------------------------------------

delta = (1/Data.number_startpoints)*ones(1,Data.dim);

X = zeros(Data.number_startpoints,Data.dim);
for j = 1:Data.dim
    for i = 1:Data.number_startpoints
        X(i,j) = ((2*i-1)/2)*delta(j);
    end
end

P = zeros(Data.number_startpoints,Data.dim);
P(:,1) = (1:Data.number_startpoints)';
if (mod(Data.number_startpoints,2) == 0)
   k = Data.number_startpoints/2;
else
   k = (Data.number_startpoints-1)/2;
   P(k+1,:) = (k+1)*ones(1,Data.dim);
end

for j = 2:Data.dim
   P(1:k,j) = randperm(k)';
   for i = 1:k
      if (rand(1) <= 0.5)
         P(Data.number_startpoints+1-i,j) = Data.number_startpoints+1-P(i,j);
      else
         P(Data.number_startpoints+1-i,j) = P(i,j);
         P(i,j) = Data.number_startpoints+1-P(i,j);
      end
   end
end

InitialPoints = zeros(Data.number_startpoints,Data.dim);
for j = 1:Data.dim
    for i = 1:Data.number_startpoints
        InitialPoints(i,j) = X(P(i,j),j);
    end
end

end%function