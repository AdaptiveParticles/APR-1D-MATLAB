function [x_p,temp_p] = cont_sample_particles(R,y)
%
%   Bevan Cheeseman 2017
%
%
%   Brute force sampling scheme.
%
%

dx = min(diff(y));

%first re-scale the algorithm in terms of pixels
R_p = round(R/dx);

N = length(R);

%empty indicator array containing particles
temp_p = zeros(N,1);

temp_p(1) = 1;

for i = 1:N
   
    lb = max(i-R_p(i),1);
    
    inside = sum(temp_p(lb:i));
    
    if(inside > 0)
       %everything is okay still there is a particle 
    else
       %we need another particle
       temp_p(i) = 1;
    end
    
end

x_p = y((temp_p)==1);


end