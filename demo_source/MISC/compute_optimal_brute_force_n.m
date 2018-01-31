function R_opt =  compute_optimal_brute_force_n(apr,y)
%
%   Bevan Cheeseman 2017
%
%   Brute force solve for the resolution function using iteration over all
%   points in the domain.
%
%

%search sampling points, we have already found the lowest resolution required.
f = apr.f(y);

N = length(y);

scale = apr.scale(1);

E = apr.E;

dx =min(diff(y));

R_opt = zeros(size(y));

for i = 1:N
    recon = 0;
    R = -1;
    
    while recon <= E
        R = R + 1;
        
        min_i = max(1,i-R);
        max_i = min(i+R,N);
        
        min_f = min(f(min_i:max_i));
        max_f = max(f(min_i:max_i));
        
        recon = max(abs(min_f - f(i)),abs(max_f - f(i)))/scale;
    end
  
    R_opt(i) = (R-1)*dx;
    
end



end