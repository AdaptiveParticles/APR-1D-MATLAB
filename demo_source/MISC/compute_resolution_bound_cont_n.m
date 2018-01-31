function [R_opt_b,y] =  compute_resolution_bound_cont_n(apr,y)
%
%   Bevan Cheeseman 2017
%
%   Brute force solve for the resolution function using iteration over all
%   points in the domain.
%
%

dx =min(diff(y));

N = length(y);

% for i = 1:N
%     
%     R = -1;
%     bound = true;
%     
%     while bound
%         R = R + 1;
%         
%         min_i = max(1,i-R);
%         max_i = min(i+R,N);
%         
%         min_L = min(apr.L_f(min_i:max_i));
%         
%         bound = R*dx<=min_L;
%     end
%     
%     R_opt_b(i) = (R-1)*dx;
%     
% end

for i = 1:N
    
    R = -1;
    bound = true;
    
    min_L = 9999999;
    
    while bound
        R = R + 1;
        
        min_i = max(1,i-R);
        max_i = min(i+R,N);
        
        min_L = min(min(min_L,apr.L_f(min_i)),apr.L_f(max_i));
        
        bound = R*dx<=min_L;
    end
    
    R_opt_b(i) = (R-1)*dx;
    
end



end