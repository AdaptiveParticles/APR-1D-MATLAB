function [R_opt_b,y] =  compute_resolution_bound_cont(apr)
%
%   Bevan Cheeseman 2017
%
%   Brute force solve for the resolution function using iteration over all
%   points in the domain.
%
%

if(isfield(apr,'a'))

    a = apr.a;
else
    a = 1;
end

higher_res = 0;

%search sampling points, we have already found the lowest resolution required.
y = (apr.s_dom(1)):((apr.s_dom(2) - apr.s_dom(1))/2^(apr.l_max+higher_res)):apr.s_dom(2);

L = apr.f(y);

N = length(y);

scale = apr.scale(1);

E = apr.E;

dx =min(diff(y));

% see multiple resolution conditions, for the general bound
for i = 1:a
    apr.L{i} = ((apr.E*apr.scale(i))./abs(apr.df{i}(y))).^(1/i);
end

apr.L_f = apr.L{1};

% now take the minimum across them
for i = 2:a
   apr.L_f = min(apr.L_f,apr.L{i}); 
end


R_opt_b = zeros(size(y));

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