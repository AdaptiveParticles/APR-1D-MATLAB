function R = calc_implied_res_func_y(i,l,dom,y)
%
%   Bevan Cheeseman 2017
%
%   Plots the original 
%
%

y_i = zeros(length(i)*2,1);
R_i = zeros(length(i)*2,1);

sig = (dom(2) - dom(1));
for n = 1:length(i)
    
    size = sig/2^l(n);
    loc = (i(n))*size + dom(1);
    
    R_i(2*(n) - 1) = size;
    R_i(2*(n)) = size;
    
     
    y_i(2*(n) - 1) = loc;
    y_i(2*(n)) = loc + .9999*size;
    
end

%y_i = [y(1);y_i;y(end)];
%R_i = [R_i(1);R_i;R_i(end)];

R = interp1(y_i,R_i,y);




end