function [R_i,y_i] = calc_implied_res_func_doubled(i,l,dom)
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
    
    R_i(2*(n) - 1) = size*2;
    R_i(2*(n)) = size*2;
    
     
    y_i(2*(n) - 1) = loc;
    y_i(2*(n)) = loc + size;
    
    
end







end