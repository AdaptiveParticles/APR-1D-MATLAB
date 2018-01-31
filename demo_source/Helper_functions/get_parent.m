function [p_i,p_l] = get_parent(i,l)
%
%   Bevan Cheeseman 2017
%
%   Return the parent particle cell of particle cell c_{i,l}
%

p_l = l-1;
p_i = floor(i/2);


end