function [i_n,l_n] = get_neigh(i,l)
%   
%   Bevan Cheeseman 2017
%
%   Input a particle cell and get a list of the neighbours, including
%   current cell
%

i_max = min(i+1,2^l-1);
i_min = max(0,i-1);

i_n = i_min:i_max;
l_n = ones(size(i_n))*l;


end