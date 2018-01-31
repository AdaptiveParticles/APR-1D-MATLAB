function [l_c,i_c] = get_children(i,l)
%
%   Bevan Cheeseman 2017
%
%   Given a current paritcle cell c_{i,l}, compute its set of children
%
%

% two children in 1D
i_c = [i*2,i*2+1];
%ones level higher
l_c = ones(size(i_c))*(l+1);


end