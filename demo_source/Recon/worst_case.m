function [max_rc,min_rc,y] =worst_case(apr)
%
%   Bevan Cheeseman 2017
%
%   Constructs the worst case re-construction based on the apr   
%
%

Sigma = apr.s_dom(2) - apr.s_dom(1);

y = apr.s_dom(1):(Sigma/2^apr.l_max):apr.s_dom(2); 



l_ = interp1(apr.y_p,apr.c_l,y,'nearest');
nign = (~isnan(l_));

l_ = l_(nign);
y = y(nign);

min_rc = zeros(size(y));
max_rc = zeros(size(y));

for i = 1:length(y)
   R = Sigma/(2^l_(i)); 
   p = y(i);
   
   neigh_ind = find((apr.y_p >= (p - R)) &  (apr.y_p <= (p + R)));
   
   min_rc(i) = min(apr.f_p(neigh_ind));
   max_rc(i) = max(apr.f_p(neigh_ind));
   
end



end