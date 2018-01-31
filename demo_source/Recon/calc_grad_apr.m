function grad = calc_grad_apr(apr)
%
%   Bevan Cheeseman 2017
%
%   Calculate Gradient
%
%

grad = zeros(apr.num_parts,1);

grad(1) = (apr.f_p(2)-apr.f_p(1))/(apr.y_p(2)-apr.y_p(1));
grad(end) = (apr.f_p(end)-apr.f_p(end-1))/(apr.y_p(end)-apr.y_p(end-1));

for p = 2:(apr.num_parts-1)
    grad(p) = 0.5*(apr.f_p(p)-apr.f_p(p-1))/(apr.y_p(p)-apr.y_p(p-1))...
      +  0.5*(apr.f_p(p+1)-apr.f_p(p))/(apr.y_p(p+1)-apr.y_p(p));
end



end