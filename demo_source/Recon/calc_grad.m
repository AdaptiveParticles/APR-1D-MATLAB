function grad = calc_grad(y,f)
%
%   Bevan Cheeseman 2017
%
%   Calculate Gradient
%
%

grad = zeros(size(f));

grad(1) = (f(2)-f(1))/(y(2)-y(1));
grad(end) = (f(end)-f(end-1))/(y(end)-y(end-1));

for p = 2:(length(f)-1)
    grad(p) = 0.5*(f(p)-f(p-1))/(y(p)-y(p-1))...
      +  0.5*(f(p+1)-f(p))/(y(p+1)-y(p));
end


end