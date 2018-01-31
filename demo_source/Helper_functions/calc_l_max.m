function l_max = calc_l_max(grad,dom,E,scale,order)
%
%   Bevan Cheeseman 2017
%
%   Calculates the maximum level l_max for sampling
%


%domain size
Sigma = (dom(2) - dom(1));

% 1/|\nabla f|
inv_grad_max = @(x) min(Sigma,1./abs(grad(x)));



% find the minimum on the domain
x_max = fminbnd(inv_grad_max,dom(1),dom(2));

brute_y = linspace(dom(2),dom(1),4000);
[bm,bi] = min(inv_grad_max(brute_y));

if(bm < inv_grad_max(x_max))
   x_max = brute_y(bi); 
end

% this then gives us our minimum value of L
L_min = ((E*scale)*inv_grad_max(x_max)).^(1/order);


%then find the l_max that rounds to that resolution
l_max = ceil(log2(Sigma/L_min));

end