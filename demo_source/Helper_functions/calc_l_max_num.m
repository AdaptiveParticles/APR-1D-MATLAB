function l_max = calc_l_max_num(grad,dom,E,scale,order)
%
%   Bevan Cheeseman 2017
%
%   Calculates the maximum level l_max for sampling
%


%domain size
Sigma = (dom(2) - dom(1));

% 1/|\nabla f|
inv_grad_max = min(Sigma,1./abs(grad));


% find the minimum on the domain
f_min = min(inv_grad_max);


% this then gives us our minimum value of L
L_min = ((E*scale)*f_min).^(1/order);

%then find the l_max that rounds to that resolution
l_max = ceil(log2(Sigma/L_min));

end