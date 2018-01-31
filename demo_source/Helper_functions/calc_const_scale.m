function scale = calc_const_scale(f,dom)
%
%   Bevan Cheeseman 2017
%
%   This computes a constant scale, setting it tot he maximum value of the
%   function over the interval [dom(1),dom(2)];
%

inv_f = @(x) (1./abs(f(x)));
f_ = @(x) abs(f(x));

% find the minimum on the domain
x_max = fminbnd(inv_f,dom(1),dom(2));

brute_y = linspace(dom(2),dom(1),4000);
[bm,bi] = min(inv_f(brute_y));

if(bm < inv_f(x_max))
   x_max = brute_y(bi); 
end

x_min = fminbnd(f_,dom(1),dom(2));

scale = abs(f(x_max)-f(x_min));

end