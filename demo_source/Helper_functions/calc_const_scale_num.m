function scale = calc_const_scale_num(f)
%
%   Bevan Cheeseman 2017
%
%   This computes a constant scale, setting it tot he maximum value of the
%   function over the interval [dom(1),dom(2)];
%


f_max = max(abs(f));

f_min = min(abs(f));

scale = abs(f_max-f_min);

end