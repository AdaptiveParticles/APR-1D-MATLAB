function [l,i] = calculate_particle_cell(y,L,dom,l_max,varargin)
%
%   For a given y and L combo it calcualtes the particle cell level l and
%   spatial coordinate i on that level
%
if ~isempty(varargin)
   type = varargin{1}; 
else
    type = 1;
end

%domain size
Sigma = max(dom(2) - dom(1));
% level
l = ceil(log2(Sigma./L));

if type >1
   l = l - 1;
end

l = max(1,l);

l = min(l,l_max);

y = y(:);
l = l(:);

% the spatial location (find the size of the box, divide y and round down, 0 indexing)
i = floor((y-dom(1))./(Sigma./(2.^l)));

i = min(i,(2.^l)-1);

end