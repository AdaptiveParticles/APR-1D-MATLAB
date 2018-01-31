function apr = get_apr_1D_sym(E,dom,f_s,a,varargin)
%
%
%   Bevan Cheeseman 2017
%
%   Takes in a simulink function and then creates the APR, must have variable x defined as the main variable.
%   maintaining
%   |df^(k-1) - \hat{df^(k-1)}| <= E for up to (k) <= a orders. a = 1 just
%   function, a = 2, dfdx ect
%
%

if(~isempty(varargin))
   %optimized up a resolution level
   type = varargin{1};
   
   if(length(varargin) > 1)
      apr.sample_depth = varargin{2}; 
   end
   
else
   type = 1; 
end

if length(a) ==2
    a_min = a(1);
    a_max = a(2);
else
    a_min = 1;
    a_max = a;
end


apr.s_dom = dom;
apr.E = E;
apr.l_min = 1;

syms x;
apr.f = matlabFunction(f_s);

%% set up the different derivative functions to be evaluated
for i = 1:a_max
    apr.df{i} = matlabFunction(diff(f_s,i)); 
end

%% compute the constant scale functions
apr.scale(1) = calc_const_scale(apr.f,apr.s_dom);

for i = 1:(a_max-1)
    apr.scale(i+1) = calc_const_scale(apr.df{i},apr.s_dom);
end

%% compute the maximum level
apr.l_max = 0;



for i = a_min:a_max
    apr.l_max = max(apr.l_max,calc_l_max(apr.df{i},apr.s_dom,apr.E,apr.scale(i),i));
end


%% now construct L, local resolution estimate

y = (apr.s_dom(1)):((apr.s_dom(2) - apr.s_dom(1))/2^apr.l_max):apr.s_dom(2);

% see multiple resolution conditions, for the general bound
for i = a_min:a_max
    apr.L{i} = ((apr.E*apr.scale(i))./abs(apr.df{i}(y))).^(1/i);
end

apr.L_f = 999999*ones(size(y));

% now take the minimum across them
for i = a_min:a_max
   apr.L_f = min(apr.L_f,apr.L{i}); 
end

%% construct the local particle cell set

L_pc = create_local_particle_set(y,apr.L_f,apr.s_dom,apr.l_max,type);

%% Calculate V the optimal valid particle cell set
tic;

V_pc = pulling_scheme_method(L_pc,apr.l_min,apr.l_max);

apr.calc_v_time = toc;

%% Now construct the APR samplign the function


 apr = sample_apr(V_pc,apr,apr.f,type);

end