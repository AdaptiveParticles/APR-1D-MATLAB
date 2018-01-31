function apr = get_apr_1D_numeric(E,dom,f_s,y_s,a,varargin)
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
   %(Equivalence Optimization) optimized up a resolution level
   type = varargin{1};
   
   if(length(varargin) > 1)
      apr.sample_depth = varargin{2}; 
   end
   
else
   type = 1; 
end

apr.s_dom = dom;
apr.E = E;
apr.l_min = 1;

apr.f = f_s;
apr.y = y_s;

%% set up the different derivative functions to be evaluated
% prev = apr.f;
% 
% for i = 1:a
%     apr.df{i} = calc_grad(apr.y,prev); 
% end

grad_order = 1;

for i = 1:a
    apr.df{i} = get_grad_k(apr.y,apr.f,i,grad_order)'; 
end

%% compute the constant scale functions
apr.scale(1) = calc_const_scale_num(apr.f);

for i = 1:(a-1)
    apr.scale(i+1) = calc_const_scale_num(apr.df{i});
end

%% compute the maximum level
apr.l_max = 0;

for i = 1:a
    apr.l_max = max(apr.l_max,calc_l_max_num(apr.df{i},apr.s_dom,apr.E,apr.scale(i),i));
end


%% now construct L, local resolution estimate

% see multiple resolution conditions, for the general bound
for i = 1:a
    apr.L{i} = ((apr.E*apr.scale(i))./abs(apr.df{i})).^(1/i);
end

apr.L_f = apr.L{1};

% now take the minimum across them
for i = 2:a
   apr.L_f = min(apr.L_f,apr.L{i}); 
end

%% construct the local particle cell set

L_pc = create_local_particle_set(apr.y,apr.L_f,apr.s_dom,apr.l_max,type);

%% Calculate V the optimal valid particle cell set
tic;

V_pc = pulling_scheme_method(L_pc,apr.l_min,apr.l_max);

apr.calc_v_time = toc;

%% Now construct the APR samplign the function

apr = sample_apr(V_pc,apr,apr.f,type,'numeric');

end