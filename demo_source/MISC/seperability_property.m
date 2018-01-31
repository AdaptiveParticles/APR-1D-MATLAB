%
%
%   Bevan Cheeseman 2017
%
%   Shows the seperability property, creating the optimal particle cell set
%   V from each individually, or the two joint.
%
%
 clear all
%% Create First 1

type = 1;
l_max = 14;
L1 = 0.0000001;
s_dom = [0,1]';
y = 0.3;

apr_1.l_max = l_max;
apr_1.l_min = 1;
apr_1.s_dom = s_dom;

syms x
f(x) =  exp(-(x-0.3)^2/5) + exp(-(x+5)^2/.1) ;

L_pc_1 = create_local_particle_set(y,L1,s_dom,l_max,type);

V_pc_1 = pulling_scheme_method(L_pc_1,1,l_max);

apr_1 = sample_apr(V_pc_1,apr_1,f);

[R_i,y_i] = calc_implied_res_func(apr_1.c_i,apr_1.c_l,s_dom);

figure;plot(y_i,R_i,'x-')
axis equal

%% Create Second 1

y = 0.6;
L1 = 0.001;

L_pc_2 = create_local_particle_set(y,L1,s_dom,l_max,type);

V_pc_2 = pulling_scheme_method(L_pc_2,1,l_max);

apr_2 = sample_apr(V_pc_2,apr_1,f);

[R_i,y_i] = calc_implied_res_func(apr_2.c_i,apr_2.c_l,s_dom);

hold on 
plot(y_i,R_i,'x-')

%% Together

figure
y = [0.3,0.6];
L1 = [0.001,0.001];

L_pc_2 = create_local_particle_set(y,L1,s_dom,l_max,type);

V_pc_2 = pulling_scheme_method(L_pc_2,1,l_max);

apr_2 = sample_apr(V_pc_2,apr_1,f);

[R_i,y_i] = calc_implied_res_func(apr_2.c_i,apr_2.c_l,s_dom);

hold on 
plot(y_i,R_i,'x-')
axis equal