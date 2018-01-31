%
%
%   Bevan Cheeseman 2017
%
%   Shows the seperability property, creating the optimal particle cell set
%   V from each individually, or the two joint.
%
%   (Thesis New one)
 clear all
%% Create First 1

type = 1;
l_max = 14;
L1 = 0.0000001;
s_dom = [0,1]';
y1 = 0.3;

apr_1.l_max = l_max;
apr_1.l_min = 1;
apr_1.s_dom = s_dom;

syms x
f(x) =  exp(-(x-0.3)^2/5) + exp(-(x+5)^2/.1) ;

L_pc_1 = create_local_particle_set(y1,L1,s_dom,l_max,type);

V_pc_1 = pulling_scheme_method(L_pc_1,1,l_max);

apr_1 = sample_apr(V_pc_1,apr_1,f);

[R_i,y_i] = calc_implied_res_func(apr_1.c_i,apr_1.c_l,s_dom);

figure;plot(y_i,R_i,'-')
axis equal

%% Create Second 1

y2 = 0.6;
L2 = 0.01;

L_pc_2 = create_local_particle_set(y2,L2,s_dom,l_max,type);

V_pc_2 = pulling_scheme_method(L_pc_2,1,l_max);

apr_2 = sample_apr(V_pc_2,apr_1,f);

[R_i,y_i] = calc_implied_res_func(apr_2.c_i,apr_2.c_l,s_dom);

hold on 
plot(y_i,R_i,'-')

%% Together


y = [y1,y2];
L1 = [L1,L2];

L_pc_2 = create_local_particle_set(y,L1,s_dom,l_max,type);

V_pc_2 = pulling_scheme_method(L_pc_2,1,l_max);

apr_2 = sample_apr(V_pc_2,apr_1,f);

[R_i,y_i] = calc_implied_res_func(apr_2.c_i,apr_2.c_l,s_dom);

hold on 
plot(y_i,R_i,'-')
axis equal

print('pulling_scheme','-depsc','-painters','-loose','-cmyk');

