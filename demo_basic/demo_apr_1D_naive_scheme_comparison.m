%
%   Bevan Cheeseman 2018
%   APR Demo in 1D
%
%   This demo, shows the comparison of computing the solution using the
%   Pulling Scheme, and instead a naive brute force approach using Particle
%   Cells. Both have a worst-case scaling of O(N), however the Pulling
%   Scheme takes advantage of the unique features of V, namely the
%   Seperability, Predictability, Redundancy, and Equivalence.
%

close all
clear all

%% Set things up
s1 = 0.05;
s2 = 0.001;
a = 0.5;
b = -0.3;

%input function
f = @(x) exp((-(x-a).^2)/s1) - exp((-(x-b).^2)/s2);
%known gradient
grad = @(x) -(2/s1)*(x-a).*exp((-(x-a).^2)/s1) + (2/s2)*x.*exp((-(x-b).^2)/s2);

%set the size of the domain
apr.s_dom = [-2,2];

%set the reletaive error
apr.E = 0.1;

%use a constant local intensity scale \sigma(y)
apr.scale = calc_const_scale(f,apr.s_dom);

%level min and max of Particle Cells
apr.l_min = 1;
apr.l_max = calc_l_max(grad,apr.s_dom,apr.E,apr.scale,1);

y = (apr.s_dom(1)):((apr.s_dom(2) - apr.s_dom(1))/2^apr.l_max):apr.s_dom(2);

%% Calculate the Local Resolution Estimate L(y) the local particle cell set

L = apr.E*apr.scale./(abs(grad(y)));

%% Calculate the Local Particle Cell Set \mathcal{L}

L_pc = create_local_particle_set(y,L,apr.s_dom,apr.l_max);

%% Calculate \mathcal{V} the optimal valid particle cell set using the naive approach
tic;

%this involves iterating over the whole tree from lowest level to highest
%level (Note: this is **NOT** the Pulling Scheme and is given for comparison)
% however, it also uses Particle Cells (For continuous solution comparison see: )
V_pc_n = naive_method(L_pc,apr.l_min,apr.l_max);

calc_naive_time = toc;

%% Calculate \mathcal{V} the optimal valid particle cell set using the Pulling Scheme
tic;

V_pc = pulling_scheme_method(L_pc,apr.l_min,apr.l_max);

calc_pulling_scheme_time = toc;

%% Check that the are both the valid solution and the same

%pulling scheme result
check_rep(V_pc,L_pc,apr.l_min,apr.l_max);

%naive scheme result
check_rep(V_pc_n,L_pc,apr.l_min,apr.l_max);

%% Now construct the APR samplign the function

apr = sample_apr(V_pc,apr,f);

disp(['Pulling Scheme Alg took: ',num2str(calc_pulling_scheme_time)]);
disp(['Naiive Alg took: ',num2str(calc_naive_time)]);



