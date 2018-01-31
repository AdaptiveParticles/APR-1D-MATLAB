%
%   Bevan Cheeseman 2018
%   APR Demo in 1D
%
%   Forms the APR for a known 1D function, and compares and plots the
%   Reconstructions using different methods. Any method that takes a
%   non-negative weighted average of particles should satisfy the
%   Reconstruction Condition, (i.e. the infinity norm of the reconstructions 
%   should be less then E, with the worst-case reconstruction representation
%   representing the worst possible.)
%

close all
clear all

%% Set things up
s1 = 0.05;
s2 = 0.001;
a = 0.5;
b = -0.3;

f = @(x) exp((-(x-a).^2)/s1) - exp((-(x-b).^2)/s2);

grad = @(x) -(2/s1)*(x-a).*exp((-(x-a).^2)/s1) + (2/s2)*x.*exp((-(x-b).^2)/s2);

%domain size
apr_ps.s_dom = [-2,2];

%set the Relative Error.
apr_ps.E = 0.05;

%use constant Local Intensity Scale
apr_ps.scale = calc_const_scale(f,apr_ps.s_dom);

%maximum and minimum level
apr_ps.l_min = 1;
apr_ps.l_max = calc_l_max(grad,apr_ps.s_dom,apr_ps.E,apr_ps.scale,1);

%input sampling
y = (apr_ps.s_dom(1)):((apr_ps.s_dom(2) - apr_ps.s_dom(1))/2^apr_ps.l_max):apr_ps.s_dom(2);

%% Calculate L the local particle cell set

L = apr_ps.E*apr_ps.scale./(abs(grad(y)));

L_pc = create_local_particle_set(y,L,apr_ps.s_dom,apr_ps.l_max);


%% Calculate V the optimal valid particle cell set
tic;

V_pc = pulling_scheme_method(L_pc,apr_ps.l_min,apr_ps.l_max);

calc_pulling_scheme_time = toc;

%% Now construct the APR samplign the function

apr_ps = sample_apr(V_pc,apr_ps,f);

figure;
plot(apr_ps.y_p,apr_ps.f_p,'x-')

title('The Adaptive Particle Representation')
xlabel('Particle Locations y_p');
ylabel('Particle Property f_p');

disp(['Pulling Scheme Alg took: ',num2str(calc_pulling_scheme_time)]);

%% Plot the particle cells
figure;plot(apr_ps.y_p,apr_ps.c_l,'x')
title('V Particle Cells')
xlabel('Particle Locations y_p');
ylabel('Particle Cell level');

%% Check Rep Satisfies the Conditions

check_rep(V_pc,L_pc,apr_ps.l_min,apr_ps.l_max)

%% Reconstruct and check

f_r = interp1(apr_ps.y_p,apr_ps.f_p,y);
figure;plot(y,f_r)
f_gt = f(y);
hold on 
plot(y,f_gt)

title('Linear Reconstruction')
ylabel('f(y)');
xlabel('y');

recon_err = max(abs(f_gt-f_r));

disp(['Desired: ',num2str(apr_ps.E),' Linear Interpolation E_obs: ',num2str(recon_err)])

%% Nearest Neighbour Reconstruction

f_r = interp1(apr_ps.y_p,apr_ps.f_p,y,'nearest');
figure;plot(y,f_r)
f_gt = f(y);
hold on 
plot(y,f_gt)

title('Nearest Neighbour Reconstruction')
ylabel('f(y)');
xlabel('y');

recon_err = max(abs(f_gt-f_r));

disp(['Desired: ',num2str(apr_ps.E),' Nearest Neighbour E_obs: ',num2str(recon_err)])

%% Worse Case Reconstruction

[max_rc,min_rc,y] = worst_case(apr_ps);

figure;plot(y,max_rc);hold on
plot(y,min_rc);

recon_max = max(abs(f(y)-max_rc));
recon_min = max(abs(f(y)-min_rc));

title('Worst-case Reconstruction')
ylabel('f(y)');
xlabel('y');

recon_worst = max(recon_max,recon_min);

disp(['Desired: ',num2str(apr_ps.E),' Worst Case E_obs: ',num2str(recon_worst)])
