%
%   Bevan Cheeseman 2018
%
%   Here we compare the the impact of having full knowledge of the input
%   function, or only knowing the function point-wise and having to compute
%   its gradient and sampling.
%
%   The plots show the Local Resolution Estimate L(y), compared to the
%   computed Implied Resolution Function R^* from V.
%

clear all
close all

%%

E = 0.1; %set the relative error
a = 1; % the number of derivative conditions you wish to set.
dom = [-10,10];

% define symulink symbolic function (so it can take arbitrary derivatives)
x_o = 0.5;
s = 0.005;

%currently the scale works for 0 centered functions to if you want to use
%offset functions, make sure you change that!
syms x
f_s(x) = exp(-(x-x_o)^2/s) - exp(-(x-0.3)^2/5) + exp(-(x+5)^2/.1) ;

% we are comparing with the analytical solution so will use its l_max to
% set y sampling search
apr_f = get_apr_1D_sym(E,dom,f_s,a,1);

A = matlabFunction(f_s);
f = matlabFunction(f_s);

y_numerical = (apr_f.s_dom(1)):((apr_f.s_dom(2) - apr_f.s_dom(1))/2^apr_f.l_max):apr_f.s_dom(2);
f_numerical = A(y_numerical);

%% Get APR

apr = get_apr_1D_numeric(E,dom,f_numerical,y_numerical,a,1);

%% Result Summary

disp(['Pulling Scheme Alg took: ',num2str(apr.calc_v_time)]);
disp(['Pulling Scheme Alg (full) took: ',num2str(apr_f.calc_v_time)]);
disp(['num_parts numeric: ',num2str(apr.num_parts),' full_res: ',num2str(2^apr.l_max)])

disp(['num_parts analytical: ',num2str(apr_f.num_parts),' full_res: ',num2str(2^apr_f.l_max)])

figure;
plot(apr.y_p,apr.f_p,'x-')
xlabel('y_p')
ylabel('f_p')
title('APR Sampling')

%% Test Accuracy

plot_flag = true;

%% Look at the Local Resolution Estimate Functions
figure; hold on;

y = (apr.s_dom(1)):((apr.s_dom(2) - apr.s_dom(1))/2^apr.l_max):apr.s_dom(2);

plot(y,apr.L{1},'x-','Displayname',['L(y)'])

plot(apr_f.y_p,(apr_f.s_dom(2) - apr_f.s_dom(1))./(2.^apr_f.c_l),'o-','Displayname','R optimal-analytical');
plot(apr.y_p,(apr.s_dom(2) - apr.s_dom(1))./(2.^apr.c_l),'*-','Displayname','R optimal-numerical');

ylim([0,apr.s_dom(2)-apr.s_dom(1)])
legend('show')
xlabel('y')
ylabel('Resolution')

%% Reconstruct and check
f_r = interp1(apr.y_p,apr.f_p,y);
f_gt = f(y);
recon_err = max(abs(f_gt-f_r));

disp(['Desired: ',num2str(apr.E),' Linear Interpolation E_obs: ',num2str(recon_err)])

%% Nearest Neighbour Reconstruction

f_r = interp1(apr.y_p,apr.f_p,y,'nearest');
f_gt = f(y);
recon_err = max(abs(f_gt-f_r));

disp(['Desired: ',num2str(apr.E),' Nearest Neighbour E_obs: ',num2str(recon_err)])

%% Worse Case Reconstruction

[max_rc,min_rc,y] = worst_case(apr);
recon_max = max(abs(f(y)-max_rc));
recon_min = max(abs(f(y)-min_rc));
recon_worst = max(recon_max,recon_min);

disp(['Desired: ',num2str(apr.E),' Worst Case E_obs: ',num2str(recon_worst)])

