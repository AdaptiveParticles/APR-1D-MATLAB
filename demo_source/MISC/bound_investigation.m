%
%
%
%   Bevan Cheeseman 2017
%
%   Emperically Investigating the bound relationship between R_b and R^*
%
%
%

%
%   Run different order of approximation and compare with continous solutions demo
%
%   Works for symbollically defined functions that are differentiable up to
%   the order you wish to bound.
%
%

clear all
close all

%%

E = 0.05;
a = 1;
dom = [-.5,.5];

% define symulink symbolic function (so it can take arbitrary derivatives)
x_o = 0.02;
s = 0.01;

%currently the scale works for 0 centered functions to if you want to use
%offset functions, make sure you change that!

syms x
%f_s(x) = exp(-(x-x_o)^2/s) - exp(-(x-0.3)^2/5) + exp(-(x+5)^2/.1) ;
f_s(x) = exp(-(x-x_o)^2/s);

plot_flag = true;

compare_optimal = true;

check_derivatives = true;

%% Get APR

N_o = 4000;

y_n = linspace(dom(1),dom(2),N_o);

f = matlabFunction(f_s(x));
f_n = f(y_n);

apr = get_apr_1D_numeric(E,dom,f_n,y_n,a,2);

apr.f = f;

%compare to the apr with only the first derivative
apr_f = get_apr_1D_sym(E,dom,f_s,1,1);

%apr = apr_f;

%% Result Summary

disp(['Pulling Scheme Alg took: ',num2str(apr.calc_v_time)]);
disp(['Pulling Scheme Alg (full) took: ',num2str(apr_f.calc_v_time)]);
disp(['num_parts: ',num2str(apr.num_parts),' full_res: ',num2str(2^apr.l_max)])

disp(['num_parts a = 1: ',num2str(apr_f.num_parts),' full_res: ',num2str(2^apr_f.l_max)])

figure;
plot(apr.y_p,apr.f_p,'x-')
xlabel('y_p')
ylabel('f_p')
title('APR Sampling')

%% Test Accuracy



check_f_recon(apr,plot_flag);

%% Look at the Local Resolution Estimate Functions
figure; hold on;

y = y_n;

for i = 2:a
    plot(y,apr.L{i},'Displayname',['L:',num2str(i)])
end
plot(y,apr.L_f,'Displayname','L final');
plot(y,apr.L{1},'x-','Displayname',['L:',num2str(1)])

plot(apr_f.y_p,(apr_f.s_dom(2) - apr_f.s_dom(1))./(2.^apr_f.c_l),'o-','Displayname','R optimal-1');
plot(apr.y_p,(apr.s_dom(2) - apr.s_dom(1))./(2.^apr.c_l),'*-','Displayname','R optimal-a');

res_max = max((apr.s_dom(2) - apr.s_dom(1))./(2.^apr.c_l));

ylim([0,apr.s_dom(2)-apr.s_dom(1)])
ylim([0,res_max])
legend('show')
xlabel('y')
ylabel('Resolution')

%% True Optimal Solution Comparisons
if (compare_optimal)
    
    apr.a = a;
    
    tic;
    % reconstruction bound brute force solution with continuous resolution
    R_opt =  compute_optimal_brute_force_n(apr,y);
    time_bf = toc;
    
    tic;
    % resolution condition bound continous resolution solution
    [R_opt_b,y_b] =  compute_resolution_bound_cont_n(apr,y);
    time_cont = toc;
    
    % timing results
    disp(['Pulling Scheme Alg took: ',num2str(apr.calc_v_time)]);
    disp(['Brute Force Cont Alg took: ',num2str(time_bf)]);
    disp(['Resolution Bound Cont took: ',num2str(time_cont)]);
    
    figure;plot(y,R_opt,'Displayname','R continous brute-force')
    hold on
    plot(y_b,R_opt_b,'Displayname','R continous bound')
    
    plot(apr.y_p,(apr.s_dom(2) - apr.s_dom(1))./(2.^apr.c_l),'*-','Displayname','R optimal-a');
    
    ylim([0,max(R_opt_b)])
    legend('show')
    xlabel('Spatial Domain (y)')
    ylabel('Resolution R(y)')
    
    %difference, in terms of integral (1/R)
    int_R_opt = sum(min(diff(y))./R_opt);
    int_R_opt_b = sum(min(diff(y))./R_opt_b);
    int_R_apr = apr.num_parts;
    
    disp(['Estimated Perfect Sampling Brute Force: ',num2str(int_R_opt)]);
    disp(['Estimated Perfect Sampling Brute Force: ',num2str(int_R_opt_b)]);
    disp(['APR number of particles: ',num2str(apr.num_parts)]);
    disp(['Ratio: ',num2str(apr.num_parts/int_R_opt)]);
    disp(['Ratio rb: ',num2str(apr.num_parts/int_R_opt_b)]);
    
    figure;plot(min(diff(y))./R_opt)
    hold on;plot(min(diff(y))./R_opt_b)
    
end

%% Compare the reconstruction of the derivatives

%drawing squares.
figure;

sig = (apr.s_dom(2) - apr.s_dom(1));
for i = 1:length(apr.c_l)
    hold on;
    size = sig/2^apr.c_l(i);
    loc = (apr.c_i(i))*size + apr.s_dom(1);
    
    rectangle('Position',[loc,0,size,size],'FaceColor',[0 .5 .5],'EdgeColor','b',...
    'LineWidth',1)
    
end

[R_i,y_i] = calc_implied_res_func(apr.c_i,apr.c_l,apr.s_dom);

plot(y,apr.L_f,'Displayname','L final');



plot(y,R_opt,'Displayname','R continous brute-force')
hold on
plot(y_b,R_opt_b,'Displayname','R continous bound')

plot(y_i,R_i,'k','Displayname','R^* implied particle cells','LineWidth',2);

%[R_i,y_i]=calc_implied_res_func_doubled(apr.c_i,apr.c_l,apr.s_dom);

%plot(y_i,R_i,'k','Displayname','R^* implied particle cells','LineWidth',2);

%plot(apr.y_p,(apr.s_dom(2) - apr.s_dom(1))./(2.^apr.c_l),'*-','Displayname','R optimal-a');

axis equal
ylim([0,0.5])
xlim([dom(1),dom(2)])

xlabel('Spatial Domain y')
ylabel('Resolution Function')
legend('show')

%% Checking the actual sampling attempts

[x_p,temp] = cont_sample_particles(R_opt_b,y_b);

figure;plot(x_p,f(x_p),'x-')

[x_p2,temp2] = cont_sample_particles(R_opt,y);

hold on;plot(x_p2,f(x_p2),'x-')

plot(apr.y_p,apr.f_p,'o-')

%% Bounding it from below

R = calc_implied_res_func_y(apr.c_i,apr.c_l,apr.s_dom,y);

figure;plot(y,R_opt_b/2,'Displayname','R continous brute-force')
hold on
plot(y_b,R_opt_b,'Displayname','R continous bound')
    
plot(apr.y_p,(apr.s_dom(2) - apr.s_dom(1))./(2.^apr.c_l),'*-','Displayname','R optimal-a');

%% Gradient of the solution

grad_R = gradient(R_opt_b,min(diff(y)));
figure;
plot(grad_R);

figure;
histogram(grad_R);

%%

figure;
plot(y,1./R)
hold on
plot(y,1./R_opt_b)
plot(y,1./(.5*R_opt_b))

%%

figure;
plot(y,R./R_opt_b)


%%

N = 5;

x = cumsum(2.^(0:N));
yq = (0:N);

figure;plot(x,yq)
hold on
%plot(x,log(x)*2)

plot(x,log2(x+1) - 1)

%plot(x,log(x))

%plot(x,log(x)/1.3699)


%%

x = linspace(1,128,1000);

rat = @(x,A,B,m) log(m*x+B)./(m*log(x+A));

figure;plot(x,rat(x,0,0,.5))

%%

x = linspace(1,1000,1000);

del = @(x,A,B,m) (log(m*x+B)/m-log(x+A))./(log(x+A));

figure;plot(x,del(x,0,0,.45))
