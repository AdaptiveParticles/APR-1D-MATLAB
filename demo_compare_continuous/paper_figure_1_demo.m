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

E = 0.1;
a = 1;
dom = [-.5,.5];

% define symulink symbolic function (so it can take arbitrary derivatives)
x_o = 0.02;
s = 0.01;

%currently the scale works for 0 centered functions to if you want to use
%offset functions, make sure you change that!

syms x
%f_s(x) = exp(-(x-x_o)^2/s) - exp(-(x-0.3)^2/.1) + exp(-(x+5)^2/.1) ;
f_s(x) = exp(-(x-x_o)^2/s);

plot_flag = true;

compare_optimal = true;

%% Get APR

N_o = 4000;

y_n = linspace(dom(1),dom(2),N_o);

f = matlabFunction(f_s(x));
f_n = f(y_n);

%duplicated just for timing
apr = get_apr_1D_numeric(E,dom,f_n,y_n,a,1);

apr = get_apr_1D_numeric(E,dom,f_n,y_n,a,1);

apr.f = f;

%compare to the apr with only the first derivative
apr_f = get_apr_1D_sym(E,dom,f_s,1,1);

%apr = apr_f;

%% Result Summary


disp(['num_parts: ',num2str(apr.num_parts),' full_res: ',num2str(2^apr.l_max)])
disp(['num_parts sym: ',num2str(apr_f.num_parts),' full_res: ',num2str(2^apr_f.l_max)])

figure;
plot(apr.y_p,apr.f_p,'x-')
xlabel('y_p')
ylabel('f_p')
title('APR Sampling')

%% Test Accuracy
plot_flag_recon = false; %if you want output of the reconstruction chance to true
check_f_recon(apr,plot_flag_recon);

%% Look at the Local Resolution Estimate Functions
figure; hold on;

y = y_n;

plot(y,apr.L_f,'Displayname','L(y)');

plot(apr_f.y_p,(apr_f.s_dom(2) - apr_f.s_dom(1))./(2.^apr_f.c_l),'o-','Displayname','R^*');

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
    
    disp(['Estimated Perfect Sampling R_c: ',num2str(int_R_opt)]);
    disp(['Estimated Perfect Sampling R_b: ',num2str(int_R_opt_b)]);
    disp(['APR number of particles: ',num2str(apr.num_parts)]);
    disp(['Ratio APR - Continuous Parts: ',num2str(apr.num_parts/int_R_opt)]);
    disp(['Ratio APR - Conditnous Parts R_b: ',num2str(apr.num_parts/int_R_opt_b)]);
    
%     figure;plot(min(diff(y))./R_opt)
%     figure;plot(min(diff(y))./R_opt_b)
    
end

%% Compare the reconstruction of the derivatives

disp(['Pulling Scheme Alg took: ',num2str(apr.calc_v_time)]);
disp(['Computing R_c took: ',num2str(time_bf)]);
disp(['Computing R_b took: ',num2str(time_cont)]);

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

plot(y,apr.L_f,'Displayname','L(y)');

plot(y,R_opt,'Displayname','R_c')
hold on
plot(y_b,R_opt_b,'Displayname','R_b')
plot(y_i,R_i,'k','Displayname','R^* implied particle cells','LineWidth',2);

axis equal
ylim([0,0.5])
xlim([dom(1),dom(2)])

xlabel('Spatial Domain y')
ylabel('Resolution Function')
legend('show')
title('Comparison between optimal implied and continuous Resolution Functions')
