%
%   Bevan Cheeseman 2018   
%
%   This demo illustrates the APR satisfying the Reconstuction Condition,
%   for varying E, and reconstruction methods.
%   
%   Note: Requires pre-computation of a scale and maximum gradient over the
%   domain to calculate the required initial sampling for each E.
%

clear all
close all

%%

a = 1;
dom = [-10,10];

% define symulink symbolic function (so it can take arbitrary derivatives)
x_o = 0.5;
s = 0.05;

%currently the scale works for 0 centered functions to if you want to use
%offset functions, make sur eyou cahnge that!

syms x
f_s(x) = exp(-(x-x_o)^2/s) - exp(-(x-0.3)^2/5) + exp(-(x+5)^2/.1) ;

f = matlabFunction(f_s);
%% need to find the maximum gradient on the domain so we can set a maximum
% sampling for a given E.

N_search = 400;

y_search = linspace(dom(1),dom(2),N_search);
f_search = f(y_search);
%domain size
Sigma = (dom(2) - dom(1));
grad_search = get_grad_k(y_search,f_search,1,1); %just computes a numerical estimate of the gradient
max_grad = max(abs(grad_search));

%compute the constant local scale
scale_estimate = calc_const_scale_num(f_search);

%% Get APR for different reconstruction errors E

E_min = 0.01;
E_max = 1;
N = 20;

E_ = linspace(E_min,E_max,N);

linear_inf = zeros(size(E_));
pc_inf = zeros(size(E_));
wc_inf = zeros(size(E_));


for i = 1:N
    % this then gives us our minimum value of L
    L_min = ((E_(i)*scale_estimate)/max_grad);

    %then find the l_max that rounds to that resolution
    l_max = ceil(log2(Sigma/L_min));
    
    %input samples
    y_numerical = linspace(dom(1),dom(2),2^(l_max));
    f_numerical = f(y_numerical);
       
    apr = get_apr_1D_numeric(E_(i),dom,f_numerical,y_numerical,1,1);
    
    %the reconstruction compatison requires the function handle
    apr.f = f;
    
    [linear_inf(i),pc_inf(i),wc_inf(i)] = check_f_recon(apr,false);
    
end

%% 

figure;
plot(E_,E_,'--','Displayname','Reconstruction Condition');hold on
plot(E_,wc_inf,'Displayname','Worst Case Recon')
plot(E_,linear_inf,'Displayname','Linear Recon')
plot(E_,pc_inf,'Displayname','Nearest Neighbour PC Recon')

xlabel('Relative Error (E)')
ylabel('Observed |E|_{\infty}')
title('Reconstruction Condition')
legend('show')

%% 