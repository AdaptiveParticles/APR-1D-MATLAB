%
%   Bevan Cheeseman 2018
%
%   Here we show that the equivalence optimization gives the same solution
%   but for a smaller computational and memory cost, reduces costs by 2^d
%   where d is the dimension, so expect reduction of cost by 2 in this
%   example
%

clear all
close all

%%

E = 0.05; %set the relative error
dom = [-10,10];

% define symulink symbolic function (so it can take arbitrary derivatives)
x_o = 0.5;
s = 0.005;

%currently the scale works for 0 centered functions to if you want to use
%offset functions, make sure you change that!
syms x
f_s(x) = exp(-(x-x_o)^2/s) - exp(-(x-0.3)^2/5) + exp(-(x+5)^2/.1) ;

%run 5 times and take the last to get a fair timing comparison

for i = 1:10
    apr = get_apr_1D_sym(E,dom,f_s,1,1);
end

for i = 1:10
    apr_equiv = get_apr_1D_sym(E,dom,f_s,1,2);
end

disp(['Pulling Scheme Alg took: ',num2str(apr.calc_v_time)]);
disp(['Pulling Scheme Alg (Equivalence) took: ',num2str(apr_equiv.calc_v_time)]);

disp(['Total number of particles: ',num2str(apr.num_parts)]);
disp(['Total number of particles (Equivalence): ',num2str(apr_equiv.num_parts)]);


%% Memory cost difference

%first without equivalence optimization
count_C = 0;
count_L = 0;

for i = 1:length(apr.L_pc)
   for j = 1:length(apr.L_pc{i})
       
       count_C = count_C + 1;
       
       if(apr.L_pc{i}(j) > 0)
          count_L = count_L + 1; 
       end
   end
end

disp(['Size of whole tree: ',num2str(count_C)])
disp(['Size of L: ',num2str(count_L)])

%with equivalence optimition

count_Ce = 0;
count_Le = 0;

for i = 1:length(apr_equiv.L_pc)
   for j = 1:length(apr_equiv.L_pc{i})
       
       count_Ce = count_Ce + 1;
       
       if(apr_equiv.L_pc{i}(j) > 0)
          count_Le = count_Le + 1; 
       end
   end
end

disp(['Size of whole tree (Equivalence Optimization): ',num2str(count_Ce)])
disp(['Size of L (Equivalence Optimization): ',num2str(count_Le)])

%% Plot the two solutions

figure;
plot(apr.y_p,apr.f_p,'x-')
hold on
plot(apr_equiv.y_p,apr_equiv.f_p,'o-')

title('Comparison of Solutions')
xlabel('Particle Locations y_p');
ylabel('Particle Property f_p');

