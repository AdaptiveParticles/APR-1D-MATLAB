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
L1 = [0.005,0.003,0.02];
s_dom = [0,1]';
y = [0.2,0.3,0.8];

apr_1.l_max = l_max;
apr_1.l_min = 1;
apr_1.s_dom = s_dom;

syms x
f(x) =  exp(-(x-0.3)^2/5) + exp(-(x+5)^2/.1) ;

L_pc_1 = create_local_particle_set(y,L1,s_dom,l_max,type);

V_pc_1 = pulling_scheme_method_ex(L_pc_1,1,l_max);

apr_1 = sample_apr_type(V_pc_1,apr_1,f);

[R_i,y_i] = calc_implied_res_func(apr_1.c_i,apr_1.c_l,s_dom);

apr = apr_1;

cm_type = 'parula(4)';

cm = colormap(cm_type);

%drawing squares.
figure1 = figure('Color',[0.2 0.2 0.2]);

%whitebg(gcf, [1,1,1])
figure1.Position = [600,600,1400,500];

subplot(2,1,1,'Parent',figure1);
axis equal
ylim([0,0.15])
xlim([s_dom(1),s_dom(2)])

axis off

subplot(2,1,2,'Parent',figure1);
xlim([s_dom(1),s_dom(2)])
ylim([0,11]);
axis off

pause on;

sig = (apr.s_dom(2) - apr.s_dom(1));

for l = apr.l_max:-1:apr.l_min
    
    indx = find(apr.c_l == l);
    
    indx_L = find(L_pc_1{l}==1);
    
    if(isempty(indx_L))
        indx_L = -1;
    end
    
    for j = 1:length(indx)
        
        i = indx(j);
        
        if(((apr.c_type(i)) == 1))
            subplot(2,1,1,'Parent',figure1);
            title('Pulling Scheme and Seperable Property of V - Particle Cells and Implied Resolution Function R^*','Color','w')
            hold on;
            size = sig/2^apr.c_l(i);
            loc = (apr.c_i(i))*size + apr.s_dom(1);
            loc_2 = (apr.c_i(i)+0.5)*size + apr.s_dom(1);
            
            rectangle('Position',[loc,0,size,size],'FaceColor',cm(2,:),'EdgeColor',cm(1,:),...
                'LineWidth',0.5)
            ylabel('Particle Cells','Color','w','FontSize',12)
            subplot(2,1,2,'Parent',figure1);
            title('Particle Cell Level l (y-axis)','Color','w')
            hold on;
            scatter(loc_2,l,30,cm(2,:),'filled')
            
            drawnow();
          
        end
        
        %pause(0.1);
    end
    
    for j = 1:length(indx)
        
        i = indx(j);
        
        if(((apr.c_type(i)) == 2))
            
            
            subplot(2,1,1,'Parent',figure1);
            hold on;
            size = sig/2^apr.c_l(i);
            loc = (apr.c_i(i))*size + apr.s_dom(1);
            loc_2 = (apr.c_i(i)+0.5)*size + apr.s_dom(1);
            
            rectangle('Position',[loc,0,size,size],'FaceColor',cm(3,:),'EdgeColor',cm(1,:),...
                'LineWidth',0.5)
        
            subplot(2,1,2,'Parent',figure1);
     
            hold on;
            scatter(loc_2,l,30,cm(3,:),'filled')
            
            drawnow();

        end
        
        pause(0.1);
    end
    
    for j = 1:length(indx)
        
        i = indx(j);
        
        if(((apr.c_type(i)) == 3))
            
            
            subplot(2,1,1,'Parent',figure1);
            hold on;
            size = sig/2^apr.c_l(i);
            loc = (apr.c_i(i))*size + apr.s_dom(1);
            loc_2 = (apr.c_i(i)+0.5)*size + apr.s_dom(1);
            ylabel('Particle Cells','Color','w','FontSize',12)
            rectangle('Position',[loc,0,size,size],'FaceColor',[.6,.6,.6],'EdgeColor',cm(1,:),...
                'LineWidth',0.5)
            
            subplot(2,1,2,'Parent',figure1);
            hold on;
            scatter(loc_2,l,30,[.6,.6,.6],'filled')
            
            drawnow();

        end
        
        pause(0.1);
    end
    
end



