function visualize_pc(apr)

a = apr.a;

higher_res = 1;

%search sampling points, we have already found the lowest resolution required.
y = (apr.s_dom(1)):((apr.s_dom(2) - apr.s_dom(1))/2^(apr.l_max+higher_res)):apr.s_dom(2);


% see multiple resolution conditions, for the general bound
for i = 1:a
    apr.L{i} = ((apr.E*apr.scale(i))./abs(apr.df{i}(y))).^(1/i);
end

apr.L_f = apr.L{1};

% now take the minimum across them
for i = 2:a
   apr.L_f = min(apr.L_f,apr.L{i}); 
end


type = 1;



L_pc = create_local_particle_set(y,apr.L_f,apr.s_dom,apr.l_max,type);

%% Calculate V the optimal valid particle cell set

V_pc = pulling_scheme_method(L_pc,apr.l_min,apr.l_max);

l_min = 1;
l_max = apr.l_max;

sig = (apr.s_dom(2) - apr.s_dom(1));

figure;
hold on

for l = l_min:l_max
    for i = 0:(length(V_pc{l})-1)
        
        full = L_pc{l}(i+1);
        
        if(full)
            size = sig/2^l;
            loc = (i)*size + apr.s_dom(1);
            height = size/2;
    
            rectangle('Position',[loc,height,size,size/2],'FaceColor',[0 .5 .5],'EdgeColor','b',...
            'LineWidth',1)
            
        end
        
    end
    
end

axis equal;
ylim([0,sig/2])
xlim([apr.s_dom(1),apr.s_dom(2)])

plot(y,apr.L_f,'x-')

figure;
hold on

for l = l_min:l_max
    for i = 0:(length(V_pc{l})-1)
        
        full = V_pc{l}(i+1);
        
        if((full < 4) && (full > 0))
            size = sig/2^l;
            loc = (i)*size + apr.s_dom(1);
            height = size/2;
    
            rectangle('Position',[loc,height,size,size/2],'FaceColor','b','EdgeColor','g',...
            'LineWidth',1)
            
        end
        
    end
    
end

axis equal
ylim([0,sig/2])
xlim([apr.s_dom(1),apr.s_dom(2)])
plot(y,apr.L_f,'x-')

figure;
hold on

for l = l_min:l_max
    for i = 0:(length(V_pc{l})-1)
        
        full = V_pc{l}(i+1);
        
        if((full < 4) && (full > 0))
            size = sig/2^l;
            loc = (i)*size + apr.s_dom(1);
            height = size/2;
    
            rectangle('Position',[loc,0,size,size],'FaceColor','b','EdgeColor','g',...
            'LineWidth',1)
            
        end
        
    end
    
end

axis equal
ylim([0,sig/2])
xlim([apr.s_dom(1),apr.s_dom(2)])
plot(y,apr.L_f)


end