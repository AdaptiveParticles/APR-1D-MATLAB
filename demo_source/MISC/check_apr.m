function check_apr(apr)

type = 1;

y = (apr.s_dom(1)):((apr.s_dom(2) - apr.s_dom(1))/2^apr.l_max):apr.s_dom(2);

L_pc = create_local_particle_set(y,apr.L_f,apr.s_dom,apr.l_max,type);

%% Calculate V the optimal valid particle cell set

V_pc = pulling_scheme_method(L_pc,apr.l_min,apr.l_max);

l_min = 1;
l_max = apr.l_max;

for l = l_min:l_max
    for i = 0:(length(V_pc{l})-1)
        
        type = V_pc{l}(i+1);
        
        if((type > 0) && (type <= 3))
            
            empty = true;
            
            %get your neighbours
            [i_n,l_n] = get_neigh(i,l);
            
            %check if your neighbour descendents are in L
            for n = 1:length(i_n)
                temp = check_descendants(i_n(n),l_n(n),L_pc,l_max);
                empty = min(temp,empty);
                
            end
            
            if(empty == false)
               disp('BROKEN'); 
            end
                       
        end
        
    end
    
end




end