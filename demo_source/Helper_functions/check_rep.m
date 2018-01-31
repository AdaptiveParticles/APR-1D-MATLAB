function check_rep(V_pc,L_pc,l_min,l_max)

counter = 0;

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
            
            counter = counter + 1;
                       
        end
        
    end
    
end

disp(['Total number of Particle Cells: ',num2str(counter)]);

end