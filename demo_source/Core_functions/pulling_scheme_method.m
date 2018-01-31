function V_pc = pulling_scheme_method(L_pc,l_min,l_max)
%
%   Bevan Cheeseman 2017
%
%   Pulling Scheme
%
%   Makes use of Theorems.3 and Theorems.4, to construct the valid optimal
%   set more or less directly
%

EMPTY = 0;
SEED = 1;
BOUNDARY = 2;
FILLER = 3;
PROPOGATE = 4;
ASCENDANT = 5;
ASCENDANT_NEIGH = 6;

l_max = length(L_pc) + l_min - 1;

%loop from highest level and resolution to lowest level and resolution
for l_current = (l_max:-1:l_min)
    if(l_current~=l_max)
        set_ascendant_neighbours(l_current);
        
        set_filler(l_current)
    end
    
    fill_neighbours(l_current)
    
end

%does it in-place so copy across
V_pc = L_pc;

    function fill_neighbours(l)
        
        for i = 0:(length(L_pc{l})-1)
            type = L_pc{l}(i+1);
            if(type == SEED || type == PROPOGATE)
                %loop over neighbours
                [i_n,~] = get_neigh(i,l);
                
                %if empty place a boundary particle cell
                for n = 1:length(i_n)
                    if(L_pc{l}(i_n(n)+1) == EMPTY)
                        L_pc{l}(i_n(n)+1) = BOUNDARY;
                    end
                end
                
                %apply theorem 4.
                set_parent(i,l)
            elseif (type == ASCENDANT)
                %apply theorem 4.
                set_parent(i,l)
            end
            
        end
        
    end
    function set_ascendant_neighbours(l)
        for i = 0:(length(L_pc{l})-1)
            type = L_pc{l}(i+1);
            
            %spread the filler cells removing redundancy by doing this at
            %the one level up
            if(type == ASCENDANT)
                %loop over neighbours
                [i_n,~] = get_neigh(i,l);
                
                %if empty place a boundary particle cell
                for n = 1:length(i_n)
                    if(L_pc{l}(i_n(n)+1) == EMPTY)
                        L_pc{l}(i_n(n)+1) = ASCENDANT_NEIGH;
                    elseif(L_pc{l}(i_n(n)+1) == SEED)
                        L_pc{l}(i_n(n)+1) =PROPOGATE;
                    end
                end
            end        
        end
        
    end
    function set_filler(l)
        %
        %
        %
        %
        
        for i = 0:(length(L_pc{l})-1)
            type = L_pc{l}(i+1);
            
            %spread the filler cells removing redundancy by doing this at
            %the one level up
            if(type == ASCENDANT_NEIGH || type == PROPOGATE)
                %loop over neighbours
                [~,i_c] = get_children(i,l);
                
                %if empty children place a filler cell
                for c = 1:length(i_c)
                    if(L_pc{l+1}(i_c(c)+1) == EMPTY)
                        L_pc{l+1}(i_c(c)+1) = FILLER;
                    end
                end
            end
        end
        
    end
    function set_parent(i,l)
        %
        %   Spreads the filler cells and removes seed cells above by Theorem.4
        %
        %
        if(l>l_min)
            [p_i,p_l] = get_parent(i,l);
            type = L_pc{p_l}(p_i+1);
            %if(type ~= SEED)
                L_pc{p_l}(p_i+1) = ASCENDANT;
            %end
            
        end
        
        
    end



end



