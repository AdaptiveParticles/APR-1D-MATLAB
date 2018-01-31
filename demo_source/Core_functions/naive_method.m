function V_pc = naive(L_pc,l_min,l_max)
%
%   Bevan Cheeseman 2017
%
%   The naiive implimentation of finding the valid particle cell set
%

V_pc = cell(l_max,1);

for l = 1:l_max
   V_pc{l} = zeros(2^l,1); 
end

%
%   Checking for valid solutions, this is direct application of Theorem.~1
%   and guranteeing satisfaction of the reoslution bound.
%
for l = l_min:l_max
    for i = 0:(length(V_pc{l})-1)
        
       empty = true;
       
       %get your neighbours
       [i_n,l_n] = get_neigh(i,l);
       
       %check if your neighbour descendents are in L
       for n = 1:length(i_n)
           temp = check_descendants(i_n(n),l_n(n),L_pc,l_max);
           empty = min(temp,empty);
           
           if(~empty)
               break;
           end
       end
       
       if(empty)
           V_pc{l}(i+1) = 1;
       end
       
    end
    
end

% then find the optimal, now going from the bottom to the top, essentially
% checking the optimality theorem, if there is a parent that is also valid
% that the current rep cannot be optimal.
for l = (l_max-1):-1:l_min
    for i = 0:(length(V_pc{l})-1)
        
        % first check is it valid
        if(V_pc{l}(i+1))
            
            temp = check_descendants(i_n(n),l_n(n),L_pc,l_max);
            empty = min(temp,empty);
            
            [l_c,i_c] = get_children(i,l);
            
            for c = 1:length(l_c)
                %set the children to 0
                V_pc{l_c(c)}(i_c(c)+1) = 0;
            end
            
        end
        
    end
    
end





end