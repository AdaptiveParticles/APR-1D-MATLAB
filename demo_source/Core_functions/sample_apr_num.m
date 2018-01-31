function apr = sample_apr_num(V_pc,apr,f,varargin)
%first get valid particle cells
if(~isempty(varargin))
   type = varargin{1}; 
    
else
    type = 1;
end


apr.c_i = [];
apr.c_l = [];


if type == 1
    %simple routine
    
    % now get the set
    for l = (apr.l_max):-1:apr.l_min
        for i = 0:(length(V_pc{l})-1)
            % first check is it valid
            if(V_pc{l}(i+1)>0 && V_pc{l}(i+1)<4)
                apr.c_i(end+1) = i;
                apr.c_l(end+1) = l;
            end
        end
    end

elseif (type ==2)
    %routine with augmented type
     % now get the set
    for l = (apr.l_max-1):-1:apr.l_min
        for i = 0:(length(V_pc{l})-1)
            
            % seed and boundary resolve at highest reoslution
            if(V_pc{l}(i+1)>0 && V_pc{l}(i+1)<=2)
                [l_c,i_c] = get_children(i,l);
                
                %place particle at higher resolution
                for c = 1:length(i_c)
                    apr.c_i(end+1) = i_c(c);
                    apr.c_l(end+1) = l+1;
                end
                     
            end
            
            %filler type resolve at this level
            if(V_pc{l}(i+1)==3)
                apr.c_i(end+1) = i;
                apr.c_l(end+1) = l;
            end
            
        end
    end
    
elseif (type ==3)
    %routine with augmented type
     % now get the set
    for l = (apr.l_max-1):-1:apr.l_min
        for i = 0:(length(V_pc{l})-1)
            
            % seed and boundary resolve at highest reoslution
            if(V_pc{l}(i+1)==1)
                [l_c,i_c] = get_children(i,l);
                
                %place particle at higher resolution
                for c = 1:length(i_c)
                    apr.c_i(end+1) = i_c(c);
                    apr.c_l(end+1) = l+1;
                end
                     
            end
            
            %filler type resolve at this level
            if(V_pc{l}(i+1)==3 || V_pc{l}(i+1)==2)
                apr.c_i(end+1) = i;
                apr.c_l(end+1) = l;
            end
            
        end
    end
    
end

%% Now construct the APR and sample the particles

apr.num_parts = length(apr.c_i);
apr.y_p = zeros(apr.num_parts,1);
apr.f_p = zeros(apr.num_parts,1);

for p = 1:(apr.num_parts)
    apr.y_p(p) = apr.s_dom(1) + (apr.c_i(p)+0.5)*(apr.s_dom(2)-apr.s_dom(1))/2^apr.c_l(p);
    apr.f_p(p) = f(apr.y_p(p));
end

% now lets sort them in terms of y for ease of plotting and calculations
[~,s_index] = sort(apr.y_p);
apr.y_p = apr.y_p(s_index);
apr.f_p = apr.f_p(s_index);
apr.c_l = apr.c_l(s_index);
apr.c_i = apr.c_i(s_index);

end