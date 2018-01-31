function apr = sample_apr(V_pc,apr,f,varargin)
%% Settings
if(~isempty(varargin))
   type = varargin{1}; 
   
   if(length(varargin) > 1)
      
      for i = 1:length(varargin)
            if(strcmp(varargin{i},'numeric'))
                sampled_type = varargin{i};
            else
                sampled_type = 'exact';
            
            end
      end
      
   else
       sampled_type = 'exact';
   end
   
else
    type = 1;
    sampled_type = 'exact';
end


apr.c_i = [];
apr.c_l = [];

%% Now Sample

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

if(isfield(apr,'sample_depth'))
    %need to resample at higher resolution set by the factor sample_depth,
    %this is requried if you want to compute higher derivatives that
    %require more particles in the support of the radius
    
    
    
    
   if(apr.sample_depth > 1)
       add_d = apr.sample_depth-1;
       
       new_c_l = zeros(length(apr.c_i)*2^add_d,1);
       new_c_i = zeros(length(apr.c_i)*2^add_d,1);
       
       add_num = 2^add_d;
       
       counter = 1;
       
       
       for i = 1:length(apr.c_i)
           new_l = apr.c_l(i) + add_d;
           new_i0 = apr.c_i(i)*add_num;
           
           new_c_l(counter:(counter+add_num-1)) = new_l;
           
           new_c_i(counter:(counter+add_num-1)) = new_i0:(new_i0 + add_num -1);
           
           counter = counter + add_num;
           
       end
       
       apr.c_i = new_c_i;
       apr.c_l = new_c_l;
       
   end
    
end


%% Now construct the APR and sample the particles

apr.num_parts = length(apr.c_i);
apr.y_p = zeros(apr.num_parts,1);
apr.f_p = zeros(apr.num_parts,1);

if(strcmp(sampled_type,'exact'))
    
    for p = 1:(apr.num_parts)
        apr.y_p(p) = apr.s_dom(1) + (apr.c_i(p)+0.5)*(apr.s_dom(2)-apr.s_dom(1))/2^apr.c_l(p);
        apr.f_p(p) = f(apr.y_p(p));
    end
    
else
    
    for p = 1:(apr.num_parts)
        apr.y_p(p) = apr.s_dom(1) + (apr.c_i(p)+0.5)*(apr.s_dom(2)-apr.s_dom(1))/2^apr.c_l(p);
    end
    
    %we need to interp to the new particles, here we are assuming non-noisy
    %input
    apr.f_p = interp1(apr.y,apr.f,apr.y_p,'linear','extrap');
    %apr.f_p = interp1(apr.y,apr.f,apr.y_p,'cubic','extrap');
end

% now lets sort them in terms of y for ease of plotting and calculations
[~,s_index] = sort(apr.y_p);
apr.y_p = apr.y_p(s_index);
apr.f_p = apr.f_p(s_index);
apr.c_l = apr.c_l(s_index);
apr.c_i = apr.c_i(s_index);

end