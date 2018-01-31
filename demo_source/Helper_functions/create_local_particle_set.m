function L_pc = create_local_particle_set(y,L,s_dom,l_max,varargin)
%
%   Bevan Cheeseman 2017
%
%   Creating the local particle cell set from L(y)
%

%% Set up the set C of all possible particle cells




if(~isempty(varargin))
    
   type = varargin{1}; 
else
   type = 1;
end

if type > 1
    L_pc = cell((l_max-1),1);
    
    for l = 1:(l_max-1)
        L_pc{l} = zeros(2^l,1);
    end
    
    l_max = l_max -1;
    
else
    L_pc = cell((l_max),1);
    
    for l = 1:(l_max)
        L_pc{l} = zeros(2^l,1);
    end
    
end


[c_l,c_i] = calculate_particle_cell(y,L,s_dom,l_max,type);

%for one level down optimizations

c_l = max(1,c_l);

for j = 1:length(c_l)

    L_pc{c_l(j)}(c_i(j)+1) = 1;

end


end