function empty = check_descendants(i,l,L,l_max)
%
%   Bevan Cheeseman 2017
%
%   Take a particle cell c_{i,l} and check its descendants too see if any
%   are occupied in L
%
%   recursively defined function
%
%   Returns value empty = true, if there are no descendents in L, false otherwise
%

%get children
empty = true;

if(l == l_max)
    empty = true;
else
    [l_c,i_c] = get_children(i,l);
    
    for c = 1:length(l_c)
        temp  = recursive_children(empty,i_c(c),l_c(c),L,l_max);
        
        if(~temp)
            empty = false;
            break;
        end
    end
end

end
function empty = recursive_children(empty,i,l,L,l_max)
%
%   Bevan Cheeseman 2017
%
%   Recursive Function for checking children
%

check = L{l}(i+1);

if(check == 1)
    % in L
    empty = false;
else
    if l ~= l_max
        [l_c,i_c] = get_children(i,l);
        
        for c = 1:length(l_c)
            temp  = recursive_children(empty,i_c(c),l_c(c),L,l_max);
            
            if(~temp)
                empty = false;
                break;
            end
        end
        
    else
        empty = true;
    end
end


end