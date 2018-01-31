function comb = get_set(temp,comb,dis_total,dim)
%
%   Bevan Cheeseman 2016
%
%   Recursive Function for monomial co-efficient generation 
%
    temp_old = temp;
    temp = zeros(size(temp_old,1)*dim,dim);
    
    for k = 1:size(temp_old,1)
        for i = 1:dim
            temp((k-1)*dim+i,:) = temp_old(k,:);
            
            temp((k-1)*dim+i,i) = temp_old(k,i) + 1;
            
        end
    end
    
    if (dis_total == sum(temp(1,:)))
       %append the results to the combinatorial total
        for k = 1:length(temp)
            comb(end+1,:) = temp(k,:); 
        end
    else
        %otherwise keep going recursive
        comb = get_set(temp,comb,dis_total,dim);
    end

end