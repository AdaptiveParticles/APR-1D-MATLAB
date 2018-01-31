function b = get_dirivative_b(dim,order,dir_ind)
%
%
%   Bevan Cheeseman 2016
%
%   Input dim, order, and dirivative index {1,0,0} for x diff
%
%



comb = zeros(1,dim);

for num = 1:order
   temp = zeros(1,dim);
    
   comb = get_set(temp,comb,num,dim);   
end

comb_u = unique(comb,'rows');

indx = ones(size(comb_u(:,1)));

for i = 1:dim
    indx = (comb_u(:,i) == dir_ind(i)).*indx;
end


co_eff = 1;

for i = 1:dim
    co_eff = co_eff*factorial(dir_ind(i)); 
end


b = zeros(size(indx));
%b(indx ==1) = co_eff*(-1).^(sum(dir_ind));
b(indx ==1) = co_eff;

end