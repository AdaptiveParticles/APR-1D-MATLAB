function grad_k = get_grad_k(y,f,k,r)
%
%   Bevan Cheeseman 2017
%
%   Calculate derivative of order k
%

order = k+r-1;
dim = 1;

dir = k;
%% set up

num_points = nchoosek(dim + order,dim);


comb = zeros(1,dim);

for num = 1:order
    temp = zeros(1,dim);
    
    comb = get_set(temp,comb,num,dim);
end

comb_u = unique(comb,'rows');

b = get_dirivative_b(dim,order,dir);

%coeffs = get_coeffs(x_p',x_q,b,comb_u,1);

grad_k = zeros(length(y),1);

offset_m = round(num_points/2);
offset_p = num_points - offset_m;

y = y(:);
f = f(:);

N = length(y);

for i = 1:length(y)
  
    min_y = max(1,i-offset_m);  
    max_y = min_y + num_points -1;
  
    if((i + offset_p) > N)
        max_y = N;
        min_y = N - (num_points-1);
    end
  
    x_q = y(i);
    x_p = y(min_y:max_y);
    
    coeffs = get_coeffs(x_p,x_q,b,comb_u,1,k);
    
    grad_k(i) = sum(coeffs'.*f(min_y:max_y));
    
end



end
function coeffs = get_coeffs(x_p,x_q,b,comb_u,dim,k)

x_pq = (x_p-repmat(x_q,length(x_p),1));

dist = sqrt(sum((x_pq).^2,2));

h = max(dist);
eps = h;

num_terms = size(comb_u,1);
num_points = size(x_p,1);

Vg = ones(num_terms,num_points);

for i = 1:num_terms
    for j = 1:num_points
        for d = 1:dim
           Vg(i,j) = Vg(i,j)*(x_pq(j,d)/h).^comb_u(i,d);   
        end 
    end
end

V = Vg';


W=-(abs(dist).^2)/(2*(eps^2));
%Wd= -(abs(dist/eps).^2);
E=diag(exp(W));
B=E*V;
A=B'*B;

eps = max(dist);

Wd= -(abs(dist/eps).^2);

aT=A\b;

coeffs=(1/h^k)*(V*aT.*exp(Wd))';


end