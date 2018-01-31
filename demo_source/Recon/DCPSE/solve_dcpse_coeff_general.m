function out = solve_dcpse_coeff_general(x_p,order,dim,x_q,varargin)
%
%
%   Bevan Cheeseman 2016
%
%   x_p local coordinates of points
%
%   Computers Vandermonde of any order or dimension
%
%
%

%% Compute the monomial basis

if isempty(varargin)

    comb = zeros(1,dim);

    for num = 1:order
        temp = zeros(1,dim);
    
        comb = get_set(temp,comb,num,dim);   
    end

    comb_u = unique(comb,'rows');

else
    comb_u = varargin{1}; 

end

num_terms = size(comb_u,1);
num_points = size(x_p,1);

%% Construct the Vandermonde Matrix

x_pq = (x_p-repmat(x_q,length(x_p),1));

dist = sqrt(sum((x_pq).^2,2));

h = max(dist);
eps = h;

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
Wd= -(abs(dist/eps).^2);
E=diag(exp(W));
B=E*V;
A=B'*B;

out.condition_numberA=condest(A);
sz = size(V);
if(any(sz~=sz(1)))
    out.condition_numberV=0;
else
    out.condition_numberV=condest(V);
end
out.A = A;
out.comb_u = comb_u;
out.V = V;


end

