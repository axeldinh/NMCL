function U_reshape = fine_to_current(U_fine, dx_fine, dx, a, b)


x = a:dx:b;
x_fine = a:dx_fine:b;
N = length(x);
U_reshape = zeros(2,N);

for i = 1:N
    
    j = 1;
    
    while x_fine(j) < x(i)
        j = j+1;
    end
    
    U_reshape(:,i) = U_fine(:,j);
    
end