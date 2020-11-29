function U = ref_to_current(U_ref, x_ref, xc)
% returns an array U which contains values of U_ref rescaled on xc

U = zeros(2,length(xc));

for i = 1:length(xc)
    j = find(abs(xc(i)-x_ref) == min(abs(xc(i)-x_ref)));
    j = j(1); % In case xc ends up exactly between two values of x_ref
    U(:,i) = U_ref(:,j);
end

