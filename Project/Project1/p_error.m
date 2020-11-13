function e = p_error(U_approx, U_ex, dx, p)

if p == inf
    e = max(abs(U_approx - U-ex), [], 2);
else
    e = (dx*sum(abs(U_approx-U_ex).^p, 2)).^(1/p);
end