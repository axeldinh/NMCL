function e = p_error(U_approx, U_ex, h, p)

if p == inf
    e = max(abs(U_approx - U_ex), [], 2);
else
    e = (h*sum(abs(U_approx-U_ex).^p, 2)).^(1/p);
end