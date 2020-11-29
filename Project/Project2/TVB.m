function psi = TVB(v,M,h);

psi = v(:,1); inds = find(psi > M*h^2);

if (size(v(:,1)>0))
    psi(inds) = MINMOD(v(inds,:));
end
