function E = erreurL2(x,w,th,wex,thex,M)
% Norme L2 de l'erreur (w,th) - (wex,thex), les fonctions EF etant affines par morceaux
E2 = 0;
for e = 1:length(x)-1
    xa = x(e); xb = x(e+1); h = xb-xa;
    fw = @(t) ( w(e)*(xb-t)/h +  w(e+1)*(t-xa)/h - wex(t)).^2;
    ft = @(t) (th(e)*(xb-t)/h + th(e+1)*(t-xa)/h - thex(t)).^2;
    E2 = E2 + quad(fw,xa,xb,M) + quad(ft,xa,xb,M);
end
E = sqrt(E2);
end
