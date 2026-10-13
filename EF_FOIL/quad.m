function S = quad(f,a,b,M)
% Integrale de f sur [a,b] par la methode des trapezes composes (M sous-intervalles)
% (meme fonction que celle vue en cours)
% ATTENTION : ce nom masque la fonction native "quad" de Matlab, c'est voulu
% pour rester coherent avec le cours.
dxt = (b-a)/M;
S = 0;
for j = 0:M-1
    xt = a + j*dxt;
    S = S + (f(xt) + f(xt+dxt))/2;
end
S = S*dxt;
end
