function [x,A,B,w,th] = assemblage(N,L,EI,S,q,m,FL,ML,M,reduit)
% Construit A et B, impose les CL, resout A U = B.
% N noeuds interieurs -> N+2 noeuds, N+1 elements, h = L/(N+1)  (comme en cours)
% Inconnue : U = (w_1,th_1, w_2,th_2, ..., w_{N+2},th_{N+2})'
% CL : encastrement en x=0 (Dirichlet), effort FL et moment ML en x=L (Neumann)
h = L/(N+1);
x = (0:N+1)*h;
A = zeros(2*(N+2));
B = zeros(2*(N+2),1);
for e = 1:N+1                          % boucle sur les elements
    idx = 2*(e-1) + (1:4);             % ddl de l'element e
    A(idx,idx) = A(idx,idx) + LHS(x(e),x(e+1),EI,S,M,reduit);
    B(idx)     = B(idx)     + RHS(x(e),x(e+1),q,m,M);
end
% Conditions de Neumann (naturelles) : ajout au second membre
B(end-1) = B(end-1) + FL;              % S(w'-th)(L) = FL
B(end)   = B(end)   + ML;              % EI th'(L)   = ML
% Conditions de Dirichlet en x=0 : w(0)=0, th(0)=0
A([1 2],:) = 0;  A(:,[1 2]) = 0;
A(1,1) = 1;  A(2,2) = 1;
B([1 2]) = 0;
% Resolution
U  = A\B;
w  = U(1:2:end);
th = U(2:2:end);
end