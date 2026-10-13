function Ke = LHS(xa,xb,EI,S,M,reduit)
% Matrice elementaire 4x4 (membre de gauche a(.,.)) sur l'element [xa,xb]
% Ordre des ddl locaux : (w_a, theta_a, w_b, theta_b)
% Lignes = fonction test (v ou eta), colonnes = inconnue (w ou theta)
%   a((w,th),(v,eta)) = Int EI th' eta' + S (w'-th)(v'-eta) dx
% EI, S : handles @(x) ; M : nb de sous-intervalles du trapeze
% reduit = 1 : termes de cisaillement integres au point milieu (1 point)
h   = xb - xa;
phi  = {@(x) (xb-x)/h ,  @(x) (x-xa)/h };   % fonctions chapeau locales
dphi = {@(x) -1/h+0*x ,  @(x)  1/h+0*x };
Qb = @(f) quad(f,xa,xb,M);                  % flexion : trapezes
if reduit
    Qs = @(f) quadm(f,xa,xb,1);             % cisaillement : point milieu
else
    Qs = @(f) quad(f,xa,xb,M);
end
Ke = zeros(4,4);
for i = 1:2                 % noeud de la fonction test
    for j = 1:2             % noeud de l'inconnue
        pv = phi{i}; dpv = dphi{i};
        pu = phi{j}; dpu = dphi{j};
        r = 2*(i-1) + [1 2];
        c = 2*(j-1) + [1 2];
        Ke(r(1),c(1)) =  Qs(@(x) S(x).*dpu(x).*dpv(x));           % w   / v
        Ke(r(1),c(2)) = -Qs(@(x) S(x).*pu(x).*dpv(x));            % th  / v
        Ke(r(2),c(1)) = -Qs(@(x) S(x).*dpu(x).*pv(x));            % w   / eta
        Ke(r(2),c(2)) =  Qb(@(x) EI(x).*dpu(x).*dpv(x)) ...       % th  / eta
            + Qs(@(x) S(x).*pu(x).*pv(x));
    end
end
end