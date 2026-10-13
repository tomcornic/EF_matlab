function Be = RHS(xa,xb,q,m,M)
% Vecteur elementaire 4x1 (membre de droite l(.)) sur l'element [xa,xb]
% Ordre : (w_a, theta_a, w_b, theta_b)
%   l((v,eta)) = Int q v + m eta dx   (+ termes de bord ajoutes dans assemblage)
h   = xb - xa;
phi = {@(x) (xb-x)/h , @(x) (x-xa)/h};
Be = zeros(4,1);
for i = 1:2
    Be(2*i-1) = quad(@(x) q(x).*phi{i}(x), xa, xb, M);
    Be(2*i)   = quad(@(x) m(x).*phi{i}(x), xa, xb, M);
end
end
