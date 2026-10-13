function simp = quadsimp(f, a, b, M)
% QUADSIMP  Formule composée de Simpson pour approcher I(f) = int_a^b f(x) dx
%   f : handle de fonction (ex. @(x) sin(x))
%   a, b : bornes de l'intervalle
%   M : nombre de sous-intervalles (chacun utilise ses 2 extrémités et son milieu)

arguments (Input)
    f
    a
    b
    M
end
arguments (Output)
    simp
end

dxt = (b - a)/M;            % pas de la subdivision
simp = 0;

for j = 0:M-1
    xt = a + j*dxt;         % extrémité gauche du sous-intervalle j
    simp = simp + f(xt) + 4*f(xt + dxt/2) + f(xt + dxt);
end

simp = simp*dxt/6;
end