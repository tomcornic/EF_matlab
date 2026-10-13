function pm = quadpm(F,a,b,M)
%qua Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    F
    a
    b
    M
end

arguments (Output)
    
    pm
end

dxt = (b - a)/M;            % pas de la subdivision
    pm = 0;

    for j = 0:M-1
        xt = a + j*dxt;         % extrémité gauche du sous-intervalle j
        pm = pm + F(xt + dxt/2);  % on évalue f au milieu
    end

    pm = pm*dxt;

end