function rg = quadrg(F,a,b,M)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    F
    a
    b
    M
end

arguments (Output)
    
    rg
end

dxt = (b - a)/M;
rg = 0;
for j = 0:M-1
    xt = a + j*dxt;
    rg = rg + F(xt);
end
rg = rg*dxt;

end