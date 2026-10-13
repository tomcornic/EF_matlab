function tpz = quadtpz(f,a,b,M)
%qua Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    f
    a
    b
    M
end

arguments (Output)
    
    tpz
end

dxt = (b - a)/M;
    tpz = 0;
    for j = 0:M-1
        xt = a + j*dxt;
        tpz = tpz + (f(xt) + f(xt + dxt))/2;
    end
    tpz = tpz*dxt;
end