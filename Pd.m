% This is the function file for Dew Pressure Function

function P = Pd(Y1,g1,g2,P1sat,P2sat)
    Y2 = 1-Y1;
    term1 = Y1/(g1*P1sat);
    term2 = Y2/(g2*P2sat);
    P = 1/(term1+term2);
end