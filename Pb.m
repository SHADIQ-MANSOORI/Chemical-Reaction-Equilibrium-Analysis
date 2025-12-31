% This is a function file for  Bubble Pressure

function p = Pb(X1,g1,g2,P1sat,P2sat)
    X2 = 1-X1;
    p = (X1*g1*P1sat)+(X2*g2*P2sat);
end