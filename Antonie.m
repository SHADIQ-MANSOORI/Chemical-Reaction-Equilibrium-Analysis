%This is a funtion file for Antonie Eq 

function Psat = Antonie(A,B,C,T)
    Psat = exp(A-(B/(T+C)));
end