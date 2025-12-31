%This is funtion file for Del for i component

function d = deliC(Vi,Vj,aij,T)
    R = 1.987;
    d = (Vi/Vj)*exp(-aij/(R*T));
end