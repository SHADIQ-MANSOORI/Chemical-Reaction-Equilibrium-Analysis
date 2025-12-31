%This is a Funtion file for del

function d = del(k,V,a,T)
    R = 8.314; %J/mol.K 
    if k == 1 %for del12
        d = (V(2)/V(1))*exp(-a/(R*T));
    elseif k == 2 %for del21
        d = (V(1)/V(2))*exp(-a/(R*T));
    end
end