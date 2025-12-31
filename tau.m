%This is a function file to calculate tau for NRTL

function t = tau(k,b,T)
    R = 8.314; %J/mol.K
    if k==12
        t = b(1)/(R*T);
    elseif k == 21
        t = b(2)/(R*T);
    else
        error("k should be equal to 12 for tau12 or 21 for tau21");
    end
end