%This is a function file to calculate Gibbs Energy for NRTL

function G = gibbs(k,alp,tau)

    if k==12
        G = exp(-alp*tau(1));
    elseif k==21
        G = exp(-alp*tau(2));
    else
        error("k should be equal to 12 for G12 or 21 for G21");
    end
end