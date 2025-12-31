%This is a function file to calcualte gibbs enenrgy for NRTL 3 component

function G = gibbs3C(alp, Tau)
    for i = 1:3
        for j = 1:3
            G(i,j) = exp(-alp(i,j) * Tau(i,j));
        end
    end
end