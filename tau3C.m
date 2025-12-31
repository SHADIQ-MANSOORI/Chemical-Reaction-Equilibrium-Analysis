%This is a function file of tau for 3 component for NRTL equation

function t = tau3C(b,T)
    R = 1.987;
    for i = 1:3
        for j = 1:3
            t(i,j) = b(i,j)/(R*T);
        end
    end
end