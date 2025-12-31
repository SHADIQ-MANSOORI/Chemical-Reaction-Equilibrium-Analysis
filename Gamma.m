%This is a function file to caalculate Gamma , Activity Coeff for wilson.

function g = Gamma(k,X1,del12,del21)

    X2 = 1-X1;
    term1 = -log(X1+(X2*del12));
    term2 = -log(X2+(X1*del21));
    term3 = del12/(X1+(X2*del12));
    term4 = del21/(X2+(X1*del21));
    if k==1
        g = exp(term1+(X2*(term3-term4)));
    elseif k==2
        g = exp(term2-(X1*(term3-term4)));
    else
        error('Only k = 1 or k = 2 is acceptable for gamma1 and gamma2');
    end

end