%This is function file of gamma function for NRTL

function g = Gamma_NRTL(k,X1,tau,G)
    X2 = 1 - X1;
    if k == 1
        lng = X2^2 * (tau(2)*(G(2)/(X1 + X2*G(2)))^2 + tau(1)*G(1)/(X2 + X1*G(1))^2);
    elseif k == 2
        lng = X1^2 * (tau(1)*(G(1)/(X2 + X1*G(1)))^2 + tau(2)*G(2)/(X1 + X2*G(2))^2);
    else
        error("k must be 1 or 2");
    end
    g = exp(lng);
end
