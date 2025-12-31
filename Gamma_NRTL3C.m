%This is function file of gamma for 3 component for NRTL

function gamma = Gamma_NRTL3C(Tau, G, x1, x2, x3)

x = [x1, x2, x3];
gamma = zeros(1,3);

for i = 1:3
    % ----- first summation term -----
    sum1 = 0;
    for j = 1:3
        denom1 = 0;
        for k = 1:3
            denom1 = denom1 + x(k) * G(k,i);
        end
        sum1 = sum1 + x(j) * G(j,i) * Tau(j,i) / denom1;
    end

    % ----- second summation term -----
    sum2 = 0;
    for j = 1:3
        denom2 = 0;
        for k = 1:3
            denom2 = denom2 + x(k) * G(k,i);
        end

        % Inner term: Σ_m x_m G_mj τ_mj / Σ_n x_n G_nj
        numer_inner = 0;
        denom_inner = 0;
        for m = 1:3
            numer_inner = numer_inner + x(m) * G(m,j) * Tau(m,j);
            denom_inner = denom_inner + x(m) * G(m,j);
        end

        sum2 = sum2 + (x(j) * G(i,j) / denom2) * (Tau(i,j) - numer_inner / denom_inner);
    end

    % final log gamma for component i
    ln_gamma = sum1 + sum2;
    gamma(i) = exp(ln_gamma);
end
end