close all;
clear;
clc;

%Acetone(1) , Water(1)
A = [14.3916 16.2630];
B = [2795.82 3799.89];
C = [-43.15 -46.80];

%NRTL Parameter
V = [0.07405 0.01807];
b = [2642.1 5013.3]; % [b12 b21]
alp = 0.5343;

P = 101.33; % KPa 
T1sat = Ti(A(1),B(1),C(1),P);
T2sat = Ti(A(2),B(2),C(2),P);
T = (T1sat+T2sat)/2;
X1 = 0:0.01:1;

for i=1:length(X1)
    x1 = X1(i);
    x2 = 1-x1;
    P1sat = Antonie(A(1),B(1),C(1),T);
    P2sat = Antonie(A(2),B(2),C(2),T);
    alpA = P1sat/P2sat;
    tau12 = tau(12,b,T);
    tau21 = tau(21,b,T);
    tau_clc = [tau12 tau21];
    G12 = gibbs(12,alp,tau_clc);
    G21 = gibbs(21,alp,tau_clc);
    G = [G12 G21];
    g1 = Gamma_NRTL(1,x1,tau_clc,G);
    g2 = Gamma_NRTL(2,x1,tau_clc,G);
    Psat = P/((g1*x1)+((g2*x2)/alpA));
    T = Ti(A(1),B(1),C(1),Psat);
    if T<T2sat && T>T1sat 
        Tf(i) = T;
        x(i) = x1;
        y(i) = (g1*Psat*x1)/P;
    end
end
figure;
Tf = [T2sat Tf(2:end-1) T1sat];
Tf = smooth(Tf, 3);    % moving average smoothing
y  = smooth(y, 3);
plot(x,Tf,'LineWidth',2);
hold on;
plot(y,Tf,'LineWidth',2);
hold off;
xlabel("X,Y------------>");
ylabel("T------------->");
legend('T-x','T-y');
title("T-x-y diagram");
grid on;
grid minor;

