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

T = 333.15; % K 
tau12 = tau(12,b,T);
tau21 = tau(21,b,T);
tau = [tau12 tau21];
G12 = gibbs(12,alp,tau);
G21 = gibbs(21,alp,tau);
G = [G12 G21];
P1sat = Antonie(A(1),B(1),C(1),T);
P2sat = Antonie(A(2),B(2),C(2),T);

X1 = 0:0.01:1;

for i=1:length(X1)
    x1 = X1(i);
    g1 = Gamma_NRTL(1,x1,tau,G);
    g2 = Gamma_NRTL(2,x1,tau,G);
    P(i) = Pb(x1,g1,g2,P1sat,P2sat);
    y(i) = (x1*P1sat*g1)/P(i);
end

figure;
plot(X1,P,'LineWidth',2);
hold on;
plot(y,P,'LineWidth',2);
hold off;
xlabel("X,Y------------>");
ylabel("P------------->");
legend('P-x','P-y',Location='northwest');
title("P-x-y diagram");
grid on;
grid minor;