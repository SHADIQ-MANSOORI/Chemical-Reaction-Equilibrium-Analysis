close all;
clear;
clc;
R = 8.314; % J/mol.K
%Acetone(1) , Water(1)
A = [14.3916 16.2630]; 
B = [2795.82 3799.89];
C = [-43.15 -46.80];

%Wilson Parameter
V = [0.07405 0.01807];
a = [1219.5 6062.5]; % [a12 a21]

X1 = 0:0.01:1;
T = 333.15; % degree C
P1sat = Antonie(A(1),B(1),C(1),T);
P2sat = Antonie(A(2),B(2),C(2),T);
del12 = del(1,V,a(1),T);
del21 = del(2,V,a(2),T);
for i=1:length(X1)
    count(i) = i;
    x1 = X1(i);
    g1 = Gamma(1,x1,del12,del21);
    g2 = Gamma(2,x1,del12,del21);
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



