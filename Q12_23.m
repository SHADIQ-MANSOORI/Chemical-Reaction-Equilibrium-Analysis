close all;
clear;
clc;

%Acetone(1) , Methanol(2) , Water(3)
A = [14.3916 16.5938 16.2620]; 
B = [2795.82 3644.30 3799.89];
C = [-43.15 -33.39 -46.80];

%NRTL Parameter
V = [0.07405 0.04073 0.01807];
a = [0 0.3084 0.5343;
     0.3084 0 0.2994;
     0.5343 0.2994 0];

b = [0 184.70 631.05;
    222.64 0 -253.88;
    1197.41 845.21 0];



X1 = 0.3; %Given
X2 = 0.4;
X3 = 1-X1-X2;
P = 101.33; %kPa , given
T1sat = Ti(A(1),B(1),C(1),P);
T2sat = Ti(A(2),B(2),C(2),P);
T3sat = Ti(A(3),B(3),C(3),P);
T_bbl = (X1*T1sat) + (X2*T2sat) + (X3*T3sat);
T_old = 0;
tol = 1e-11;
while true
    P1sat = Antonie(A(1),B(1),C(1),T_bbl);
    P2sat = Antonie(A(2),B(2),C(2),T_bbl);
    P3sat = Antonie(A(3),B(3),C(3),T_bbl);
    alpA = P1sat/P2sat;
    alpA_p = P1sat/P3sat;

    G = gibbs3C(a,tau3C(b,T_bbl));
    g = Gamma_NRTL3C(tau3C(b,T_bbl),G,X1,X2,X3);
    Psat = P/((X1*g(1))+((X2*g(2))/alpA)+((X3*g(3))/alpA_p));
    T_bbl = Ti(A(1),B(1),C(1),Psat);
    if abs(T_bbl-T_old) < tol
        break;
    else
        T_old = T_bbl;
    end
end
fprintf("a. BUbble Temp = %.3f K\n",T_bbl);
%---------Part B---------------------
g=[1 1 1];
T_dew = (X1*T1sat) + (X2*T2sat) + (X3*T3sat);
T_old = 0;
Y1 = 0.3;
Y2 = 0.4;
Y3 = 1-Y1-Y2;
while true
    P1sat = Antonie(A(1),B(1),C(1),T_dew);
    P2sat = Antonie(A(2),B(2),C(2),T_dew);
    P3sat = Antonie(A(3),B(3),C(3),T_dew);
    alpA = P1sat/P2sat;
    alpA_p = P1sat/P3sat;
    G = gibbs3C(a,tau3C(b,T_dew));
    Psat = P*((Y1/g(1))+((Y2*alpA)/g(2))+((Y3*alpA_p)/g(3)));
    T_dew = Ti(A(1),B(1),C(1),Psat);
    x1 = (Y1*P)/(g(1)*Psat);
    x2 = (Y2*P)/(g(2)*Psat);
    g = Gamma_NRTL3C(tau3C(b,T_bbl),G,x1,x2,X3);
    if abs(T_dew-T_old) < tol
        break;
    else
        T_old = T_dew;
    end
end
fprintf("b. Dew Temp = %.3f K\n",T_dew);
%-------Part C--------------------
T = (T_dew+T_bbl)/2;
Z1 = 0.3;
Z2 = 0.4;
Z3 = 1-Z1-Z2;
tol = 1e-6;
g = [1 1 1];
y_old = 0;
x_old = 0;
V_old = 0;
G = gibbs3C(a,tau3C(b,T));
options = optimoptions('fsolve','Display','off');
while true
    
    K1 = (g(1)*P1sat)/P;
    K2 = (g(2)*P2sat)/P;
    K3 = (g(3)*P3sat)/P;
    f = @(V) ((Z1*K1)/(1+(V*(K1-1)))) + ((Z2*K2)/(1+(V*(K2-1)))) + ((Z3*K3)/(1+(V*(K3-1)))) - 1;
    V_sol = fsolve(f,0.3,options);
    y1 = (Z1*K1)/(1+(V_sol*(K1-1)));
    y2 = (Z2*K2)/(1+(V_sol*(K2-1)));
    x1 = y1/K1;
    x2 = y2/K2;
    g = Gamma_NRTL3C(tau3C(b,T),G,x1,x2,X3);
    if abs(V_old-V_sol)<tol && abs(y_old-y1)<tol && abs(x_old-x1)<tol
        break;
    else
        y_old = y1;
        x_old = x1;
        V_old = V_sol;
    end
end

fprintf("c. Flash Cal Output-----------------------------------------------------------------------------------\n");
fprintf("\t x1 = %.3f\t X2 = %.3f\t Y1 = %.3f\t Y2 = %.3f\t V = %.3f\t L = %.3f \n",x1,1-x1,y1,1-y1,V_sol,1-V_sol);
fprintf("-------------------------------------------------------------------------------------------------------\n");
