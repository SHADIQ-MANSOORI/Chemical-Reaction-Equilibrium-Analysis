%%
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

X1 = 0.3;
X2 = 0.4;
X3 = 1-X1-X2;
T = 338.15;
NoOfCompo = 3;
Tau = tau3C(b,T);
G = gibbs3C(a,Tau);

g = Gamma_NRTL3C(Tau,G,X1,X2,X3);
P1sat = Antonie(A(1),B(1),C(1),T);
P2sat = Antonie(A(2),B(2),C(2),T);
P3sat = Antonie(A(3),B(3),C(3),T);

P_bbl = (g(1)*X1*P1sat) + (g(2)*X2*P2sat) + (g(3)*X3*P3sat);
fprintf("a. Bubble Pressure = %.3f kPa \n",P_bbl);

%------Part B----------------
Y1 = 0.3;
Y2 = 0.4;
Y3 = 1-Y1-Y2;
tol = 0.01;
g = [1 1 1];
P_old = 0;
while true
    P_dew = 1/((Y1/(g(1)*P1sat))+(Y2/(g(2)*P2sat))+(Y3/(g(3)*P3sat)));
    X1 = (Y1*P_dew)/(g(1)*P1sat);
    X2 = (Y2*P_dew)/(g(2)*P2sat);
    g = Gamma_NRTL3C(Tau,G,X1,X2,X3);
    if abs(P_dew-P_old) < tol
        break;
    else
        P_old = P_dew;
    end
end
fprintf("b. Dew Pressure = %.3f KPa\n",P_dew);
%-----Part C---------------
Z1 = 0.3;
Z2 = 0.4;
Z3 = 1-Z1-Z2;
P = (P_bbl+P_dew)/2;
tol = 1e-10;
g = [1 1 1];
y_old = 0;
x_old = 0;
V_old = 0;
options = optimoptions('fsolve','Display','off');
i = 0;
while true
    
    K1 = (g(1)*P1sat)/P;
    K2 = (g(2)*P2sat)/P;
    K3 = (g(3)*P3sat)/P;

    f = @(V) ((Z1*K1)/(1+(V*(K1-1)))) + ((Z2*K2)/(1+(V*(K2-1)))) + ((Z3*K3)/(1+(V*(K3-1)))) - 1;
    V_sol = fsolve(f,0.5,options);
    y1 = (Z1*K1)/(1+(V_sol*(K1-1)));
    y2 = (Z2*K2)/(1+(V_sol*(K2-1)));
    x1 = y1/K1;
    x2 = y2/K2;
    g = Gamma_NRTL3C(Tau,G,X1,X2,X3);
    if abs(V_old-V_sol)<tol && abs(y_old-y1)<tol && abs(x_old-x1)<tol
        break;
    else
        y_old = y1;
        x_old = x1;
        V_old = V_sol;
    end
    i = i+1;
end

fprintf("c. Flash Cal Output-----------------------------------------------------------------------------------\n");
fprintf("\t x1 = %.3f\t X2 = %.3f\t X3 = %.3f\t Y1 = %.3f\t Y2 = %.3f\t Y3 = %.3f\t V = %.3f\t L = %.3f \n",x1,x2,1-x1-x2,y1,y2,1-y1-y2,V_sol,1-V_sol);
fprintf("-------------------------------------------------------------------------------------------------------\n");


