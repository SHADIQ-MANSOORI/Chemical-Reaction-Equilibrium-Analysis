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

X1 = 0.3; %Given
X2 = 1-X1;
P = 101.33; %kPa , given
T1sat = Ti(A(1),B(1),C(1),P);
T2sat = Ti(A(2),B(2),C(2),P);

T_bbl = (X1*T1sat) + (X2*T2sat);
T_old = 0;
tol = 1e-11;
while true
    P1sat = Antonie(A(1),B(1),C(1),T_bbl);
    P2sat = Antonie(A(2),B(2),C(2),T_bbl);
    alpA = P1sat/P2sat;
    tau12 = tau(12,b,T_bbl);
    tau21 = tau(21,b,T_bbl);
    tau_cal = [tau12 tau21];
    G12 = gibbs(12,alp,tau_cal);
    G21 = gibbs(21,alp,tau_cal);
    G = [G12 G21];
    g1 = Gamma_NRTL(1,X1,tau_cal,G);
    g2 = Gamma_NRTL(2,X1,tau_cal,G);
    Psat = P/((X1*g1)+((X2*g2)/alpA));
    T_bbl = Ti(A(1),B(1),C(1),Psat);
    if abs(T_bbl-T_old) < tol
        break;
    else
        T_old = T_bbl;
    end
end

fprintf("BUbble Temp = %.3f K\n",T_bbl);

%---------Part B---------------------
g1 = 1;
g2 = 1;
T_dew = (X1*T1sat) + (X2*T2sat);
T_old = 0;
Y1 = 0.3;
Y2 = 1-Y1;
while true
    P1sat = Antonie(A(1),B(1),C(1),T_dew);
    P2sat = Antonie(A(2),B(2),C(2),T_dew);
    alpA = P1sat/P2sat;
    tau12 = tau(12,b,T_dew);
    tau21 = tau(21,b,T_dew);
    tau_cal = [tau12 tau21];
    G12 = gibbs(12,alp,tau_cal);
    G21 = gibbs(21,alp,tau_cal);
    G = [G12 G21];
    Psat = P*((Y1/g1)+((Y2*alpA)/g2));
    T_dew = Ti(A(1),B(1),C(1),Psat);
    x1 = (Y1*P)/(g1*Psat);
    g1 = Gamma_NRTL(1,x1,tau_cal,G);
    g2 = Gamma_NRTL(2,x1,tau_cal,G);
    if abs(T_dew-T_old) < tol
        break;
    else
        T_old = T_dew;
    end
end
fprintf("Dew Temp = %.3f K\n",T_dew);

%-------Part C--------------------
T = (T_dew+T_bbl)/2;
Z1 = 0.3;
Z2 = 1-Z1;
tol = 1e-6;
g1 = 1;
g2 = 1;
y_old = 0;
x_old = 0;
V_old = 0;
tau12 = tau(12,b,T);
tau21 = tau(21,b,T);
tau_cal = [tau12 tau21];
G12 = gibbs(12,alp,tau_cal);
G21 = gibbs(21,alp,tau_cal);
G = [G12 G21];
options = optimoptions('fsolve','Display','off');
while true
    
    K1 = (g1*P1sat)/P;
    K2 = (g2*P2sat)/P;

    f = @(V) ((Z1*K1)/(1+(V*(K1-1)))) + ((Z2*K2)/(1+(V*(K2-1)))) - 1;
    V_sol = fsolve(f,0.5,options);
    y1 = (Z1*K1)/(1+(V_sol*(K1-1)));
    x1 = y1/K1;
    g1 = Gamma_NRTL(1,x1,tau_cal,G);
    g2 = Gamma_NRTL(2,x1,tau_cal,G);
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

%-------Part D----------------------
Tb1 = Ti(A(1),B(1),C(1),P);
Tb2 = Ti(A(2),B(2),C(2),P);
P1sat = Antonie(A(1),B(1),C(1),Tb2);
P2sat = Antonie(A(2),B(2),C(2),Tb1);
tau12 = tau(12,b,Tb2);
tau21 = tau(21,b,Tb1);
tau_cal = [tau12 tau21];
G12 = gibbs(12,alp,tau_cal);
G21 = gibbs(21,alp,tau_cal);
G = [G12 G21];
alp12 = (Gamma_NRTL(1,0,tau_cal,G)*P1sat)/(Gamma_NRTL(2,0,tau_cal,G)*P2sat);
alp21 = (Gamma_NRTL(1,1,tau_cal,G)*P1sat)/(Gamma_NRTL(2,1,tau_cal,G)*P2sat);
flag = NaN;
if ( (alp12 < 1 && 1 < alp21) || (alp21 < 1 && 1 < alp12) )
    flag = 1;% 1 lie b/w alp12 and alp21
else
    flag = 0;
end

if flag == 1
    func = @(x) ...
        log(Gamma_NRTL(1,x,tau_cal,G) / Gamma_NRTL(2,x,tau_cal,G)) ...
        - log(P2sat / P1sat);
    Xaz = fsolve(func, 0.5, options);
    g_az = Gamma_NRTL(1,Xaz,tau,G);
    Paz = g_az*P1sat;
    
else
    fprintf("d. Azeotrope does not exist at the given condition\n");
end
