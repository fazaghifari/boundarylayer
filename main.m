%% Laminar Viscous Solver
% name    : LVS 2D Airfoil 
% author  : Irsyad L, Ghifari A. F, Rashid
% date    : December 2018
% version : 1.0

%% clearing display and variables 
clc; clear all; 
fprintf('\t \t Laminar viscous solver for 2D airfoil\n')
fprintf('\t \t Starting program...\n')
deltadelta = 100;

%% Main Program
U = input('Input Freestream Velocity :');
aoa = input('Input Angle of Attack :');
load naca1408.txt;
%Reverse indexing, panel 1 begin from lower section trailing edge
foilcoord = flip(naca1408);
Xb = foilcoord(:,1)'; 
Yb = foilcoord(:,2)'; 
M = length(Xb)-1;
G = zeros(1,M+1);
delta = zeros(M,1);
transp1 = 0;
transp2 = 0;
i=1;

while i < 10
    deltaimin1 = delta; %thickness at previous iter
    [Vtan,X,Y,Cp,Xb,Yb] = VortexPanelMethod(aoa,G,Xb,Yb); %solve vortex panel
    [Cf1,Cf2,delta,G,YBL,transp1,transp2,su,sl] = boundarylayer(U,Vtan,X,Y); %solve BL
    deltadelta = sum(abs(delta-deltaimin1)./delta)*100; %delta(difference) of delta star (BL thickness)
    if i == 1
        CP1 = Cp;
    end
    disp(deltadelta);
    i=i+1;
end
%% Eliminate Cp Errors at trailing edge
%Lower Cp
    for i = 1:M/2
        if abs(Cp(i))> 1 
            Cp(i) = CP1(i);
        end
    end
%Upper Cp
    for i = 1+M/2 : M
        if abs(Cp(i))> 1 
            Cp(i) = CP1(i);
        end
    end
%% Calculate Cl and Cd
CP_u = 0;CP_l=0;
CF_u = 0;CF_l=0;
%Coefficient of lift
%upper
for i=su:M-1
    CP_u = CP_u + (Cp(i+1)+Cp(i))*(X(i+1)-X(i))/2;  
end
%lower
for i=sl:-1:2
    CP_l = CP_l + (Cp(i-1)+Cp(i))*(X(i-1)-X(i))/2;
end
Cl = (CP_l - CP_u)*cos(aoa);
%Coefficient of drag
%upper
for i=su:transp1-1
    CF_u = CF_u + (Cf1(i+1)+Cf1(i))*(X(i+1)-X(i))/2;
end
%lower
for i=sl:-1:transp2+1
    CF_l = CF_l + (Cf1(i-1)+Cf1(i))*(X(i-1)-X(i))/2;
end
Cd = CF_u + CF_l ;

%% Plotting
disp('Lift and Drag Coefficients:')
disp('Cl')
disp(Cl)
disp('Cd')
disp(Cd)

figure(1)
CpPlot = -Cp./max(Cp);
plot(Xb,Yb)
plot(X,CpPlot,'-o')
ylim([-1,1])
grid on

figure(2)
plot(X,Y,X,YBL)
legend('airfoil','boundary layer')
title('Boundary Layer')
axis([0 1 -0.3 0.3]);
xlabel('X')
ylabel('Y')            
grid on; 
