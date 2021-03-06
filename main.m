%% Laminar Viscous Solver
% name    : LVS 2D Airfoil 
% author  : Irsyad L, Ghifari A. F, Rashid
% date    : December 2018
% version : 1.0

%% clearing display and variables 
clc; clear all; 
fprintf('\t \t Laminar viscous solver for 2D airfoil\n')
fprintf('\t \t Starting program...\n')
errordelta = 100;

%% Main Program
U = input('Input Freestream Velocity :');
aoa = input('Input Angle of Attack :');
Tfs = input('Input Freestream Temperature (K) :');
Twall = input('Input Airfoil Surface Temperature (K) :');
load naca2410.txt;
%Reverse indexing, panel 1 begin from lower section trailing edge
foilcoord = flip(naca2410);
Xb = foilcoord(:,1)'; 
Yb = foilcoord(:,2)'; 
M = length(Xb)-1;
G = zeros(1,M+1);
delta = zeros(M,1);
transp1 = 0;
transp2 = 0;
i=1;

while errordelta >= 1e-05
    deltaimin1 = delta; %thickness at previous iter
    [Vtan,X,Y,Cp,Xb,Yb] = VortexPanelMethod(aoa,G,Xb,Yb); %solve vortex panel
    [Cf1,Cf2,delta,G,YBL,transp1,transp2,su,sl] = boundarylayer(U,Vtan,X,Y); %solve BL
    errordelta = sum(abs(delta-deltaimin1)./delta); %delta(difference) of delta star (BL thickness)
    if i == 1
        CP1 = Cp;
    end
    disp(errordelta);
    i=i+1;
end
[qwall,delc,Ythermal] = thermal(Twall,Tfs,U,Vtan,X,Y);

xtransp1 = X(transp1); ytransp1 = Y(transp1);
xtransp2 = X(transp2); ytransp2 = Y(transp2);
% %% Eliminate Cp Errors at trailing edge
% %Lower Cp
%     for i = 1:M/2
%         if abs(Cp(i))> 1 
%             Cp(i) = CP1(i);
%         end
%     end
% %Upper Cp
%     for i = 1+M/2 : M
%         if abs(Cp(i))> 1 
%             Cp(i) = CP1(i);
%         end
%     end
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
CpPlot = -Cp;
plot(Xb,Yb)
plot(X,CpPlot,'-o')
ylim([-1,1])
xlabel('X')
ylabel('-Cp')
grid on

figure(2)
plot(X,Y,'k',X,YBL,'b',X,Ythermal,'r',xtransp1,ytransp1,'x',xtransp2,ytransp2,'^')
legend('airfoil','momentum boundary layer', 'thermal boundary layer','transition 1','transition 2')
title('Boundary Layer')
axis([0 1 -0.3 0.3]);
xlabel('X')
ylabel('Y')            
grid on; 

figure(3)
plot(Xb,Yb)
plot(X,qwall,'-o')
xlabel('X')
ylabel('q wall (W/m2)')
grid on