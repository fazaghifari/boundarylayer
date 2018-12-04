%% Laminar Viscous Solver
% name    : LVS 2D Airfoil 
% author  : Irsyad L, Ghifari A. F, Rashid
% date    : December 2018
% version : 1.0

%% clearing display and variables 
clc; clear all; 
fprintf('\t \t Laminar viscous solver for 2D airfoil\n')
fprintf('\t \t Starting program...\n')
delta = 0;
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
G = zeros(1,M);
delta = zeros(M,1);

for i = 1:1
    deltaimin1 = delta;
    [Vtan,X,Y,Cp,Xb,Yb] = VortexPanelMethod(aoa,G,Xb,Yb);
    [Cf,delta,G,YBL,transp1,transp2] = boundarylayer(U,Vtan,X,Y);
    deltadelta = sum(abs(delta-deltaimin1)./delta)*100;
end


%% Plotting
figure(1)
CpPlot = -Cp./max(Cp);
plot(Xb,Yb)
plot(X,CpPlot,'-o')
grid on

figure(2)
plot(X,Y,X,YBL)
legend('airfoil','boundary layer')
title('Boundary Layer')
axis([0 1 -0.3 0.3]);
xlabel('X')
ylabel('Y')            
grid on; 
