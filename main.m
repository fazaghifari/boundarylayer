%% Laminar Viscous Solver
% name    : LVS 2D Airfoil 
% author  : Irsyad L, Ghifari A. F, Rashid
% date    : December 2018
% version : 1.0

%% clearing display and variables 
clc; clear all; 
fprintf('\t \t Laminar viscous solver for 2D airfoil\n')
fprintf('\t \t Starting program...\n')
G = 0;

%% Main Program
U = input('Input Freestream Velocity :');
aoa = input('Input Angle of Attack :');
[Vtan,X,Y,Cp,Xb,Yb] = VortexPanelMethod(aoa,G);
[Cf,delta,G] = boundarylayer(U,Vtan,X,Y);

%% Plotting
CpPlot = -Cp./max(Cp);
plot(Xb,Yb)
plot(X,CpPlot,'-o')
grid on
