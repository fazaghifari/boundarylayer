clc
clear
hold on
grid on
%axis([0 1 -1 1])
% Xb = [1,0.933,0.75,0.5,0.25,0.067,0,0.067,0.25,0.5,0.75,0.933,1]
% Yb = [0,-0.005,-0.017,-0.033,-0.042,-0.033,0,0.045,0.076,0.072,0.044,0.013,0]
%Input Yb dari JavaFoil harus dibalik
%Jangan lupa ubah M
% Xb = [0,0.067,0.25,0.5,0.75,0.933,1];
% Yb = [0,0.045,0.076,0.072,0.044,0.013,0];
XbIn = fscanf(fopen('Xbody.txt'),'%f');
YbIn = fscanf(fopen('Ybody.txt'),'%f');
Xb = XbIn';
Yb = YbIn';

M = length(Xb)-1;
MP1 = M + 1;
Alpha = 0 * pi/180;


for i=1:M
    ip1 = i + 1;
    X(i) = 0.5 * (Xb(i) + Xb(ip1));
    Y(i) = 0.5 * (Yb(i) + Yb(ip1));
    S(i) = sqrt ((Xb(ip1)-Xb(i))^2 + (Yb(ip1)-Yb(i))^2 );
    Theta(i) = atan2(Yb(ip1)-Yb(i), Xb(ip1)-Xb(i));
    Sine(i) = sin (Theta(i));
    Cosine(i) = cos (Theta(i));
    RHS(i) = sin (Theta(i) - Alpha);
    RHSC(i) = cos (Theta(i) - Alpha);
end
i=0;

for i=1:M
    for j=1:M
        if (i==j)
            CN1(i,j) = -1;
            CN2(i,j) = 1;
            CT1(i,j) = 0.5*pi;
            CT2(i,j) = 0.5*pi;
        else
            A = -(X(i)-Xb(j))*Cosine(j) - (Y(i)-Yb(j))*Sine(j);
            B = (X(i)-Xb(j))^2 + (Y(i)-Yb(j))^2;
            C = sin(Theta(i)-Theta(j));
            D = cos(Theta(i)-Theta(j));
            E = (X(i)-Xb(j))*Sine(j) - (Y(i) - Yb(j))*Cosine(j);
            F = log (1 + S(j)*(S(j)+2*A)/B);
            G = atan2 (E*S(j), B+A*S(j));
            P = (X(i)-Xb(j))*sin(Theta(i)-2.*Theta(j))+(Y(i)-Yb(j))*cos(Theta(i)-2.*Theta(j));
            Q = (X(i)-Xb(j))*cos(Theta(i)-2.*Theta(j))-(Y(i)-Yb(j))*sin(Theta(i)-2.*Theta(j));
            CN2(i,j) = D + 0.5*Q*F/S(j) - (A*C+D*E)*G/S(j);
            CN1(i,j) = 0.5*D*F + C*G - CN2(i,j);
            CT2(i,j) = C + 0.5*P*F/S(j) +(A*D-C*E)*G/S(j);
            CT1(i,j) = 0.5*C*F - D*G -CT2(i,j);
        end
    end
end

for i=1:M
    AN(i,1) = CN1(i,1);
    AN(i,MP1) = CN2(i,M);
    AT(i,1) = CT1(i,1);
    AT(i,MP1)=CT2(i,M);
    for j=2:M
        AN(i,j) = CN1(i,j)+CN2(i,j-1);
        AT(i,j) = CT1(i,j)+CT2(i,j-1);
    end
    AN(MP1,1) = 1;
    AN(MP1,MP1) = 1;
    for j=2:M
        AN(MP1,j) = 0;
    end
    RHS(MP1) = 0;
end

Gam = inv(AN)*RHS';
Vtan = RHSC' + AT*Gam; Vplot = Vtan./max(Vtan);
Cp = 1-Vtan.^2;
CpPlot = -Cp./max(Cp);
plot(Xb,Yb)
plot(X,Vplot)
plot(X,CpPlot)
Xp = X';