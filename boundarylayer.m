function [Cf1,Cf2,delta,Gee,YBL,transp1,transp2,su,sl] = boundarylayer(U,Vtan,X,Y)
    %% ________________________________THWAITES__________________________________
    %Referensi : Dinamika Fluida (Lavi R Zuhal), Viscous Flow (Frank White),Intro
                % to theoretical and computational aerodynamics (Moran)

    miu = 0.0000181206; %Dynamics Viscosity
    rho = 1.225;        %Density
    nu = miu/rho;       %kinematics Viscosity
    Uin = Vtan*U;
    M = length(X);
    MP1 = M+1;
    %% .......................Calculate momentum thickness........................
    teta = zeros(M,1);
    % mencari titik stagnasi
    sl = find(Uin==min(Uin));
    su = sl;

    %upper airfoil
    teta(su) = sqrt(0.075*nu/(abs((Uin((su)+1)-Uin(su))/(X((su)+1)-X((su))))));
    for j=su+1:M
        integral1 = 0;
        for i=su+1:j
            integral1 = integral1 + (Uin(i)^5 + Uin(i-1)^5)*abs((X(i)-X(i-1)))/2 ;
        end
        teta(j) = sqrt(0.45*nu*(integral1)/(Uin(j)^6));
    end
    %lower airfoil
    teta(sl) = sqrt(0.075*nu/(abs((Uin((sl)-1)-Uin(sl))/(X((sl)-1)-X((sl))))));
    for j=sl-1:-1:1
        integral1 = 0;
        for i=sl-1:-1:j
            integral1 = integral1 + (Uin(i)^5 + Uin(i+1)^5)*abs((X(i)-X(i+1)))/2 ;
        end
        teta(j) = sqrt(0.45*nu*(integral1)/(Uin(j)^6));
    end


    %% ..............Calculate pressure gradient parameter (lamda)...............
    lamda = zeros(M,1);
    %upper
    lamda(su) = ((teta(su)^2)/nu)*(Uin(su+1)-Uin(su))/abs((X(su+1)-X(su)));
    for i=su+1:M
        lamda(i) = ((teta(i)^2)/nu)*(Uin(i)-Uin(i-1))/abs((X(i)-X(i-1)));
    end
    %lower
    lamda(sl) = ((teta(sl)^2)/nu)*(Uin(sl-1)-Uin(sl))/abs((X(sl-1)-X(sl)));
    for i=sl-1:-1:1
        lamda(i) = ((teta(i)^2)/nu)*(Uin(i)-Uin(i+1))/abs((X(i)-X(i+1)));
    end

    %% ............Calculate wall shear-stress and displacement thickness.........
    L = zeros(M,1);
    S = zeros(M,1);
    H = zeros(M,1);
    Cf1 = zeros(M,1); %Obtained using l
    Cf2 = zeros(M,1); %Obtained using tau wall, Cf1 and Cf2 have the same physical meaning, only different in the way to obtain them
    tauw = zeros(M,1);
    delta = zeros(M,1);
    YBL = zeros(M,1);
    %shear stress at wall
    %upper
    for i = su:M
        z = 0.25 - lamda(i);
        if lamda(i)<0.1 && lamda(i)>0
            L(i) = 0.22 + 1.57*lamda(i) -1.8*lamda(i)^2;
            H(i) = 2.61 - 3.75*lamda(i) + 5.24*lamda(i)^2;
            S(i) = (lamda(i)+0.09)^0.62;
        else if lamda(i)<=0 && lamda(i)>-0.1
            L(i) = 0.22 + 1.402*lamda(i) + 0.018*lamda(i)/(lamda(i)+0.107);
            H(i) = 2.088 + 0.0731/(lamda(i)+0.14);
            S(i) = (lamda(i)+0.09)^0.62;
        else if 0.1<= lamda <= 0.25
            H(i)= 2.0 + 4.14*z - 83.5*z^2 + 854*z^3 - 3337*z^4 + 4576*z^5 ;
            S(i) = (lamda(i)+0.09)^0.62;
            L(i) = L(i-1);
        else
            L(i) = L(i-1);
            H(i) = H(i-1);
            S(i) = S(i-1);
            end
            end
        end
        if teta(i)==0
            Cf1(i) = 0;
        else
            Cf1(i) = 2*L(i)*nu/(Uin(i)*teta(i));
        end
        tauw(i) = S(i)*miu*Uin(i)/teta(i);
        delta(i) = teta(i)*H(i);
        Cf2(i) = tauw(i)/(0.5*rho*Uin(i)^2);
    end
    
    %lower
    z = 0.25 - lamda(i);
    for i = sl:-1:1
        if lamda(i)<0.1 && lamda(i)>0
            L(i) = 0.22 + 1.57*lamda(i) -1.8*lamda(i)^2;
            H(i) = 2.61 - 3.75*lamda(i) + 5.24*lamda(i)^2;
            S(i) = (lamda(i)+0.09)^0.62;
        else if lamda(i)<=0 && lamda(i)>-0.1
            L(i) = 0.22 + 1.402*lamda(i) + 0.018*lamda(i)/(lamda(i)+0.107);
            H(i) = 2.088 + 0.0731/(lamda(i)+0.14);
            S(i) = (lamda(i)+0.09)^0.62;
        else if 0.1<= lamda <= 0.25
            H(i)= 2.0 + 4.14*z - 83.5*z^2 + 854*z^3 - 3337*z^4 + 4576*z^5 ;
            S(i) = (lamda(i)+0.09)^0.62;
            L(i) = L(i+1);
        else
            L(i) = L(i+1);
            H(i) = H(i+1);
            S(i) = S(i+1);
            end
            end
        end
        if teta(i)==0
            Cf1(i) = 0;
        else
            Cf1(i) = 2*L(i)*nu/(Uin(i)*teta(i));
        end
        tauw(i) = S(i)*miu*Uin(i)/teta(i);
        delta(i) = teta(i)*H(i);
        Cf2(i) = tauw(i)/(0.5*rho*Uin(i)^2);
    end
    %evaluate separation point
    %upper
    for i=su:M
        if lamda(i)<=-0.09
            xsepu = X(i);
            k1 = i;
            break
        else
            xsepu = X(M);
            k1 = M;
        end
    end
    %lower
    for i=sl:-1:1
        if lamda(i)<=-0.09
            xsepl = X(i);
            k2 = i;
            break
        else
            xsepl = X(1);
            k2 = 1;
        end
    end

    %boundary layer thickness
    for i=1:M
       if Y(i)>0
           YBL(i) = Y(i) + delta(i);
       else
           YBL(i) = Y(i) - delta(i);
       end
    end
    
    Gee = (Uin .* delta)';
    Gee(MP1) = Gee(1);
     %% ..............Calculate Transition Location...............
     
     %Upper Section
     for i=su:M
     Reteta = Uin(i)* teta(i)/nu;
     Rex = Uin(i)* X(i)/nu;
     RR = 2.8*Rex^0.4;
        if Reteta >= RR
            transp1 = i;
            break;
        else
            transp1 = M-1;
        end
     end
     %Lower Section
     for i=sl:-1:1
     Reteta = Uin(i)* teta(i)/nu;
     Rex = Uin(i)* X(i)/nu;
     RR = 2.8*Rex^0.4;
        if Reteta >= RR
            transp2 = i;
            break;
        else
            transp2 = 1;
        end
     end