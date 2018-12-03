function [Cf,delta] = boundarylayer(U,Vtan,X,Y)
    %% ________________________________THWAITES__________________________________
    %Referensi : Dinamika Fluida (Lavi R Zuhal), Viscous Flow (Frank White),Intro
                % to theoretical and computational aerodynamics (Moran)

    Uin = Vtan*U;
    %% .......................Calculate momentum thickness........................
    teta = zeros(M,1);
    % mencari titik stagnasi
    for i=1:M
        if Uin(i) == min(Uin)
            sl = i;
        end
    end
    su = sl+1;

    %upper airfoil
    teta(su) = sqrt(0.075*nu/(abs((Uin((su)+1)-Uin(su))/(X((su)+1)-X((su))))));
    for j=su+1:M
        integral1 = 0;
        for i=su+1:j
            integral1 = integral1 + (Uin(i)^5 + Uin(i-1)^5)*abs((X(i)-X(i-1)))/2 ;
        end
        teta(j) = sqrt(0.45*nu*integral1/(Uin(j)^6));
    end
    %lower airfoil
    teta(sl) = sqrt(0.075*nu/(abs((Uin((sl)-1)-Uin(sl))/(X((sl)-1)-X((sl))))));
    for j=sl-1:-1:1
        integral1 = 0;
        for i=sl-1:-1:j
            integral1 = integral1 + (Uin(i)^5 + Uin(i+1)^5)*abs((X(i)-X(i+1)))/2 ;
        end
        teta(j) = sqrt(0.45*nu*integral1/(Uin(j)^6));
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
    H = zeros(M,1);
    Cf = zeros(M,1);
    delta = zeros(M,1);
    YBL = zeros(M,1);
    %shear stress at wall
    %upper
    for i = su:M
        if lamda(i)<0.1 && lamda(i)>0
            L(i) = 0.22 + 1.57*lamda(i) -1.8*lamda(i)^2;
            H(i) = 2.61 - 3.75*lamda(i) + 5.24*lamda(i)^2;
        else if lamda(i)<=0 && lamda(i)>-0.1
            L(i) = 0.22 + 1.402*lamda(i) + 0.018*lamda(i)/(lamda(i)+0.107);
            H(i) = 2.088 + 0.0731/(lamda(i)+0.14);
            else
            L(i) = L(i-1);
            H(i) = H(i-1);
            end
        end
        if teta(i)==0
            Cf(i) = 0;
        else
            Cf(i) = 2*L(i)*nu/(Uin(i)*teta(i));
        end
        delta(i) = teta(i)*H(i);
    end
    %lower
    for i = sl:-1:1
        if lamda(i)<0.1 && lamda(i)>0
            L(i) = 0.22 + 1.57*lamda(i) -1.8*lamda(i)^2;
            H(i) = 2.61 - 3.75*lamda(i) + 5.24*lamda(i)^2;
        else if lamda(i)<=0 && lamda(i)>-0.1
            L(i) = 0.22 + 1.402*lamda(i) + 0.018*lamda(i)/(lamda(i)+0.107);
            H(i) = 2.088 + 0.0731/(lamda(i)+0.14);
            else
            L(i) = L(i+1);
            H(i) = H(i+1);
            end
        end
        if teta(i)==0
            Cf(i) = 0;
        else
            Cf(i) = 2*L(i)*nu/(Uin(i)*teta(i));
        end
        delta(i) = teta(i)*H(i);
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