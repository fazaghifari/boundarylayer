function [qwall,delc,Ythermal] = thermal(Twall,Tfs,U,Vtan,X,Y)
    miu = 0.0000181206; %Dynamics Viscosity
    rho = 1.225;        %Density
    nu = miu/rho;       %kinematics Viscosity
    k = 0.024;          %Thermal Conductivity
    alpha = 1.9e-5;     %Thermal Diffusivity
    
    Pr = nu/alpha;      %Prandtl Number
    a = (1/(0.332*Pr^0.35))^2; 
    b = 2.95*Pr^0.07;
    M = length(X);
    MP1 = M+1;
    Uin = Vtan*U;
    
    Ythermal = zeros(M,1);
    qwall    = zeros(M,1);
    

    %%%............Calculate Conduction Thickness............%%%
    
    delc = zeros(M,1);
    % mencari titik stagnasi
    sl = find(Uin==min(Uin));
    su = sl;
    
    %upper airfoil
    delc(su) = sqrt(3.96*nu/(abs((Uin((su)+1)-Uin(su))/(X((su)+1)-X((su))))));
    for j=su+1:M
        integral1 = 0;
        for i=su+1:j
            integral1 = integral1 + (Uin(i)^(b-1) + Uin(i-1)^(b-1))*abs((X(i)-X(i-1)))/2 ;
        end
        delc(j) = sqrt(a*nu*(integral1)/(Uin(j)^b));
    end    
    
    %lower airfoil
    delc(sl) = sqrt(3.96*nu/(abs((Uin((sl)-1)-Uin(sl))/(X((sl)-1)-X((sl))))));
    for j=sl-1:-1:1
        integral1 = 0;
        for i=sl-1:-1:j
            integral1 = integral1 + (Uin(i)^(b-1) + Uin(i+1)^(b-1))*abs((X(i)-X(i+1)))/2 ;
        end
        delc(j) = sqrt(a*nu*(integral1)/(Uin(j)^b));
    end
    
    for i=1:M
       if Y(i)>0
           Ythermal(i) = Y(i) + delc(i);
       else
           Ythermal(i) = Y(i) - delc(i);
       end
       qwall(i) = k*(Twall-Tfs)/delc(i);
    end
end