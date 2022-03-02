clear
close

%% ELEC4700 Assignment 2 Part 2(a)

VoL = 100; % Left side of the Area has Boundary Voltage = VoL
VoR = 0; % Right side of the Area has Boundary Voltage = VoR
VoT = 0; % Top side of the Area has Boundary Voltage = VoT
VoB = 0; % Bottom side of the Area has Boundary Voltage = VoB
Pixel = 10; % Number of mesh per unit length/width

Lb = 2; % Box Length
Wb = 2; % Box Width
xBox = Lb*Pixel; % Box length Pixel Number
yBox = Wb*Pixel; % Box width Pixel Number
L = 9; % Area Length
W = 6; % Area Width
nx = L*Pixel; % Area length Pixel Number
ny = W*Pixel; % Area wdith Pixel Number
sigma = 1; % Outside Box Area Conductivity
BoxSigma = 0.01; % Inside Box Area Conductivity
G = sparse(nx*ny); % G matrix has size(nx*ny,nx*ny)
B = zeros(nx*ny,1); % B is the product of G matrix * V
Conductivity = sigma*ones(nx,ny); % Conductivity of the entire area
%% (a) Calculate the current flow at the two contacts. Generate plots of Cond(x,y), V(x,y), E(x,y), J(x,y)

% At the left side of the area contact, follow euqation R =
% rho*Length/Area. In this case is a 2D model, so the area will be just
% Width. R = sigma*Length/Width = 1*9/6 = 1.5. Also, Current = Voltage/R =
% VoL/1.5 = 1/1.5 = 0.6667A.

% Set Inside Box area conductivity to variable BoxSigma
for iRow = 1:nx
    for jColumn = 1:ny

        if iRow>=(nx-xBox)/2 && iRow<=((nx+xBox)/2) && jColumn<=yBox
            Conductivity(iRow,jColumn) = BoxSigma;
        elseif iRow>=(nx-xBox)/2 && iRow<=((nx+xBox)/2) && jColumn<=ny && jColumn>ny-yBox
            Conductivity(iRow,jColumn) = BoxSigma;
        end
    end
end

% for jColumn = 1:ny
%     for iRow = 1:nx
for iRow = 1:nx
    for jColumn = 1:ny
        n = jColumn+(iRow-1)*ny;
        % Left side Boundary Condition
       
        if iRow == 1      
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = VoL;
        % Right side Boundary Condition
           
        elseif iRow == nx
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = VoR;

        % Bottom side Boundary Condition
        elseif jColumn == 1    
            nxm = jColumn+((iRow-1)-1)*ny;
            nxp = jColumn+((iRow+1)-1)*ny;
            nyp = (jColumn+1)+(iRow-1)*ny;

            rxm = (Conductivity(iRow,jColumn) + Conductivity(iRow-1,jColumn))/2;
            rxp = (Conductivity(iRow,jColumn) + Conductivity(iRow+1,jColumn))/2;
            ryp = (Conductivity(iRow,jColumn) + Conductivity(iRow,jColumn+1))/2;
            
            G(n,n) = -(rxm+rxp+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nyp) = ryp;
            
        % Top side Boundary Condition
        elseif jColumn == ny 
            nxm = jColumn+((iRow-1)-1)*ny;
            nxp = jColumn+((iRow+1)-1)*ny;
            nym = (jColumn-1)+(iRow-1)*ny;

            rxm = (Conductivity(iRow,jColumn) + Conductivity(iRow-1,jColumn))/2;
            rxp = (Conductivity(iRow,jColumn) + Conductivity(iRow+1,jColumn))/2;
            rym = (Conductivity(iRow,jColumn) + Conductivity(iRow,jColumn-1))/2;
            
            G(n,n) = -(rxm+rxp+rym);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;

        else 
            nxm = jColumn+((iRow-1)-1)*ny;
            nxp = jColumn+((iRow+1)-1)*ny;
            nym = (jColumn-1)+(iRow-1)*ny;
            nyp = (jColumn+1)+(iRow-1)*ny;

            rxm = (Conductivity(iRow,jColumn) + Conductivity(iRow-1,jColumn))/2;
            rxp = (Conductivity(iRow,jColumn) + Conductivity(iRow+1,jColumn))/2;
            rym = (Conductivity(iRow,jColumn) + Conductivity(iRow,jColumn-1))/2;
            ryp = (Conductivity(iRow,jColumn) + Conductivity(iRow,jColumn+1))/2;
            
            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
        end

    end

end


% 3D Plot of Voltage V(x,y)
figure('name','Voltage V(x,y)')

Vn = G\B; % Find (ny*ny:1) size, V = G\B
% Mapping the V to the matrix size of nx*ny
for iRow = 1: nx
    for jColumn = 1:ny
         n = jColumn+(iRow-1)*ny;
         V(iRow,jColumn) = Vn(n);    
    end
end

x = linspace(0,L,nx);
y = linspace(0,W,ny);
[X,Y] = meshgrid(x,y);

surf(x,y,V')
view(31,23)
xlim([0 L])
ylim([0 W])
title("3D Plot of Voltage V(x,y)")
xlabel("Length")
ylabel("Width")
zlabel("Voltage (V)")

% 3D Plot of Electric Field E(x,y)
[Ex,Ey] = gradient(V');
figure('name','Electric Field E(x,y)')
quiver(X,Y,-Ex,-Ey)
xlim([0 L])
ylim([0 W])
title("3D Plot of Electric Field E(x,y)")
xlabel("Length")
ylabel("Width")
zlabel("Electric Field (V/m)")

figure('name','Electric Field E(x)')
surf(X,Y,-Ex)
xlim([0 L])
ylim([0 W])
title("3D Plot of Electric Field E(x)")
xlabel("Length")
ylabel("Width")
zlabel("Electric Field (V/m)")
view(2)

figure('name','Electric Field E(y)')
surf(X,Y,-Ey)
xlim([0 L])
ylim([0 W])
title("3D Plot of Electric Field E(y)")
xlabel("Length")
ylabel("Width")
zlabel("Electric Field (V/m)")
view(2)

%E = sqrt((Ex.^2)+ (Ey.^2));

% 3D Plot of Conductivity Cond(x,y)
figure('name','Conductivity Cond(x,y)')
cond = Conductivity';
surf(X,Y,cond)
xlim([0 L])
ylim([0 W])
title("3D Plot of Conductivity Cond(x,y)")
xlabel("Length")
ylabel("Width")
zlabel("Conductivity (ohm*m)")

% 3D Plot of Current Density J(x,y)
figure('name','Current density J(x,y)')
Jx = cond.*(-Ex);
Jy = cond.*(-Ey);


quiver(X,Y,Jx,Jy)
xlim([0 L])
ylim([0 W])
title("3D Plot of Current Density J(x,y)")
xlabel("Length")
ylabel("Width")
zlabel("Current Density (A/Area)")