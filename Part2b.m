clear
%% ELEC4700 Assignment 2 Part 2(b)

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

meshDensityIndex = [1:1:5];
current = zeros(length(meshDensityIndex),1);
current1 = zeros(length(meshDensityIndex),1);
current2 = zeros(length(meshDensityIndex),1);
%currentDensity = zeros(meshDensityIndex,1);

%% (b) Graph of current vs mesh size. Investigate mesh density.

for meshSize = 1:1:length(meshDensityIndex)
    
%     nx=L*meshSize*Pixel;
%     ny=W*meshSize*Pixel;
%     xBox = Lb*meshSize*Pixel;
%     yBox = Wb*meshSize*Pixel;
    nx = L*Pixel *meshDensityIndex(meshSize);
    ny = W*Pixel*meshDensityIndex(meshSize);
   
    xBox = Lb*Pixel*meshDensityIndex(meshSize);
    yBox = Wb*Pixel*meshDensityIndex(meshSize);
    
    Conductivity = sigma*ones(nx,ny);
    %From Before
    
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

    [Ex,Ey] = gradient(V');
    cond = Conductivity';
    Jx = cond.*(-Ex);
    J = Conductivity'.*gradient(-(V'));
    
%     a = sum(abs(Jx(:,45))).*W;
%     a1 = sum(abs(Jx(:,1))).*W;
%     a2 = sum(abs(Jx(:,41))).*W;

    current(meshSize,1) = max(J,[],'all');
%     current1(meshSize,1) = a1;
%     current2(meshSize,1) = a2;
    %CurrentDensity = max(J,[],"all");
    
end
%surf(V);
figure(1)
plot(meshDensityIndex,current)
xlabel('mesh density');
ylabel('Max Current Density (A/m^2)');
    
    
