clear
close

%% ELEC4700 Assignment 2 Part 1(b)

VoL = 100; % Left side of the Area has Boundary Voltage = VoL
VoR = 100; % Right side of the Area has Boundary Voltage = VoR
VoT = 0; % Top side of the Area has Boundary Voltage = VoT
VoB = 0; % Bottom side of the Area has Boundary Voltage = VoB

L = 5; % Area Length
W = 4; % Area Width
GridSize = 0.1; % Each mesh grid has the Length and Width of 0.01
nx = L/GridSize; % Length mesh size
ny = W/GridSize; % Width mesh size
index = 150; % index # of analytical series
delta = 1; % Delta x and Delta y of FD condition 
G = sparse(nx*ny); % G matrix has size(nx*ny,nx*ny)
B = zeros(nx*ny,1); % B is the product of G matrix * V
VTheory = zeros(nx,ny); % Voltage Analytical series solution

%% （b) Solve the case where V = V0 at x = 0, x = L and V = 0 at y = 0, y = W. Compare the solution of a bunch of mesh sizes to the analytical series solution


for iRow = 1: nx
    for jColumn = 1:ny
        n = jColumn+(iRow-1)*ny;
        % Left side Boundary Condition
        if iRow == 1     
            G(n,:) = 0;
            G(n,n) = 1/delta;
            B(n) = VoL;
        % Right side Boundary Condition
        elseif iRow == nx
            G(n,:) = 0;
            G(n,n) = 1/delta;
            B(n) = VoR;

        % Bottom side Boundary Condition
        elseif jColumn == 1   
            G(n,n) = 1/delta;
            B(n) = VoB;
        % Top side Boundary Condition
        elseif jColumn == ny
            G(n,n) = 1/delta;
            B(n) = VoT;     
        else 
            nxm = jColumn+((iRow-1)-1)*ny;
            nxp = jColumn+((iRow+1)-1)*ny;
            nym = (jColumn-1)+(iRow-1)*ny;
            nyp = (jColumn+1)+(iRow-1)*ny;
            G(n,n) = -4/delta;
            G(n,nxm) = 1/delta;
            G(n,nxp) = 1/delta;
            G(n,nym) = 1/delta;
            G(n,nyp) = 1/delta;
        end

    end

end
fullG = full(G);

% Find ny*ny:1 size V = G\B
Vn = G\B; 
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
% Creates a 3D Voltage numerical surface plot

figure('name','3D Plot of Numerical V(x,y)')
set(surf(X',Y',V),'linestyle','none');
title("3D Plot of Numerical V(x,y)")
xlabel("Length")
ylabel("Width")
zlabel("Voltage (V)")


% 
% Gfull = full(G);
% 
% figure(1)
% spy(G)

% Analytical solution
a = nx; % set a
b = ny; % set b
x1 = linspace(-b,b,nx); % set x of analytical solution
y1 = linspace(0,a,ny); % set y of analytical solution
for k=1:2:index % odd index only
    for i = 1:nx
        for j = 1:ny
            % Analytical series solution
            VTheory(i,j) = VTheory(i,j) + 4*VoL/pi*1/k*cosh(k*pi*x1(i)/a)/cosh(k*pi*b/a)*sin(k*pi*y1(j)/a);                           
        end 
    end
end

% Creates a 3D Voltage Analytical surface plot
figure('name','3D Plot of Analytical V(x,y)')
set(surf(X',Y',VTheory),'linestyle','none');
title("3D Plot of Analytical V(x,y)")
xlabel("Length")
ylabel("Width")
zlabel("Voltage (V)")

% Determine the numerical and analytical solution's differences
error = V - VTheory;
figure('name','3D Plot of V(x,y) Differences')
set(surf(X',Y',error),'linestyle','none');
title("3D Plot of V(x,y) Differences")
xlabel("Length")
ylabel("Width")
zlabel("Voltage (V)")


%% 
% Conclusions on meshing and comments on numerical vs analytical
    % As shown in the plot of the error between the analytical and
    % numerical Voltage solutions, the difference between two solutions 
    % are very small at the middle of the area. Howeve, at each corner of
    % the area, the differences are significant. The reason causes this happen
    % is because the corners are relatively hard to converge for the
    % analytical solutions. It usually takes more terms (close to infinite)
    % to obtain the solution with the least errors compare to the numerical one.

    % Numerical and analytical method are both great ways to modeling FD.
    % The numerical solution have the advantage like: easy to calculate,
    % and when there is a condition that there are no analytical 
    % equations, the numerical will still be useful.
    % It's disadvantage is: It usually needs high mesh density to get an
    % accurate solution.

    % For Analytical method, it is usually more accurate than numerical
    % method, becasue Analytical method are presented as math expressions.
% 