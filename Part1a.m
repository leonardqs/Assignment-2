clear
close

%% ELEC4700 Assignment 2 Part 1(a)

VoL = 1; % Left side of the Area has Boundary Voltage = VoL
VoR = 0; % Right side of the Area has Boundary Voltage = VoR
L = 3; % Area Length
W = 2; % Area Width
GridSize = 0.1; % Each mesh grid has the Length and Width of 0.01
nx = L/GridSize; % Length mesh size
ny = W/GridSize; % Width mesh size
delta = 1; % Delta x and Delta y of FD condition 
G = sparse(nx*ny,nx*ny); % G matrix has size(nx*ny,nx*ny)
B = zeros(nx*ny,1); % B is the product of G matrix * V


%% （a) Solve the simple case where V = V0 at x = 0 and V = 0 at x = L.Note that in this case the top/bottom BC are not fixed

% 
for Horizontal = 1: nx
    for Vertical = 1:ny
        n = Vertical+(Horizontal-1)*ny;
        % Left side Boundary Condition
        if Horizontal == 1 
            G(n,:) = 0;
            G(n,n) = 1/delta;
            B(n) = VoL;
        % Right side Boundary Condition
        elseif Horizontal == nx
            G(n,:) = 0;
            G(n,n) = 1/delta;
            B(n) = VoR;
        
        else 
            nxm = Vertical+((Horizontal-1)-1)*ny;
            nxp = Vertical+((Horizontal+1)-1)*ny;
            G(n,n) = -2/delta;
            G(n,nxm) = 1/delta;
            G(n,nxp) = 1/delta;
        end
    end
end

Vn = G\B; % Find (ny*ny:1), size V = G\B
% Mapping the V to the matrix size of nx*ny
for irow = 1: nx 
    for jcolumn = 1:ny
         n = jcolumn+(irow-1)*ny;
         V(irow,jcolumn) = Vn(n);    
    end
end
GMatrix = full(G);
figure('name','2D Plot of V(x)')

% Creates a 2D Voltage surface plot
x = linspace(0,L,nx);
y = linspace(0,W,ny);
H = surf(x,y,V');

set(H,'linestyle','none');
title("2D Plot of V(x)")
xlabel("Length")
ylabel("Width")
zlabel("Voltage (V)")
xlim([0 L])
ylim([0 W])
view(2)
