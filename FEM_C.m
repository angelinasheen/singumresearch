%% Clear workspace
clear; clc; close all;

%% 1. Create PDE model (general scalar PDE)
model = createpde();

%% 2. Define geometry
% Outer rectangle (2m x 2m)
R1 = [3,4,0,2,2,0,0,0,2,2]';

% Inner rectangle (center region with different conductivity)
inner_length = 0.01;
inner_width  = 0.02;
center_x = 1;
center_y = 1;
x_left   = center_x - inner_length/2;
x_right  = center_x + inner_length/2;
y_bottom = center_y - inner_width/2;
y_top    = center_y + inner_width/2;
R2 = [3,4, x_left, x_right, x_right, x_left, y_bottom, y_bottom, y_top, y_top]';

gd = [R1,R2];
ns = char('R1','R2')';
sf = 'R1+R2';

[dl,bt] = decsg(gd,sf,ns);
geometryFromEdges(model, dl);

%% 3. Plot geometry with face & edge labels 
figure
pdegplot(model,'FaceLabels','on','EdgeLabels','on')
axis equal
title('Geometry with Outer Rectangle and Inner Rectangle (edge labels shown)')

%% 4. Assign material properties
% Face 1: outer rectangle, isotropic
Kxx1 = 6e4; Kxy1 = 0; Kyx1 = 0; Kyy1 = 6e4;
c_tensor1 = [Kxx1; Kxy1; Kyx1; Kyy1];
specifyCoefficients(model, 'm', 0, 'd', 0, 'c', c_tensor1, 'a', 0, 'f', 0, 'Face', 1);

% Face 2: inner rectangle, orthotropic
Kxx2 = 2e4; Kxy2 = 0; Kyx2 = 0; Kyy2 = 4e4;
c_tensor2 = [Kxx2; Kxy2; Kyx2; Kyy2];   
specifyCoefficients(model, 'm', 0, 'd', 0, 'c', c_tensor2, 'a', 0, 'f', 0, 'Face', 2);

% 5. Set voltage boundary conditions on rectangle edges
applyBoundaryCondition(model,'dirichlet','Edge',2,'u',5);   % left
applyBoundaryCondition(model,'dirichlet','Edge',4,'u',-5);  % right

%% 6. Generate mesh
generateMesh(model,'Hmax',0.02);

%% 6a. Visualize FEM mesh
figure('Name','FEM Mesh','NumberTitle','off')
pdemesh(model.Mesh)
axis equal
xlabel('x (m)')
ylabel('y (m)')
title('Finite Element Mesh')

%% 7. Solve PDE
result = solvepde(model);
phi = result.NodalSolution;    % (numNodes x 1)

%% 8. Plot electric potential distribution
figure
pdeplot(model,'XYData',phi,'Contour','on')
xlabel('x (m)')
ylabel('y (m)')
axis equal
title('Electric Potential \phi (V)')
colorbar

%% 9. Compute current density J = -sigma * grad(phi) at nodes
nodes = model.Mesh.Nodes;      % 2 x numNodes
numNodes = size(nodes,2);

% gradients at nodes (numNodes x 1)
[gradx, grady] = evaluateGradient(result);  % already column vectors

% construct nodal sigma：外域 6e4，内域 4e4
sigma_out = 6e4;
sigma_in  = 4e4;
x_nodes = nodes(1,:)';
y_nodes = nodes(2,:)';
sigma_nodes = sigma_out * ones(numNodes,1);

in_inner_mask = (x_nodes >= x_left - 1e-12) & (x_nodes <= x_right + 1e-12) & ...
                (y_nodes >= y_bottom - 1e-12) & (y_nodes <= y_top + 1e-12);
sigma_nodes(in_inner_mask) = sigma_in;

Jx = - sigma_nodes .* gradx(:);
Jy = - sigma_nodes .* grady(:);

%% 10. Plot Jx and Jy
figure
pdeplot(model,'XYData',Jx)
title('Current Density J_x (nodal)')
xlabel('x (m)'); ylabel('y (m)'); axis equal; colorbar;

figure
pdeplot(model,'XYData',Jy)
title('Current Density J_y (nodal)')
xlabel('x (m)'); ylabel('y (m)'); axis equal; colorbar;

%% 11. Potential along center/quarter lines
F_phi = scatteredInterpolant(x_nodes, y_nodes, phi(:), 'natural', 'none');
x_line = linspace(0,2,200)';    % column vector
y_mid = 1;
y_quarter = 0.5;

phi_line = F_phi(x_line, y_mid*ones(size(x_line)));
phi_quarter = F_phi(x_line, y_quarter*ones(size(x_line)));

figure
plot(x_line, phi_line,'b-','LineWidth',2)
xlabel('x (m)'); ylabel('\phi (V)'); title('Electric Potential Along y = 1'); grid on

figure
plot(x_line, phi_quarter,'r-','LineWidth',2)
xlabel('x (m)'); ylabel('\phi (V)'); title('Electric Potential Along y = 0.5'); grid on

%% 12. Jx along those lines
F_Jx = scatteredInterpolant(x_nodes, y_nodes, Jx(:), 'natural', 'none');
Jx_center = F_Jx(x_line, y_mid*ones(size(x_line)));
Jx_quarter = F_Jx(x_line, y_quarter*ones(size(x_line)));

figure
plot(x_line, Jx_center,'b-','LineWidth',2,'DisplayName','y=1 (center)')
hold on
plot(x_line, Jx_quarter,'r-','LineWidth',2,'DisplayName','y=0.5 (quarter)')
xlabel('x (m)'); ylabel('J_x (A/m^2)'); title('J_x Along Lines'); legend('Location','best'); grid on

%% 13. Vector field (downsampled arrows)
step = 10;
idx = 1:step:numNodes;
x_sub = x_nodes(idx); y_sub = y_nodes(idx);
Jx_sub = Jx(idx); Jy_sub = Jy(idx);
mag = sqrt(Jx_sub.^2 + Jy_sub.^2);

figure('Position',[100,100,1200,900]); hold on;
colormap jet; cmap = colormap; nC = size(cmap,1);
if max(mag)-min(mag) < eps
    C_idx = ones(size(mag));
else
    C_idx = ceil((mag - min(mag))/(max(mag)-min(mag))*(nC-1))+1;
end
scaleFactor = 3e-7;
for k = 1:length(x_sub)
    quiver(x_sub(k), y_sub(k), Jx_sub(k)*scaleFactor, Jy_sub(k)*scaleFactor, 0, ...
        'Color', cmap(C_idx(k),:), 'LineWidth', 1, 'MaxHeadSize',0.8);
end
axis equal; xlabel('x (m)'); ylabel('y (m)');
title('Current density vector field'); colorbar; c = colorbar; c.Label.String = 'J magnitude (A/m^2)';

%% 14. Electric Current Density Streamslice with magnitude-colored background

% 1. Interpolant functions for Jx, Jy, and magnitude
F_Jx_grid  = scatteredInterpolant(x_nodes, y_nodes, Jx, 'natural', 'none');
F_Jy_grid  = scatteredInterpolant(x_nodes, y_nodes, Jy, 'natural', 'none');
F_mag_grid = scatteredInterpolant(x_nodes, y_nodes, sqrt(Jx.^2 + Jy.^2), 'natural', 'none');

% 2. Generate uniform grid
numGrid = 60;
[Xgrid, Ygrid] = meshgrid(linspace(min(x_nodes), max(x_nodes), numGrid), ...
                          linspace(min(y_nodes), max(y_nodes), numGrid));

% 3. Evaluate Jx, Jy and magnitude on grid
U = F_Jx_grid(Xgrid, Ygrid);  
V = F_Jy_grid(Xgrid, Ygrid);   
Mag = F_mag_grid(Xgrid, Ygrid);

% 4. Plot background (magnitude)
figure('Position', [250,250,1000,800]); hold on;
pcolor(Xgrid, Ygrid, Mag);
shading interp;
colormap jet;
caxis([min(Mag(:)), max(Mag(:))]);
cb = colorbar;
cb.Label.String = 'Current Density Magnitude (A/m^2)';

% 5. Overlay streamlines
streams = streamslice(Xgrid, Ygrid, U, V, 2); 
set(streams, 'Color','k', 'LineWidth', 1.2)  

axis equal
axis tight
xlabel('x (m)')
ylabel('y (m)')
title('Electric Current Density Streamslice (Colored Background)')
grid on


% add label on curves
%text(xq_col(end), Jx_center(end), ' y = 0.866', 'VerticalAlignment','bottom', 'HorizontalAlignment','left')
%text(xq_col(end), Jx_quarter(end), ' y = 0.433', 'VerticalAlignment','bottom', 'HorizontalAlignment','left')
%% 14. Save CSV (Jx along center & quarter)
data_center = [x_line, y_mid*ones(size(x_line)), Jx_center];
writematrix(data_center,'current_density_Jx_center.csv');
data_quarter = [x_line, y_quarter*ones(size(x_line)), Jx_quarter];
writematrix(data_quarter,'current_density_Jx_quarter.csv');

disp('Done.');
