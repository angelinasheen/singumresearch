clear all;
% 1. Model type: DirectCurrent
model = femodel(AnalysisType="dcConduction");

% 2. Define Geometry
R1 = [3,4,0,2,2,0,0,0,1.732,1.732]'; % R1：rectangle
C2 = [1,0.8,0.866,0.05,0,0,0,0,0,0]';% C2：left electrode/center (0.8, 0.866), radius 0.05
C3 = [1,1.2,0.866,0.05,0,0,0,0,0,0]';% C3：right electrode/center (1.2, 0.866), radius 0.05

gd = [R1, C2, C3];                   % Combine geometry
ns = char('R1','C2','C3')';
sf = 'R1 + C2 + C3'; 
[dl,bt] = decsg(gd,sf,ns);

% 3. Generate geometry object and assign it to the model
g = geometryFromEdges(dl);
model.Geometry = g;

% 4. Plot geometry and show edge labels to help set boundary conditions
figure
pdegplot(model.Geometry, 'EdgeLabels', 'on')
axis equal

% 5. Material properties：electrical conductivity 6e4 S/m
model.MaterialProperties = materialProperties(ElectricalConductivity=6e4);

% 6. Set voltage boundary conditions
model.EdgeBC(5) = edgeBC(Voltage=-5);     % left electrode -5 V
model.EdgeBC(6) = edgeBC(Voltage=-5);     % left electrode -5 V
model.EdgeBC(7) = edgeBC(Voltage=-5);     % left electrode -5 V
model.EdgeBC(8) = edgeBC(Voltage=-5);     % left electrode -5 V

model.EdgeBC(9) = edgeBC(Voltage=5);      % right electrode +5 V
model.EdgeBC(10) = edgeBC(Voltage=5);     % right electrode +5 V
model.EdgeBC(11) = edgeBC(Voltage=5);     % right electrode +5 V
model.EdgeBC(12) = edgeBC(Voltage=5);     % right electrode +5 V


% 7. Generate mesh with maximum element size 0.02
model = generateMesh(model, 'Hmax', 0.02);

% 7a. Visualize FEM mesh
figure
pdemesh(model.Mesh)
title('Finite Element Mesh')
axis equal


% 8. Solve model
R = solve(model);

% 9. Visualize electric potential distribution
figure
pdeplot(model.Mesh, 'XYData', R.ElectricPotential, 'Contour', 'on')
xlabel('x(m)')
ylabel('y(m)')
axis equal
cb = colorbar;                     
cb.Label.String = 'V';  

% 10. Visualize current density Jx
figure
pdeplot(model.Mesh, 'XYData', R.CurrentDensity.Jx)
title('Current Density J_x')
xlabel('x (m)')
ylabel('A/m^2')
axis equal
cb = colorbar;
cb.Label.String = 'A/m^2'; 

% 11. Visualize current density Jy
figure
pdeplot(model.Mesh, 'XYData', R.CurrentDensity.Jy)
title('Current Density J_y')
xlabel('x (m)')
ylabel('A/m^2')
axis equal
cb = colorbar;
cb.Label.String = 'A/m^2';  

% 12. Visualize Electric potential in 3D
figure
pdeplot(model.Mesh, 'XYData', R.ElectricPotential, 'ZData', R.ElectricPotential)
title('Electric Potential \phi (Surface Plot)')
xlabel('x(m)')
ylabel('y(m)')
zlabel('\phi (V)')
axis tight
view(45, 30)
cb = colorbar;
cb.Label.String = 'V'; 

% 13. Plot electric potential along the mid horizontal line 
nodes = model.Mesh.Nodes; 
p = nodes;
u = R.ElectricPotential;
F = scatteredInterpolant(p(1,:)', p(2,:)', u);

xq = linspace(0, 2, 200);       % generate coordinate(X)
yq = ones(size(xq)) * 0.866;    % generate coordinate(Y)
Vq = F(xq, yq);
figure                          % plot
plot(xq, Vq, 'LineWidth', 2)
xlabel('x (m)')
ylabel('Electric Potential \phi (V)')
%title('Electric Potential Along y = 0.866')
grid on

% 14. obtain Electric potential along quarter line & plot
yq2 = ones(size(xq)) * 0.433;  % generate coordinate(Y)
Vq2 = F(xq, yq2);          
figure                         % plot
plot(xq, Vq2, 'LineWidth', 2)
xlabel('x (m)')
ylabel('Electric Potential \phi (V)')
%title('Electric Potential Along y = 0.433')
grid on

% 15. obtain current density Jx/Jy along center line and quarter line, plot
Jx_vals = R.CurrentDensity.Jx;
F_Jx = scatteredInterpolant(p(1,:)', p(2,:)', Jx_vals);
xq_col = xq(:);
yq_col = yq(:);
yq2_col = yq2(:);
Jx_center = F_Jx(xq_col, yq_col);     % Jx(center line)  from Interpolation
Jx_quarter = F_Jx(xq_col, yq2_col);   % Jy(quarter line) from Interpolation

figure                                % plot
plot(xq_col, Jx_center, 'LineWidth', 2, 'DisplayName', 'c_1 - c_2 FEM')
hold on
plot(xq_col, Jx_quarter, 'LineWidth', 2, 'DisplayName', 'c_3 - c_4 FEM')
xlabel('x (m)')
ylabel('Current Density J_x (A/m^2)')
%title('Current Density J_x Along Lines')
legend('Location','best')
grid on

% add label on curves
%text(xq_col(end), Jx_center(end), ' y = 0.866', 'VerticalAlignment','bottom', 'HorizontalAlignment','left')
%text(xq_col(end), Jx_quarter(end), ' y = 0.433', 'VerticalAlignment','bottom', 'HorizontalAlignment','left')

% save as CSV
data_Jx_center = [xq_col, yq_col, Jx_center];
writematrix(data_Jx_center, 'current_density_Jx_center_fig2.csv');
data_Jx_quarter = [xq_col, yq2_col, Jx_quarter];
writematrix(data_Jx_quarter, 'current_density_Jx_quarter_fig3.csv');
