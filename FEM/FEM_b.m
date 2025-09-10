clear all;
% 1. Model type: DirectCurrent
model = femodel(AnalysisType="dcConduction");

% 2. Define Geometry
R1 = [3,4, 0, 2, 2, 0, 0, 0, 2, 2]';     
hole_length = 0.01;                       
hole_width = 0.01;                        
center_x = 1;
center_y = 1;
x_left   = center_x - hole_length/2;
x_right  = center_x + hole_length/2;
y_bottom = center_y - hole_width/2;
y_top    = center_y + hole_width/2;
R2 = [3,4, x_left, x_right, x_right, x_left, y_bottom, y_bottom, y_top, y_top]';   

gd = [R1,R2];                           
ns = char('R1','R2')';
sf = 'R1-R2';
[dl,bt] = decsg(gd,sf,ns);

% 3. Generate geometry object and assign it to the model
g = geometryFromEdges(dl);
model.Geometry = g;

% 4. Plot geometry and show edge labels to help set boundary conditions
figure
pdegplot(model.Geometry, EdgeLabels="on")
title("2D Geometry with Rectangular Hole")
axis equal

% 5. Set material property: electrical conductivity 6e4 S/m
model.MaterialProperties = materialProperties(ElectricalConductivity=6e4);

% 6. Set voltage BCs on the OUTER GEOMETRY left/right sides (no mesh tricks)

% Get geometry bounds (option A: via vertices)
V = model.Geometry.Vertices;          % [x y] for all geometry vertices
xmin = min(V(:,1));  xmax = max(V(:,1));
ymin = min(V(:,2));  ymax = max(V(:,2));
figure
pdegplot(model.Geometry,'EdgeLabels','on')
axis equal
title('Check edge numbers for outer boundary')
model.EdgeBC(4) = edgeBC(Voltage=-5);   % left boundary
model.EdgeBC(2) = edgeBC(Voltage=+5);   % right boundary


% 7. Generate mesh with maximum element size 0.02
model = generateMesh(model,'Hmax',0.02);

% 7a. Visualize FEM mesh
figure
pdemesh(model.Mesh)
title('Finite Element Mesh')
axis equal

% 8. Solve model
R = solve(model);

% 9. Plot electric potential distribution
figure
pdeplot(model.Mesh, 'XYData', R.ElectricPotential, 'Contour', 'on')
hold on
y_mid = 0.5;
xlabel('x(m)')
ylabel('y(m)')
axis equal
cb = colorbar;                     
cb.Label.String = 'V'; 

% 10. visualize current density Jx
figure
pdeplot(model.Mesh,'XYData',R.CurrentDensity.Jx)
title('Current Density J_x')
xlabel('x (m)')
ylabel('A/m^2')
axis equal
cb = colorbar;
cb.Label.String = 'A/m^2';

% 11. Plot y-component of current density
figure
pdeplot(model.Mesh,'XYData',R.CurrentDensity.Jy)
title('Current Density J_y')
axis equal
cb = colorbar;
cb.Label.String = 'A/m^2';

% 12. Visualize Electric potential in 3D
figure
pdeplot(model.Mesh, 'XYData', R.ElectricPotential, 'ZData', R.ElectricPotential)
title('Electric Potential \phi (Surface Plot)')
xlabel('x')
ylabel('y')
zlabel('\phi (V)')
axis tight
view(45, 30)
colorbar
cb = colorbar;
cb.Label.String = 'V'; 

% 13. Plot electric potential along the mid horizontal line 
p = model.Mesh.Nodes;
phi = R.ElectricPotential(:); 
F_interp = scatteredInterpolant(p(1,:)', p(2,:)', phi);
x_line = linspace(0, 1, 200);
y_line = y_mid * ones(size(x_line));
phi_line = F_interp(x_line', y_line');

figure
plot(x_line, phi_line, 'b-', 'LineWidth', 2);
xlabel('x coordinate (m)');
ylabel('Electric Potential \phi (V)');
title('Electric Potential Along Mid Horizontal Line');
grid on

% 14. obtain Electric potential along quarter line & plot
nodes = model.Mesh.Nodes; 
p = nodes;
u = R.ElectricPotential;
F = scatteredInterpolant(p(1,:)', p(2,:)', u);
xq = linspace(0, 2, 200);           
yq = ones(size(xq)) * 1;            
Vq = F(xq, yq);

figure                              
plot(xq, Vq, 'LineWidth', 2)
xlabel('x (m)')
ylabel('Electric Potential \phi (V)')
grid on

% 14. obtain Electric potential along quarter line & plot
yq2 = ones(size(xq)) * 0.5;         
Vq2 = F(xq, yq2);

figure                             
plot(xq, Vq2, 'LineWidth', 2)
xlabel('x (m)')
ylabel('Electric Potential \phi (V)')
grid on

% % 15. obtain current density Jx/Jy along center line and quarter line, plot
Jx_vals = R.CurrentDensity.Jx;
F_Jx = scatteredInterpolant(p(1,:)', p(2,:)', Jx_vals);
xq_col = xq(:);
yq_col = yq(:);
yq2_col = yq2(:);
Jx_center = F_Jx(xq_col, yq_col);     
Jx_quarter = F_Jx(xq_col, yq2_col);   

% ----------hole area set to NA------------------
x_hole_left = center_x - hole_length/2;
x_hole_right = center_x + hole_length/2;
Jx_center_masked = Jx_center;
Jx_center_masked(xq_col >= x_hole_left & xq_col <= x_hole_right) = NaN;

% if quarter line cross the hole
Jx_quarter_masked = Jx_quarter;
%Jx_quarter_masked(xq_col >= x_hole_left & xq_col <= x_hole_right) = NaN;
% -----------------------------------------------

figure                                  % plot
plot(xq_col, Jx_center_masked, 'LineWidth', 2, 'DisplayName', 'y = 1 (center line)')
hold on
plot(xq_col, Jx_quarter_masked, 'LineWidth', 2, 'DisplayName', 'y = 0.5 (quarter line)')
xlabel('x (m)')
ylabel('Current Density J_x (A/m^2)')
title('Current Density J_x Along Lines (with Hole Masked)')
legend('Location','best')
grid on

% save as CSV
data_Jx_center = [xq_col, yq_col, Jx_center_masked];
writematrix(data_Jx_center, 'current_density_Jx_center_fig3.csv');

data_Jx_quarter = [xq_col, yq2_col, Jx_quarter_masked];
writematrix(data_Jx_quarter, 'current_density_Jx_quarter_fig3.csv');
