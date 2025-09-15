%% Parameters 
L   = 2.000;               
W   = 1.732;               
r_e = 0.05;                % electrode radius
Vn  = -5;                   % left disk potential
Vp  = 5;                  % right disk potential
%R   = 1.443e-5;            % link resistance (Ω)
K   = 60000;                 % bond conductance (S)

% lattice spacing
lp0 = 0.005;       
dx  = 2*lp0;       
dy  = sqrt(3)*lp0; 

nx = round(L/dx) + 1;
ny = round(W/dy) + 1;

%% Generate nodes & node_map (unchanged)
nodes    = zeros(nx*ny,3);   % [id, x, y]
node_map = zeros(ny,nx);
id = 1;
for j = 1:ny
  for i = 1:nx
    x = (i-1)*dx + mod(j,2)*lp0;
    y = (j-1)*dy;
    nodes(id,:)   = [id, x, y];
    node_map(j,i) = id;
    id = id+1;
  end
end
coords = nodes(:,2:3);
N      = size(nodes,1);

%% Build conductance matrix (same) + collect edges for diagnostics
M = sparse(N,N);
edges = [];                 % [n m] once per undirected edge
for j = 1:ny
  for i = 1:nx
    n = node_map(j,i);
    if i < nx
      m = node_map(j,i+1);
      M([n m],[n m]) = M([n m],[n m]) + [ K -K; -K  K ];
      edges(end+1,:) = [n m]; %#ok<AGROW>
    end
    if j < ny
      m1 = node_map(j+1,i);
      M([n m1],[n m1]) = M([n m1],[n m1]) + [ K -K; -K  K ];
      edges(end+1,:) = [n m1]; %#ok<AGROW>
      if mod(j,2)==1 && i<nx
        m2 = node_map(j+1,i+1);
        M([n m2],[n m2]) = M([n m2],[n m2]) + [ K -K; -K  K ];
        edges(end+1,:) = [n m2]; %#ok<AGROW>
      elseif mod(j,2)==0 && i>1
        m2 = node_map(j+1,i-1);
        M([n m2],[n m2]) = M([n m2],[n m2]) + [ K -K; -K  K ];
        edges(end+1,:) = [n m2]; %#ok<AGROW>
      end
    end
  end
end

%% Dirichlet BCs on two circular electrodes (unchanged)
c1 = [ (L-0.4)/2, W/2 ];   % [0.8, 0.866]
c2 = [ (L+0.4)/2, W/2 ];   % [1.2, 0.866]
d1 = vecnorm(coords - c1,2,2);
d2 = vecnorm(coords - c2,2,2);
left_nodes  = find(d1 <= r_e);
right_nodes = find(d2 <= r_e);

bc_nodes = [left_nodes; right_nodes];
bc_vals  = [Vn*ones(numel(left_nodes),1);
            Vp*ones(numel(right_nodes),1)];

b = zeros(N,1);
for k = 1:numel(bc_nodes)
  n = bc_nodes(k);
  M(n,:) = 0;  M(n,n) = 1;   % φ(n) = bc_vals(k)
  b(n)   = bc_vals(k);
end

%% Solve potentials (unchanged)
V = M\b;

%% Interpolate φ to a regular grid
ng = 300;
xg = linspace(0, L, ng);
yg = linspace(0, W, ng);
[X, Y] = meshgrid(xg, yg);
Z = griddata(coords(:,1), coords(:,2), V, X, Y, 'cubic');

%% -------- Current density J = -sigma * grad(phi) --------
% Use physical sheet conductivity to match FEM magnitudes:
sigma_target = 6e4;                  % S/m

dxg = xg(2)-xg(1);  dyg = yg(2)-yg(1);
[dVdx, dVdy] = gradient(Z, dxg, dyg);  % first = ∂φ/∂x, second = ∂φ/∂y

Jx = -sigma_target * dVdx;            % A/m
Jy = -sigma_target * dVdy;            % A/m


%% Heatmap of Jx
figure('Color','w');
pcolor(X, Y, Jx); shading interp; axis equal tight; box on
colormap('cool'); cb=colorbar; cb.Label.String='A/m';
xlabel('x (m)'); ylabel('y (m)');
title('Current density J_x (A/m)');

%% -------- Jx along midline (y=W/2) and bottom quarter line (y=W/4) ------
y_mid = W/2;   y_quarter = W/4;
JxI = griddedInterpolant({yg,xg}, Jx, 'linear','none');  % Jx is (y,x)

x_line = xg;
Jx_mid     = JxI(y_mid*ones(size(x_line)),     x_line);
Jx_quarter = JxI(y_quarter*ones(size(x_line)), x_line);

% % mask points inside disks on each line
% mask_mid = ((x_line - c1(1)).^2 + (y_mid     - c1(2)).^2 <= r_e^2) | ...
%            ((x_line - c2(1)).^2 + (y_mid     - c2(2)).^2 <= r_e^2);
% mask_qua = ((x_line - c1(1)).^2 + (y_quarter - c1(2)).^2 <= r_e^2) | ...
%            ((x_line - c2(1)).^2 + (y_quarter - c2(2)).^2 <= r_e^2);
% Jx_mid(mask_mid)       = NaN;
% Jx_quarter(mask_qua)   = NaN;

figure('Color','w');
plot(x_line, Jx_mid,     'LineWidth',2, 'DisplayName','midline  y = W/2'); hold on
plot(x_line, Jx_quarter, 'LineWidth',2, 'DisplayName','quarter y = W/4');
grid on; box on; axis tight
xlabel('x (m)'); ylabel('J_x (A/m)');
title('Current density J_x along midline and bottom quarter line');
legend('Location','best');

% Save CSVs (x, y, Jx)
writematrix([x_line(:) repmat(y_mid,    numel(x_line),1) Jx_mid(:)],     'Jx_midline_SigNum.csv');
writematrix([x_line(:) repmat(y_quarter,numel(x_line),1) Jx_quarter(:)], 'Jx_quarterline_SigNum.csv');
