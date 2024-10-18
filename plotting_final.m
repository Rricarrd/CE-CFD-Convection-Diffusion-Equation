clear;
% Some parameters
M = 32;   
N = 2 * M;     
H = 1;     
L = 2 * H;  
dy = H / M; 
dx = L / N; 
t = 900;
len = (M-2)*(N-2);

n1= t*len + 1;
n2= (t+1)*len;

% Data input and preprocessing
data = readtable('output_t1000.csv');
t = table2array(data(n1:n2,1));
X = table2array(data(n1:n2,2));
Y = table2array(data(n1:n2,3));
U = table2array(data(n1:n2,4));
V = table2array(data(n1:n2,5));
phi = table2array(data(n1:n2,6));
Vtot = sqrt(V.^2 + U.^2);

%% %% PHI CONTOUR PLOT
figure(1)
%Contour plotting
X = reshape(flip(X),[M-2,N-2])
Y = reshape(flip(Y),[M-2,N-2])
phi = reshape(flip(phi),[M-2,N-2])
contourf(X,Y,phi)

% [x_grid,y_grid] = meshgrid(linspace(dx,L-dx,N-2),linspace(dy,H-dy,M-2));
% phi_grid = griddata(X, Y, phi ,x_grid,y_grid); %interpolates surface from  mesh and streamline values (cubic interpolation)
% contourf(x_grid,y_grid,phi_grid);

% Colorbar
c_bar = colorbar;
c_bar.Label.String = 'Phi';

%Plot parameters
xlabel('X-axis [m]');
ylabel('Y-axis [m]');
title('Phi variable','Interpreter','latex');
grid on
xlim([dx, L-dx])
ylim([dy, H-dy])
colormap cool

% %% %% VELOCITY CONTOUR PLOT
% figure(2)
% %Contour plotting
% [x_grid,y_grid] = meshgrid(linspace(dx,L-dx,M-2),linspace(dy,H-dy,N-2)); 
% v_grid = griddata(X, Y, Vtot ,x_grid,y_grid); %interpolates surface from  mesh and streamline values (cubic interpolation)
% contourf(x_grid,y_grid,v_grid);
% 
% % Colorbar
% c_bar = colorbar;
% c_bar.Label.String = 'Velocity [m/s]';
% 
% %Plot parameters
% xlabel('X-axis [m]');
% ylabel('Y-axis [m]');
% title('Velocity field','Interpreter','latex');
% xlim([dx, L-dx]);
% ylim([dy, H-dy]);
% grid on
% colormap cool

