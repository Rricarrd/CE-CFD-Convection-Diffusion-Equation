function plot_outlets(Pe,M,type,location, n_fig)

file = sprintf('%s/SCHEMES_Pe_%i_S_COMPARISON_M_%i_type_%s.csv',location,Pe,M,type);
data = readtable(file);

% Preprocesing

H = 1;
N = 2 * M;     


phi = table2array(data(:,4));


% Reshaping vectors
X = linspace(0,1,M);

phi = reshape(flip(phi),[M-2,N-2]);

% Getting outlet nodes



% PHI CONTOUR PLOT
figure(n_fig)

%Contour plotting
plot(X(0,N/2:end), phi(0,N/2:end))

%Plot parameters
xlabel('X-axis [m]');
ylabel('Y-axis [m]');
title(sprintf('Phi for Pe = %i', Pe),'Interpreter','latex');
grid off
colormap hot
end
