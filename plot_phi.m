function plot_phi(Pe,M,ts,type,scheme,location, n_fig)

file = sprintf('%s/output_Pe_%i_S_%s_M_%i_type_%s.csv',location,Pe,scheme,M,type)
data = readtable(file);

% Preprocesing

H = 1;

if type == "diagonal"
    N = M;     
    L = H;
else
    N = 2 * M;     
    L = 2 * H;
    
end

add_x = L/2;
add_y = H/2;

dy = H / M; 
dx = L / N; 

len = (M-2)*(N-2);

% Reading values
try
    if ts == "last"
        X = table2array(data(:,2));
        Y = table2array(data(:,3));
        phi = table2array(data(:,4));
    end
catch
    n1= ts*len + 1;
    n2= (ts+1)*len;

    X = table2array(data(n1:n2,2));
    Y = table2array(data(n1:n2,3));
    phi = table2array(data(n1:n2,4));
end


% Reshaping vectors
X = reshape(flip(X),[M-2,N-2])
Y = reshape(flip(Y),[M-2,N-2])
phi = reshape(flip(phi),[M-2,N-2])


% PHI CONTOUR PLOT
figure(n_fig)

%Contour plotting
pcol = pcolor(X+add_x,Y+add_y,phi);
set(pcol,'EdgeColor', 'none');


% Colorbar
c_bar = colorbar;
c_bar.Label.String = 'Phi';

%Plot parameters
xlabel('X-axis [m]');
ylabel('Y-axis [m]');
title(sprintf('Phi for Pe = %i', Pe),'Interpreter','latex');
grid off
colormap hot
end
