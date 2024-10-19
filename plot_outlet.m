% Parameters
M = 64;
type="smith-hutton";
location = "schemes";
n_fig = 1;


name="COMPARISON";
fileref = 'reference.csv';
dataref = table2array(readtable(fileref));


for Pe = [10,1000,1000000]
    x_ref = dataref(:,1);
    phi_ref = dataref(:,n_fig+1);
    plot_phi_outlets(Pe,M,type,location,name, n_fig, x_ref, phi_ref);
    n_fig = n_fig + 1;

end

% Comparacions fetes amb delta t = 0.01
% Deltatlow amb deltat = 0.001


