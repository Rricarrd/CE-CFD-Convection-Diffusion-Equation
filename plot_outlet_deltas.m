% Parameters
M = 32;
type="smith-hutton";
location = "deltat";
n_fig = 1;


name="COMPARISON_delta";

fileref = 'reference.csv';
dataref = table2array(readtable(fileref));
Pe = 1000;

for delta = [0.1, 0.01, 0.001, 0.0001]
    x_ref = dataref(:,1);
    phi_ref = dataref(:,3);
    plot_phi_outlets_delta(Pe,delta,M,type,location,name, n_fig, x_ref, phi_ref);
    n_fig = n_fig + 1;

end

% Comparacions fetes amb delta t = 0.01
% Deltatlow amb deltat = 0.001


