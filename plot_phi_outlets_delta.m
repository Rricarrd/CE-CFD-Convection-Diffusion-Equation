function plot_phi_outlets_delta(Pe,delta,M,type,location,name, n_fig, x_ref, phi_ref)

file = sprintf('%s/SCHEMES_Pe_%i_S_%s_M_%i_type_%s_t_%.1g.csv',location,Pe,name,M,type,delta);
data = readtable(file);

% Preprocesing
H = 1;
N = 2 * M;     

tab = rmmissing(data(:,2:end));

phi = table2array(tab);


% Getting data
X = linspace(0,1,size(phi,2));



% PHI CONTOUR PLOT
figure(n_fig)


plot(X, phi(1,:))
hold on
plot(X, phi(2,:))
plot(X, phi(3,:))

if delta < 0.005
    plot(X, phi(4,:))
end
plot(x_ref, phi_ref,"-o") 

%Plot parameters
xlabel('X-axis [m]');
ylabel('Phi$[\phi]$','Interpreter','latex');
title(sprintf('Phi for Pe = %i and  dt =%.1g ', Pe, delta),'Interpreter','latex');

if delta < 0.005
    legend("UDS", "CDS", "HDS", "PDS", "Reference")
else
    legend("UDS", "HDS", "PDS", "Reference")
end

saveas(figure(n_fig),sprintf('plots/outlet_Pe_%i_S_%s_M_%i_type_%s_t_%d.png',Pe,name,M,type,delta));


end
