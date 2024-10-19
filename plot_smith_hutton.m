% Parameters
M = 64;
t = "last"
type="smith-hutton";
scheme="PDS";
location = "smith_hutton"
n_fig = 1;


for Pe = [10,1000,1000000]
    plot_phi(Pe,M,t,type,scheme,location, n_fig);
    n_fig = n_fig + 1;
end