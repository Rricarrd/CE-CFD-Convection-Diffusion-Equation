

// Dimensions
int M = 64;              // number of columns
int N = 2 * M;           // number of rows
double H = 1;            // height of the channel
double L = 2 * H;        // length of the channel
const double dx = L / N; // x step
const double dy = L / M; // y step

// Numerical
const int time_steps = 10000;             // number of time steps
const double delta_convergence = 0.00001; // maximum delta for the error
const double initial_phi = 0;             // initial value
const double delta_t = 0.001;             // time step

// Physical
const double rho = 1;   // inlet density
const double S_phi = 0; // source term depenant on phi
const double S_c = 0;   // constant source term

// Diagonal flow
const double u0 = 50;   // inlet velocity for diagonal flow
const double phi_h = 2; // phi value for high side
const double phi_l = 0; // phi value for low side
