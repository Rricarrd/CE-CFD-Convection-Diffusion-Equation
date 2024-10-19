

// Dimensions
int M = 128;             // number of columns
int N = 2 * M;           // number of rows
double H = 1;            // height of the channel
double L = 2 * H;        // length of the channel
const double dx = L / N; // x step
const double dy = L / M; // y step

// Numerical
const int time_steps = 1500;              // number of time steps
const double delta_convergence = 0.00001; // maximum delta for the error
const double initial_phi = 0;             // initial value
const double delta_t = 0.001;             // time step
const double relaxation_factor = 1;       // relaxation factor

// Physical
const double rho = 1;   // inlet density
const double S_phi = 0; // source term depenant on phi
const double S_c = 0;   // constant source term

// Diagonal flow
const double u0 = 1;     // inlet velocity for diagonal flow
const double phi_h = 1;  // phi value for high side
const double phi_l = -1; // phi value for low side
