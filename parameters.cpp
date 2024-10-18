// Dimensions
const int M = 32;        // number of columns
const int N = 2 * M;     // number of rows
const double H = 1;      // height of the channel
const double L = 2 * H;  // length of the channel
const double dx = L / N; // x step
const double dy = L / M; // y step

// Numerical
const int time_steps = 1000;                // number of time steps
const double delta_convergence = 0.0000001; // maximum delta for the error
const double initial_value = 0.5;           // initial value
const double delta_t = 0.001;               // time step
const double relaxation_factor = 1;         // relaxation factor

// Physical
const double rho = 1;   // inlet density
const double S_phi = 0; // source term depenant on phi
const double S_c = 0;   // constant source term
