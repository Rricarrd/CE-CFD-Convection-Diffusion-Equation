// Dimensions
const int M = 16;        // number of columns
const int N = 2 * M;     // number of rows
const double H = 1;      // height of the channel
const double L = 2;      // length of the channel
const double dx = L / N; // x step
const double dy = L / M; // ystep

// Numerical
const int time_steps = 1000;                // number of time steps
const double delta_convergence = 0.0000001; // maximum delta for the error
const double initial_value = 0;             // initial value
const double delta_t = 0.01;                // time step
const double relaxation_factor = 1;         // relaxation factor

// Physical

const double p_in = 100000;  // inlet pressure
const double t_in = 298;     // inlet temperature
const double v_in = 1;       // inlet velocity
const double rho_in = 1.225; // inlet density
