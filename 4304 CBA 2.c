#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Simulation Parameters & Constants
#define C_VEL 4.0        // Wave propagation speed (m/s)
#define X_MIN 0.0      // Left boundary of the domain
#define X_MAX 12.0       // Right boundary of the domain
#define U_LEFT 2.0       // Left boundary condition (m/s)
#define U_RIGHT 1.0      // Right boundary condition (m/s)
#define T_FINAL 1.5      // Final simulation time (seconds)

// Enumeration for selecting the numerical scheme
typedef enum { UPWIND, LAX_WENDROFF, CRANK_NICOLSON } SchemeType;

// Function Prototypes
void initialize_grid(double *x, double *u, int nx, double dx);
void exact_solution(double *u_exact, double *x, int nx, double t);
void solve_upwind(double *u, double *u_new, int nx, double C_courant);
void solve_lax_wendroff(double *u, double *u_new, int nx, double C_courant);
void solve_crank_nicolson(double *u, double *u_new, int nx, double C_courant, double *a, double *b, double *c_diag, double *d);
void tdma(double *a, double *b, double *c, double *d, double *x, int n);
void write_output(const char *filename, double *x, double *u, double *u_exact, int nx);

int main() {
    // Grid and Time Step Configuration

    double dx = 0.025;          // Spatial step size (x) 
    double dt = 0.005;          // Time step size (s)  
    int nx = (int)((X_MAX - X_MIN) / dx);  // Number of spatial nodes
    int nt = (int)(T_FINAL / dt);          // Number of spatial nodes
    double Courant = (C_VEL*dt)/dx;
    printf("========================================================\n");
    printf("|  1D Linear Convection Solver - AERO4304 Assignment 2 | \n");
    printf("========================================================\n");
    printf("Grid nodes (nx): %d\n", nx);
    printf("Grid nodes (nt): %d\n", nt);
    printf("Spatial step (dx): %.4f m\n", dx);
    printf("Time step (dt): %.5f s\n", dt);
    printf("Courant Number (C): %.4f\n", Courant);

    // Allocate dynamic arrays
    double *x = (double *)malloc(nx * sizeof(double));
    double *u = (double *)malloc(nx * sizeof(double));
    double *u_new = (double *)malloc(nx * sizeof(double));
    double *u_exact = (double *)malloc(nx * sizeof(double));
    
    // Additional arrays for the Crank-Nicolson TDMA solver
    double *a = (double *)malloc(nx * sizeof(double));
    double *b = (double *)malloc(nx * sizeof(double));
    double *c_diag = (double *)malloc(nx * sizeof(double));
    double *d = (double *)malloc(nx * sizeof(double));

    // Loop through all three schemes
    for (SchemeType scheme = UPWIND; scheme <= CRANK_NICOLSON; scheme++) {
        initialize_grid(x, u, nx, dx);
        
        const char *scheme_name;
        const char *file_name;
        
        switch (scheme) {
            case UPWIND:
                scheme_name = "Euler Explicit First-Order Upwind";
                file_name = "output_euler_upwind.csv";
                break;
            case LAX_WENDROFF:
                scheme_name = "Explicit Lax-Wendroff";
                file_name = "output_lax_wendroff.csv";
                break;
            case CRANK_NICOLSON:
                scheme_name = "Implicit Crank-Nicolson";
                file_name = "output_crank_nicolson.csv";
                break;
        }
        
        printf("Running scheme: %s...\n", scheme_name);
        
        // Time integration loop
        double current_time = 0.0;
        for (int n = 0; n < nt; n++) {
            
            // Apply numerical schemes to compute inner nodes
            if (scheme == UPWIND) {
                solve_upwind(u, u_new, nx, Courant);
            } else if (scheme == LAX_WENDROFF) {
                solve_lax_wendroff(u, u_new, nx, Courant);
            } else if (scheme == CRANK_NICOLSON) {
                solve_crank_nicolson(u, u_new, nx, Courant, a, b, c_diag, d);
            }
            
            // Apply Dirichlet Boundary Conditions at the edges
            u_new[0] = U_LEFT;
            u_new[nx - 1] = U_RIGHT;
            
            // Update the solution for the next time step
            for (int i = 0; i < nx; i++) {
                u[i] = u_new[i];
            }
            
            current_time += dt;
        }
        
        // Compute analytical solution at final time for comparison
        exact_solution(u_exact, x, nx, current_time);
        
        // Output results to file
        write_output(file_name, x, u, u_exact, nx);
        printf("Results saved to %s\n\n", file_name);
    }

    // Free all allocated memory
    free(x); free(u); free(u_new); free(u_exact);
    free(a); free(b); free(c_diag); free(d);

    printf("Simulations complete.\n");
    return 0;
}

//Function Implementations

// Initializes grid coordinates and applies initial conditions (Step function)
void initialize_grid(double *x, double *u, int nx, double dx) {
    for (int i = 0; i < nx; i++) {
        x[i] = X_MIN + i * dx;
        if (x[i] <= 0.0) {
            u[i] = U_LEFT;
        } else {
            u[i] = U_RIGHT;
        }
    }
}

// Computes the exact analytical solution u(x,t) = F(x - ct)
void exact_solution(double *u_exact, double *x, int nx, double t) {
    for (int i = 0; i < nx; i++) {
        double x_initial = x[i] - C_VEL * t;
        if (x_initial <= 0.0) {
            u_exact[i] = U_LEFT;
        } else {
            u_exact[i] = U_RIGHT;
        }
    }
}

// 1. Euler Explicit Discretization & First-Order Upwind Differencing
void solve_upwind(double *u, double *u_new, int nx, double C_courant) {
    for (int i = 1; i < nx - 1; i++) {
        u_new[i] = u[i] - C_courant * (u[i] - u[i - 1]);
    }
}

// 2. Explicit Lax-Wendroff Scheme
void solve_lax_wendroff(double *u, double *u_new, int nx, double C_courant) {
    for (int i = 1; i < nx - 1; i++) {
        u_new[i] = u[i] - 0.5 * C_courant * (u[i + 1] - u[i - 1]) 
                   + 0.5 * pow(C_courant, 2) * (u[i + 1] - 2.0 * u[i] + u[i - 1]);
    }
}

// 3. Implicit Crank-Nicolson Scheme setup
void solve_crank_nicolson(double *u, double *u_new, int nx, double C_courant, double *a, double *b, double *c_diag, double *d) {
    int n_inner = nx - 2; // Number of unknown internal nodes
    double alpha = C_courant / 4.0;

    // Build the Tridiagonal System for inner nodes i = 1 to nx-2
    for (int i = 0; i < n_inner; i++) {
        int grid_idx = i + 1; // Map TDMA index to actual grid index
        
        a[i] = -alpha;       // Sub-diagonal
        b[i] = 1.0;          // Main diagonal
        c_diag[i] = alpha;   // Super-diagonal
        
        // Right Hand Side (RHS) vector based on known values at time 'n'
        d[i] = u[grid_idx] - alpha * (u[grid_idx + 1] - u[grid_idx - 1]);
    }

    // Apply Boundary Conditions into the RHS vector (d)
    d[0] += alpha * U_LEFT;               // Left boundary affects first inner equation
    d[n_inner - 1] -= alpha * U_RIGHT;    // Right boundary affects last inner equation

    // Solve the tridiagonal system using TDMA
    double *temp_sol = (double *)malloc(n_inner * sizeof(double));
    tdma(a, b, c_diag, d, temp_sol, n_inner);

    // Map the temporary solution back to the main u_new array
    for (int i = 0; i < n_inner; i++) {
        u_new[i + 1] = temp_sol[i];
    }
    
    free(temp_sol);
}

// Thomas Algorithm (TDMA) for solving Tridiagonal systems
// a = sub-diagonal, b = main-diagonal, c = super-diagonal, d = RHS, x = solution
void tdma(double *a, double *b, double *c, double *d, double *x, int n) {
    double *c_prime = (double *)malloc(n * sizeof(double));
    double *d_prime = (double *)malloc(n * sizeof(double));

    // Forward sweep
    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];

    for (int i = 1; i < n; i++) {
        double m = 1.0 / (b[i] - a[i] * c_prime[i - 1]);
        c_prime[i] = c[i] * m;
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) * m;
    }

    // Back substitution
    x[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    free(c_prime);
    free(d_prime);
}

// Writes spatial distribution data to a tab-separated text file
void write_output(const char *filename, double *x, double *u, double *u_exact, int nx) {
    FILE *fp;
    char path[256];
    sprintf(path, "C:\\Users\\J-for\\OneDrive\\Documents\\C - 4304\\%s", filename);

    // create directory first if needed
    // on Windows, you can do it manually before running

    fp = fopen(path,"w");
    if (fp == NULL) {
        printf("Error opening file %s for writing.\n", filename);
        return;
    }
    
    fprintf(fp, "X_Coordinate,Numerical_U,Exact_U\n");
    for (int i = 0; i < nx; i++) {
        fprintf(fp, "%.6f,%.6f,%.6f\n", x[i], u[i], u_exact[i]);
    }
    
    fclose(fp);
}