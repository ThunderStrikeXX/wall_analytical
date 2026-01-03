/**
 * Analytical solution of the 1D transient heat conduction equation
 * in a solid wall with internal heat generation.
 *
 * Domain:
 *   z ∈ [0, L]
 *
 * Governing equation:
 *   ρ c_p ∂T/∂t = k ∂²T/∂z² + Q
 *
 * Boundary conditions:
 *   - z = 0 : zero heat flux (adiabatic)
 *             ∂T/∂z = 0
 *
 *   - z = L : prescribed temperature
 *             T = 300 K
 *
 * Initial condition:
 *   - Uniform temperature field
 *             T(z, 0) = 300 K
 *
 * Source term:
 *   - Uniform, constant volumetric heat generation
 *             Q = const [W/m³]
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <omp.h>

int main() {

    // Physics and domain
    constexpr int N = 100;                  // Number of wall cells
    constexpr double L = 1.0;               // Wall length [m]
	constexpr double dz = L / (N - 1);      // Wall cell size [m]
    constexpr double dt = 1e-3;             // Time step [s]
    constexpr int time_iter = 1000;         // Number of time iterations
    constexpr int harm = 100;               // Number of harmonics [-]
    const double pi = acos(-1.0);

    constexpr double k = 20.0;              // Steel thermal conductivjy [W/mK]
    constexpr double rho = 7850.0;          // Steel density [kg/m3]
    constexpr double cp = 500.0;            // Steel specific heat [J/kgK]
    constexpr double T_amb = 300.0;         // Ambient temperature [K]
    constexpr double Q = 1e8;               // Heat pipe volumetric source term [W/m3]

	// Temperature vector
    std::vector<double> T(N, 300.0);
    
    // Output file
    std::ofstream file("wall_analytical.dat");

    // Coefficient A_n from the initial conditions
    auto A_n = [&](int n) {
        double sign = std::pow(-1.0, n);
        double val = sign * -16 * Q * L * L / (k * std::pow(pi, 3) * std::pow(2 * n + 1, 3));
        return val;
    };

	std::vector<double> lambda_vec(harm);
	std::vector<double> An_vec(harm);
	std::vector<double> xx_vec(N);

    double start = omp_get_wtime();

    for (int i = 0; i < N; ++i) xx_vec[i] = i * dz;

    for (int n = 0; n < harm; ++n) {
        
        lambda_vec[n] = (2.0 * n + 1.0) * pi / (2.0 * L);
        An_vec[n] = A_n(n);
    }

    // Time loop
    for (int j = 0; j < time_iter; ++j) {

        double t = dt * j;

        // Node loop
		// Note: parallelizazion here is not useful and slows down the execution
        for (int i = 0; i < N; ++i) {

            double Ts = T_amb + Q / (2.0 * k) * (L * L - xx_vec[i] * xx_vec[i]);

            // Transient solution
            double Tt = 0.0;
            for (int n = 0; n < harm; ++n) {

                double cosine = std::cos(lambda_vec[n] * xx_vec[i]);
                double expo = std::exp(-k / (rho * cp) * lambda_vec[n] * lambda_vec[n] * t);
                Tt += An_vec[n] * cosine * expo;
            }

			// Superposition of stationary and transient solution
            T[i] = Ts + Tt;
        }

        // Output
        for (int i = 0; i < N; ++i)
            file << T[i] << " ";

        file << "\n";
    }

    double end = omp_get_wtime();
    std::cout << "Execution time: " << end - start;

    file.flush();
    file.close();

    return 0;
}
