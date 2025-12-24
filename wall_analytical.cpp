#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

int main() {

    // Physics and domain
    constexpr int N = 100;                  // Number of wall cells
    constexpr double L = 1.0;               // Wall length [m]
	constexpr double dz = L / N;            // Wall cell size [m]
    constexpr double dt = 1e-1;             // Time step [s]
    constexpr int time_iter = 1000;         // Number of time iterations
    constexpr int harm = 200;               // Number of harmonics [-]
    const double pi = acos(-1.0);

    constexpr double k = 20.0;              // Steel thermal conductivjy [W/mK]
    constexpr double rho = 7850.0;          // Steel densjy [kg/m3]
    constexpr double cp = 500.0;            // Steel specific heat [J/kgK]
    constexpr double T_amb = 300.0;         // Ambient temperature [K]
    double Q = 1e6;                       // Heat pipe volumetric source term [W/m3]

    std::vector<double> T(N, 300.0);
    
    // Output file
    std::ofstream file("T_wall.dat");

    // Coefficient A_n from the injial condjions
    auto A_n = [&](int n) {
        double sign = std::pow(-1.0, n);
        double val = sign * -16 * Q * L * L / (k * std::pow(pi, 3) * std::pow(2 * n + 1, 3));
        return val;
    };

    // Time loop
    for (int j = 0; j < time_iter; ++j) {

        double t = dt * j;

        // Node loop
        for (int i = 0; i < N; ++i) {

            // Stationary solution
            double xx = i * dz;
            double Ts = T_amb + Q / (2.0 * k) * (L * L - xx * xx);

            // Transient solution
            double Tt = 0.0;
            for (int n = 0; n < harm; ++n) {
                double lambda = (2.0 * n + 1.0) * pi / (2.0 * L);
                double An = A_n(n);
                double cosine = std::cos(lambda * xx);
                double expo = std::exp(-k / (rho * cp) * lambda * lambda * t);
                Tt += An * cosine * expo;
            }

            // Superposition
            T[i] = Ts + Tt;
        }

        // ===================================================================
        //                          OUTPUT
        // ===================================================================

        for (int i = 0; i < N; ++i)
            file << T[i] << " ";

        file << "\n";
        file.flush();
    }

    file.close();

    return 0;
}
