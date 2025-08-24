#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

const int Euler = 0;
const int VelocityVerlet = 1;
const int Beeman = 2;
const int RK4 = 3;

std::vector<std::vector<std::string>> filename = {{"Euler0.txt", "VelocityVerlet0.txt", "Beeman0.txt", "RK40.txt"},
                                                  {"Euler1.txt", "VelocityVerlet1.txt", "Beeman1.txt", "RK41.txt"}};

double xsol(double t, double F0, double k, double m, double nu = 0.0) {
  if (nu == 0.0) {
    double omega = sqrt(k / m);
    return (F0 / k) * (1.0 - cos(omega * t));
  } else if (nu * nu > 4.0 * m * k) {
    double sqrDelta = sqrt(nu * nu - 4.0 * m * k);
    double r1 = (-nu + sqrDelta) / (2.0 * m);
    double r2 = (-nu - sqrDelta) / (2.0 * m);
    double C1 = (-F0 * sqrDelta) / (k * (nu + sqrDelta));
    double C2 = -F0 / k - C1;
    return C1 * exp(r1 * t) + C2 * exp(r2 * t) + (F0 / k);
  }
  // here is the case for nu <= 2*sqrt(m*k) (nu != 0)
  double sqrDelta = sqrt(4.0 * m * k - nu * nu);
  double alpha = -nu / (2.0 * m);
  double beta = (sqrDelta) / (2.0 * m);
  double Ccos = -F0 / k;
  double Csin = (-alpha / beta) * Ccos;
  return (Csin * sin(beta * t) + Ccos * cos(beta * t)) * exp(alpha * t) + (F0 / k);
}

double acc(double k, double m, double F0, double x, double v, double nu = 0.0) { return (-k * x - nu * v + F0) / m; }

double calculateTotalEnergy(double k, double m, double v, double x, double F0) {
  double potentialEnergy = 0.5 * k * x * x;
  double kineticEnergy = 0.5 * m * v * v;
  return potentialEnergy + kineticEnergy;
}

int main(int argc, char const* argv[]) {

  double m = 1.0;
  double k = 1.0;
  double nu0 = 0.5 * 2.0 * sqrt(m * k);
  double nu = 0.0;
  double F0 = 1.0;

  double t = 0.0;
  double x = 0.0;
  double v = 0.0;
  double vsave;
  double a = 0.0;
  double aprev, anext;

  double sqrt_m_k = sqrt(m / k);

  std::vector<std::vector<std::ofstream*>> ofs(2, std::vector<std::ofstream*>(4));
  for (int model = 0; model < 2; model++) {
    for (int integrator = 0; integrator < 4; integrator++) {
      ofs[model][integrator] = new std::ofstream(filename[model][integrator]);
    }
  }

  for (int model = 0; model < 2; model++) {
    nu = model * nu0;
    for (int integrator = 0; integrator < 4; integrator++) {
      double alpha0 = 1.0;
      if (integrator == Beeman) {
        alpha0 *= 2;
      } else if (integrator == RK4) {
        alpha0 *= 4;
      }
      for (double alpha = alpha0; alpha >= 0.01; alpha *= 0.95) {
        double dt = alpha * sqrt_m_k;
        double dt_2 = 0.5 * dt;
        double dt2_2 = dt * dt_2;
        double dt2_8 = dt * dt / 8.0;
        double dt_6 = dt / 6.0;
        double dt2_6 = dt * dt / 6.0;

        x = 0.0;
        v = 0.0;
        a = acc(k, m, F0, x, v, nu);

        t = 0.0;
        double err = 0.0;
        double nb = 0.0;
        double initialTotalEnergy = calculateTotalEnergy(k, m, v, x, F0);
        double totalEnergy;
        std::vector<double> energyErrors;
        double maxError = 0.0;
        double finalError = 0.0;
        //double lastPositiveRef = 0.0;
        double prevRef = 0.0;
        double energyDrift = 0.0;

        auto start_time = std::chrono::high_resolution_clock::now();
        while (t < 20.0) {

          switch (integrator) {
            case Euler: {
              x += v * dt;
              v += a * dt;
              a = acc(k, m, F0, x, v, nu);
            } break;
            case VelocityVerlet: {
              x += v * dt + a * dt2_2;
              v += a * dt_2;
              a = acc(k, m, F0, x, v, nu);
              v += a * dt_2;
            } break;
            case Beeman: {
              aprev = a;
              a = acc(k, m, F0, x, v, nu);
              x += v * dt + (4.0 * a - aprev) * dt * dt / 6.0;
              vsave = v;
              v = v + (1.5 * a - 0.5 * aprev) * dt;  // predicted
              anext = acc(k, m, F0, x, v, nu);
              v = vsave + (5.0 * anext + 8.0 * a - aprev) * dt / 12.0;  // corrected
            } break;
            case RK4: {
              double x0 = x;
              double v0 = v;
              a = acc(k, m, F0, x, v, nu);
              double k1 = a;
              x = x0 + dt_2 * v0 + dt2_8 * k1;
              v = v0 + dt_2 * k1;
              a = acc(k, m, F0, x, v, nu);
              double k2 = a;
              v = v0 + dt_2 * k2;
              a = acc(k, m, F0, x, v, nu);
              double k3 = a;
              x = x0 + dt * v0 + dt2_2 * k3;
              v = v0 + dt * k3;
              a = acc(k, m, F0, x, v, nu);
              double k4 = a;
              x = x0 + dt * v0 + dt2_6 * (k1 + k2 + k3);
              v = v0 + dt_6 * (k1 + 2 * k2 + 2 * k3 + k4);
            } break;
          }

          t += dt;

          double ref = xsol(t, F0, k, m, nu);
          err += fabs(x - ref);
          nb += 1.0;

          // Calculate total energy and energy error
          totalEnergy = calculateTotalEnergy(k, m, v, x, F0);
          double energyError = std::abs(totalEnergy - initialTotalEnergy);
          energyErrors.push_back(energyError);
          maxError = std::max(maxError, energyError);
          
          /*
          if (ref > 0.0) {
            lastPositiveRef = ref;
          }
          if (ref <= 0.0 && lastPositiveRef > 0.0) {
            finalError = energyError;
          }
          */
          // Check if ref is crossing F0/k from above or below
          if ((ref > F0/k && prevRef <= F0/k) || (ref < F0/k && prevRef >= F0/k)) {
            finalError = energyError;
          }

          // Update prevRef to the current value of ref
          prevRef = ref;
          
          
          energyDrift += energyError;
        }
        auto end_time = std::chrono::high_resolution_clock::now();
        auto tics = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();

        // Calculate RMSE of energy errors
        double sumSquaredErrors = 0.0;
        for (double error : energyErrors) {
          sumSquaredErrors += error * error;
        }
        double rmse = std::sqrt(sumSquaredErrors / energyErrors.size());

        // Calculate average error
        double averageError = sumSquaredErrors / energyErrors.size();

        // Calculate energy drift
        energyDrift /= t;

        err /= nb;
        ofs[model][integrator]->precision(15);
        (*ofs[model][integrator]) << alpha << ' ' << err << ' ' << tics / (20.0e9) << ' ' << rmse << ' ' << maxError
                                  << ' ' << averageError << ' ' << finalError << ' ' << energyDrift << '\n';
      }

    }  // for integrator
  }  // for model

  for (int model = 0; model < 2; model++) {
    for (int integrator = 0; integrator < 4; integrator++) {
      ofs[model][integrator]->close();
      delete ofs[model][integrator];
      ofs[model][integrator] = nullptr;
    }
  }

  return 0;
}
