#pragma once
#include <cassert>
#include <cmath>
#include "BlackScholesSolver.h"

class ErrorAnalysis
{
    public:
        // Linear interpolation in x at the existing time level j.
        static double interpolate(const BlackScholesSolver& solution, double x, int j)
        {
            assert(j >= 0 && j <= solution.M);
            assert(x >= solution.x[0] && x <= solution.x[solution.N]);
            if(x == solution.x[solution.N]) return solution.value(solution.N, j);

            int left = 0;
            int right = solution.N;
            while(right - left > 1)
            {
                int middle = (left + right) / 2;
                if(solution.x[middle] <= x) left = middle;
                else right = middle;
            }
            double weight = (x - solution.x[left])
                          / (solution.x[left + 1] - solution.x[left]);
            return (1.0 - weight) * solution.value(left, j)
                 + weight * solution.value(left + 1, j);
        }

        // The paper's error: max over ALL spatial and time nodes, including j=M.
        // coarse and reference must have the same time mesh.
        static double maxError(const BlackScholesSolver& coarse,
                               const BlackScholesSolver& reference,
                               int* iMax = nullptr, int* jMax = nullptr)
        {
            assert(coarse.M == reference.M);
            assert(fabs(coarse.T - reference.T) < 1e-12);
            assert(coarse.x[0] >= reference.x[0]);
            assert(coarse.x[coarse.N] <= reference.x[reference.N]);

            double maximum = 0.0;
            int maxI = 0;
            int maxJ = 0;
            for(int j = 0; j <= coarse.M; j++)
            {
                int k = 0;
                for(int i = 0; i <= coarse.N; i++)
                {
                    while(k < reference.N - 1 && reference.x[k + 1] < coarse.x[i]) k++;
                    double weight = (coarse.x[i] - reference.x[k])
                                  / (reference.x[k + 1] - reference.x[k]);
                    double interpolated = (1.0 - weight) * reference.value(k, j)
                                        + weight * reference.value(k + 1, j);
                    double error = fabs(coarse.value(i, j) - interpolated);
                    if(error > maximum)
                    {
                        maximum = error;
                        maxI = i;
                        maxJ = j;
                    }
                }
            }
            if(iMax != nullptr) *iMax = maxI;
            if(jMax != nullptr) *jMax = maxJ;
            return maximum;
        }

        static double minimumValue(const BlackScholesSolver& solution)
        {
            double minimum = solution.U[0];
            for(int k = 1; k < (solution.N + 1) * (solution.M + 1); k++)
                if(solution.U[k] < minimum) minimum = solution.U[k];
            return minimum;
        }

        static double constantCall(double x, double t, double K, double T,
                                   double r, double sigma)
        {
            assert(x >= 0.0 && K > 0.0 && sigma > 0.0 && t <= T);
            if(x == 0.0) return 0.0;
            if(t == T)
            {
                if(x > K) return x - K;
                return 0.0;
            }

            double time = T - t;
            double d1 = (log(x / K) + (r + sigma * sigma / 2.0) * time)
                      / (sigma * sqrt(time));
            double d2 = d1 - sigma * sqrt(time);
            double Phi1 = 0.5 * erfc(-d1 / sqrt(2.0));
            double Phi2 = 0.5 * erfc(-d2 / sqrt(2.0));
            return x * Phi1 - K * exp(-r * time) * Phi2;
        }

        // Compare with the unsmoothed, infinite-domain Black-Scholes call formula.
        static double maxConstantError(const BlackScholesSolver& solution, int j = 0)
        {
            assert(j >= 0 && j <= solution.M);
            double maximum = 0.0;
            for(int i = 0; i <= solution.N; i++)
            {
                double exact = constantCall(solution.x[i], j * solution.tau,
                                           solution.K, solution.T, solution.r,
                                           solution.sigmaConstant);
                double error = fabs(solution.value(i, j) - exact);
                if(error > maximum) maximum = error;
            }
            return maximum;
        }
};
