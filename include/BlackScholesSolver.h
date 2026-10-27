#pragma once
#include <cassert>
#include <cmath>
#include "TridiagonalSolver.h"

// N spatial intervals, M time intervals; U_i^j is U[j*(N+1)+i].
// The terminal condition is at j=M, and stepping proceeds toward j=0.
class BlackScholesSolver
{
    public: // Member
        int N;
        int M;
        double alpha;
        double betaStar;
        double epsilon;
        double T;
        double K;
        double Smax;
        double tau;
        double baseH; // paper mesh parameter h; first gap for the resolved mesh
        double r;
        double sigmaConstant;

        double* x;
        double* h; // h[i] = x[i]-x[i-1], with h[0]=0
        double* U;
        double* A;
        double* B;
        double* C;
        double* sigmaSquared;

    private:
        TridiagonalSolver solver;

    public: // Constructor
        BlackScholesSolver(int N_, int M_, double alpha_, double betaStar_,
                           double epsilon_ = 1e-4, double T_ = 1.0,
                           double K_ = 25.0, double Smax_ = 100.0)
            : solver(N_ - 1)
        {
            assert(N_ >= 8 && N_ % 4 == 0 && M_ > 0);
            assert(alpha_ > 0.0 && betaStar_ > 0.0 && T_ > 0.0);
            assert(epsilon_ > 0.0 && epsilon_ < K_);
            assert(K_ + epsilon_ < Smax_);

            N = N_;
            M = M_;
            alpha = alpha_;
            betaStar = betaStar_;
            epsilon = epsilon_;
            T = T_;
            K = K_;
            Smax = Smax_;
            tau = T / M;
            r = 0.06;
            sigmaConstant = 0.2;

            x = new double[N + 1];
            h = new double[N + 1];
            U = new double[(N + 1) * (M + 1)]{};
            A = new double[N + 1]{};
            B = new double[N + 1]{};
            C = new double[N + 1]{};
            sigmaSquared = new double[N + 1]{};
            generatePiecewiseMesh();
        }

        ~BlackScholesSolver()
        {
            delete[] x;
            delete[] h;
            delete[] U;
            delete[] A;
            delete[] B;
            delete[] C;
            delete[] sigmaSquared;
        }

        BlackScholesSolver(const BlackScholesSolver&) = delete;
        BlackScholesSolver& operator=(const BlackScholesSolver&) = delete;

    public:
        double& value(int i, int j)
        {
            return U[j * (N + 1) + i];
        }

        double value(int i, int j) const
        {
            return U[j * (N + 1) + i];
        }

        void generatePiecewiseMesh()
        {
            int quarter = N / 4;
            baseH = (K - epsilon) / (1.0 + alpha / betaStar * (quarter - 2));
            x[0] = 0.0;
            for(int i = 1; i < quarter; i++)
                x[i] = baseH * (1.0 + alpha / betaStar * (i - 1));

            x[quarter - 1] = K - epsilon;
            x[quarter] = K;
            x[quarter + 1] = K + epsilon;
            for(int i = quarter + 2; i <= N; i++)
                x[i] = K + epsilon + (Smax - K - epsilon)
                     * (i - quarter - 1) / (3 * N / 4 - 1);
            x[N] = Smax;

            h[0] = 0.0;
            for(int i = 1; i <= N; i++) h[i] = x[i] - x[i - 1];
        }

        void generateUniformMesh()
        {
            baseH = Smax / N;
            for(int i = 0; i <= N; i++) x[i] = i * baseH;
            x[N] = Smax;
            h[0] = 0.0;
            for(int i = 1; i <= N; i++) h[i] = x[i] - x[i - 1];
        }

        // All three regions refine with their interval counts. The transition
        // half-width is independent of epsilon, which remains the payoff width.
        void generateResolvedPiecewiseMesh(double transitionHalfWidth,
                                           int leftIntervals, int transitionIntervals)
        {
            assert(transitionHalfWidth > 0.0 && transitionHalfWidth < K);
            assert(K + transitionHalfWidth < Smax);
            assert(leftIntervals > 0 && leftIntervals < N);
            assert(transitionIntervals > 0 && transitionIntervals < N - leftIntervals);
            assert(transitionIntervals % 2 == 0);

            int rightIntervals = N - leftIntervals - transitionIntervals;
            double leftEnd = K - transitionHalfWidth;
            double rightStart = K + transitionHalfWidth;
            for(int i = 0; i <= leftIntervals; i++)
                x[i] = leftEnd * i / leftIntervals;
            x[leftIntervals] = leftEnd;
            for(int i = 1; i <= transitionIntervals; i++)
                x[leftIntervals + i] = leftEnd
                    + 2.0 * transitionHalfWidth * i / transitionIntervals;
            x[leftIntervals + transitionIntervals / 2] = K;
            x[leftIntervals + transitionIntervals] = rightStart;
            for(int i = 1; i <= rightIntervals; i++)
                x[leftIntervals + transitionIntervals + i] = rightStart
                    + (Smax - rightStart) * i / rightIntervals;
            x[N] = Smax;

            h[0] = 0.0;
            for(int i = 1; i <= N; i++) h[i] = x[i] - x[i - 1];
            baseH = h[1];
        }

        void generateResolvedPiecewiseMesh(double transitionHalfWidth)
        {
            assert(N % 16 == 0);
            generateResolvedPiecewiseMesh(transitionHalfWidth, N / 4, N / 8);
        }

        double payoff(double y) const
        {
            if(y <= -epsilon) return 0.0;
            if(y >= epsilon) return y;

            // The paper's C^4 polynomial, evaluated with z=y/epsilon.
            double z = y / epsilon;
            double z2 = z * z;
            return epsilon * (35.0 / 256.0 + z / 2.0
                 + z2 * (35.0 / 64.0 + z2 * (-35.0 / 128.0
                 + z2 * (7.0 / 64.0 - 5.0 / 256.0 * z2))));
        }

        void setTerminalPayoff()
        {
            for(int i = 0; i <= N; i++) value(i, M) = payoff(x[i] - K);

            // Both experiments use the constant interest rate r=0.06.
            for(int j = 0; j <= M; j++)
            {
                double t_j = j * tau;
                value(0, j) = 0.0;
                value(N, j) = Smax - K * exp(-r * (T - t_j));
            }
        }

        void assembleCoefficients(int j)
        {
            assert(j >= 0 && j < M);
            for(int i = 1; i < N; i++)
            {
                double sumH = h[i] + h[i + 1];
                A[i] = -x[i] / sumH * (sigmaSquared[i] * x[i] / h[i] - r);
                B[i] = 1.0 / tau + sigmaSquared[i] * x[i] * x[i]
                     / (h[i] * h[i + 1]) + r;
                C[i] = -x[i] / sumH * (sigmaSquared[i] * x[i] / h[i + 1] + r);
            }
        }

        void step(int j)
        {
            assert(j >= 0 && j < M);
            for(int i = 1; i < N; i++)
            {
                solver.a[i - 1] = A[i];
                solver.b[i - 1] = B[i];
                solver.c[i - 1] = C[i];
                solver.d[i - 1] = value(i, j + 1) / tau;
            }

            // Boundary nodes are known; only i=1,...,N-1 are unknowns.
            solver.d[0] -= A[1] * value(0, j);
            solver.d[N - 2] -= C[N - 1] * value(N, j);
            solver.a[0] = 0.0;
            solver.c[N - 2] = 0.0;
            solver.solve(U + j * (N + 1) + 1);
        }

        void solveTest1()
        {
            setTerminalPayoff();
            for(int j = M - 1; j >= 0; j--)
            {
                double t_j = j * tau;
                for(int i = 0; i <= N; i++)
                {
                    double sigma = 0.2 * (1.0 + t_j * exp(-x[i]));
                    sigmaSquared[i] = sigma * sigma;
                }
                assembleCoefficients(j);
                step(j);
            }
        }

        void solveTest2()
        {
            setTerminalPayoff();
            for(int i = 0; i <= N; i++)
            {
                double sigma = 0.4 * (2.0 + sin(x[i]));
                sigmaSquared[i] = sigma * sigma;
            }
            for(int j = M - 1; j >= 0; j--)
            {
                assembleCoefficients(j);
                step(j);
            }
        }

        void solveConstant()
        {
            setTerminalPayoff();
            for(int i = 0; i <= N; i++)
                sigmaSquared[i] = sigmaConstant * sigmaConstant;

            for(int j = M - 1; j >= 0; j--)
            {
                assembleCoefficients(j);
                step(j);
            }
        }

        // These diagnostics use the currently assembled time level.
        double largestInteriorOffDiagonal() const
        {
            double largest = A[2];
            for(int i = 1; i < N; i++)
            {
                if(i >= 2 && A[i] > largest) largest = A[i];
                if(i <= N - 2 && C[i] > largest) largest = C[i];
            }
            return largest;
        }

        double minimumDiagonalMargin() const
        {
            double minimum = B[1] - fabs(C[1]);
            for(int i = 1; i < N; i++)
            {
                double margin = B[i];
                if(i >= 2) margin -= fabs(A[i]);
                if(i <= N - 2) margin -= fabs(C[i]);
                if(margin < minimum) minimum = margin;
            }
            return minimum;
        }
};
