#pragma once
#include <cassert>

// a: lower diagonal, b: diagonal, c: upper diagonal, d: right-hand side.
class TridiagonalSolver
{
    public: // Member
        int n;
        double* a;
        double* b;
        double* c;
        double* d;

    private:
        double* cp;
        double* dp;

    public: // Constructor
        TridiagonalSolver(int n_)
        {
            assert(n_ > 0);
            n = n_;
            a = new double[n];
            b = new double[n];
            c = new double[n];
            d = new double[n];
            cp = new double[n];
            dp = new double[n];
        }

        ~TridiagonalSolver()
        {
            delete[] a;
            delete[] b;
            delete[] c;
            delete[] d;
            delete[] cp;
            delete[] dp;
        }

        // These classes own their arrays and must not be copied.
        TridiagonalSolver(const TridiagonalSolver&) = delete;
        TridiagonalSolver& operator=(const TridiagonalSolver&) = delete;

    public:
        void setSystem(const double* a_, const double* b_,
                       const double* c_, const double* d_)
        {
            for(int i = 0; i < n; i++)
            {
                a[i] = a_[i];
                b[i] = b_[i];
                c[i] = c_[i];
                d[i] = d_[i];
            }
        }

        void solve(double* x)
        {
            // Forward elimination. Scratch arrays are reused at every time step.
            assert(b[0] != 0.0);
            cp[0] = c[0] / b[0];
            dp[0] = d[0] / b[0];

            for(int i = 1; i < n; i++)
            {
                double deno = b[i] - a[i] * cp[i - 1];
                assert(deno != 0.0);
                cp[i] = 0.0;
                if(i < n - 1) cp[i] = c[i] / deno;
                dp[i] = (d[i] - a[i] * dp[i - 1]) / deno;
            }

            // Back substitution.
            x[n - 1] = dp[n - 1];
            for(int i = n - 2; i >= 0; i--)
                x[i] = dp[i] - cp[i] * x[i + 1];
        }
};
