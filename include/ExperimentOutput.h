#pragma once
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <stdexcept>
#include <string>
#include "BlackScholesSolver.h"

// Only output/plotting uses strings and filesystem; numerical arrays stay in C++.
inline void prepareOutput()
{
    std::filesystem::create_directories("../results");
}

inline FILE* openOutput(const char* path)
{
    FILE* file = fopen(path, "w");
    if(!file) throw std::runtime_error("Cannot open output file");
    return file;
}

inline void saveSolution(const BlackScholesSolver& solver, const char* path)
{
    FILE* file = openOutput(path);
    fprintf(file, "# x t U_i^j\n");
    for(int j=0; j<=solver.M; j++)
    {
        for(int i=0; i<=solver.N; i++)
            fprintf(file, "%.17g %.17g %.17g\n", solver.x[i],
                    j*solver.T/solver.M, solver.value(i,j));
        fprintf(file, "\n");
    }
    fclose(file);
}

// Run a standard plotting program. All script paths are repository-owned literals.
inline void runGnuplot(const char* script)
{
    std::string path(script);
    if(path.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_./-")
       != std::string::npos)
        throw std::invalid_argument("Use a plain relative plot script path");
    std::string command = "gnuplot '" + path + "'";
    if(std::system(command.c_str()) != 0)
        throw std::runtime_error("gnuplot failed; install gnuplot and rerun this cell");
}

inline void plotSolution(const BlackScholesSolver& solver, const char* stem, const char* title)
{
    std::string base = std::string("../results/") + stem;
    saveSolution(solver, (base + "_surface.dat").c_str());
    FILE* slices = openOutput((base + "_slices.csv").c_str());
    fprintf(slices, "x,U_t0,U_t025,U_t05,U_t075,U_T\n");
    for(int i=0; i<=solver.N; i++)
        fprintf(slices, "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n", solver.x[i],
                solver.value(i,0),solver.value(i,solver.M/4),
                solver.value(i,solver.M/2),solver.value(i,3*solver.M/4),
                solver.value(i,solver.M));
    fclose(slices);

    FILE* script = openOutput((base + ".gp").c_str());
    fprintf(script, "set terminal svg size 1000,600 font 'sans,13'\n");
    fprintf(script, "set output '%s_surface.svg'\n", base.c_str());
    fprintf(script, "set title '%s: N=M=%d'\n", title, solver.N);
    fprintf(script, "set xlabel 'x'; set ylabel 't'; set zlabel 'U_i^j' rotate by 90\n");
    fprintf(script, "set view 58,315; set xyplane at 0; set key off\n");
    fprintf(script, "set palette defined (0 '#243b83', 0.5 '#2a9d8f', 1 '#e9c46a')\n");
    fprintf(script, "splot '%s_surface.dat' using 1:2:3 with pm3d\n", base.c_str());
    fprintf(script, "set terminal pngcairo size 1600,960 font 'sans,18'\n");
    fprintf(script, "set output '%s_surface.png'; replot\n", base.c_str());
    fprintf(script, "set terminal svg size 1000,540 font 'sans,13'\n");
    fprintf(script, "set output '%s_slices.svg'\n", base.c_str());
    fprintf(script, "set title '%s: time slices near K=25'\n", title);
    fprintf(script, "set xlabel 'x'; set ylabel 'U_i^j'; set xrange [15:40]; set yrange [0:20]\n");
    fprintf(script, "set key top left; set grid; set datafile separator ','\n");
    fprintf(script, "plot '%s_slices.csv' using 1:2 with lines lw 2 title 't=0', '' using 1:4 with lines lw 2 title 't=0.5', '' using 1:6 with lines dt 2 lw 2 title 't=T'\n", base.c_str());
    fprintf(script, "set terminal pngcairo size 1600,864 font 'sans,18'\n");
    fprintf(script, "set output '%s_slices.png'; replot\n", base.c_str());
    fclose(script);
    runGnuplot((base + ".gp").c_str());
}

inline void plotConvergence(const char* stem, const char* title)
{
    std::string base = std::string("../results/") + stem;
    FILE* script = openOutput((base + "_convergence.gp").c_str());
    fprintf(script, "set terminal svg size 950,550 font 'sans,13'\n");
    fprintf(script, "set output '%s_convergence.svg'\n", base.c_str());
    fprintf(script, "set title '%s: spatial refinement, M=1024'\n", title);
    fprintf(script, "set xlabel 'N'; set ylabel 'e^{N,M} (all space-time nodes)'\n");
    fprintf(script, "set logscale xy 2; set grid; set key bottom left; set datafile separator ','\n");
    fprintf(script, "stats '%s_errors.csv' using 2 nooutput\n", base.c_str());
    fprintf(script, "anchor=STATS_min\n");
    fprintf(script, "plot '%s_errors.csv' using 1:2 with linespoints lw 2 pt 7 title 'C++ result', '' using 1:4 with points pt 6 ps 1.6 title 'Paper table', anchor*(1024.0/x)**2 with lines dt 2 title 'N^{-2} guide'\n", base.c_str());
    fprintf(script, "set terminal pngcairo size 1520,880 font 'sans,18'\n");
    fprintf(script, "set output '%s_convergence.png'; replot\n", base.c_str());
    fclose(script);
    runGnuplot((base + "_convergence.gp").c_str());
}
