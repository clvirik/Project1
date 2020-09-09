#include "project1.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <string>
#include <armadillo>

using namespace std;
using namespace arma;

void Poisson::Initialize(double exp)
{
    m_exp = exp;
    m_N = round(pow(10.0, exp)) + 1; // Number of gridpoints 

    m_u = new double[m_N]; // Exact
    m_v = new double[m_N]; // Numerical
    m_x = new double[m_N]; // x-values

    //Boundary conditions
    m_u[0] = 0, m_u[m_N-1] = 0;
    m_v[0] = 0, m_v[m_N-1] = 0;
    m_x[0] = 0, m_x[m_N-1] = 1;

    m_h = 1.0/(m_N-1.0); // Step length

    // setting x-values, and the exact solution u
    for (int i = 1; i < m_N-1; i++) {
        m_x[i] = m_x[i-1] + m_h;
        m_u[i] = Poisson::exact(m_x[i]);
    }

}

double Poisson::g(double x)
{
    //Equation (16)
    return m_h*m_h*100.0*exp(-10.0*x);
}

double Poisson::exact(double x)
{
    // Equation (17)
    return 1.0 - (1.0 - exp(-10.0))*x - exp(-10.0*x);
}


void Poisson::Integrate_proper()
{
    // l and r array from equation (8)
    double* l, * r;
    l = new double[m_N-3];
    r = new double[m_N-3];

    // Matrix of the values for B[0] = d, B[1] = \breve d, B[2] = g, and B[3] = \Tilde g.
    double** B;
    B = Poisson::Dmatrix(4, m_N-2);

    for (int i = 1; i < m_N-1; i++) {
        // setting d and g values
        B[0][i-1] = 2;
        B[2][i-1] = Poisson::g(m_x[i]);
    }

    // setting \Tilde g_0 and \breve d_0 values
    B[3][0] = B[2][0];
    B[1][0] = 1.0/B[0][0];

    for (int i = 0; i < m_N-3; i++) {
        // setting l and r values.
        l[i] = -1;
        r[i] = -1;
    }

    for (int i = 1; i < m_N-2; i++) {
        // Forward sub: equation (10) and (11)
        B[1][i] = 1.0/(B[0][i] - l[i-1]*B[1][i-1]*r[i-1]);
        B[3][i] = B[2][i] - l[i-1]*B[1][i-1]*B[3][i-1];
    }

    m_v[m_N-2] = B[3][m_N-3]*B[1][m_N-3];

    for (int i = m_N-3; i > 0; i--) {
        // Backward sub: equation (12)
        m_v[i] = (B[3][i-1] - m_v[i+1]*r[i-1])*B[1][i-1];
    }

    // Freeing memory
    delete[] l;
    delete[] r;
    Poisson::delete_Dmatrix(B, 4, m_N-2);
}

void Poisson::Integrate_efficient()
{
    double d_breve;
    double* g_tilde;
    g_tilde = new double[m_N-2];
    g_tilde[0] = Poisson::g(m_x[1]);

    for (int i = 1; i < m_N-2; i++) {
        // Forward sub, equations (13) and (14)
        d_breve = i/(i+1.0);
        g_tilde[i] = Poisson::g(m_x[i+1]) + d_breve*g_tilde[i-1];
    }

    m_v[m_N-2] = g_tilde[m_N-3]*(m_N-2.0)/(m_N-1.0);

    for (int i = m_N-3; i > 0; i--) {
        // Backward sub, equations (13) and (15)
        d_breve = i/(i+1.0); 
        m_v[i] = (g_tilde[i-1] + m_v[i+1])*d_breve;
    }

    //Freeing memory
    delete[] g_tilde;
}

void Poisson::LU_decomp()
{
    // LU-decomposition
    mat A(m_N-2, m_N-2);
    vec g(m_N-2);
    
    for (int i = 0; i < m_N-2; i++) {
        g(i) = Poisson::g(m_x[i+1]);

        for (int j = 0; j < m_N-2; j++) {
            if (i == j)
                A(i, j) = 2;
            else if (i+1==j or i-1==j)
                A(i, j) = -1;
            else
                A(i, j) = 0;
        }
    }

    mat L, U;

    lu(L, U, A);

    vec w, v;

    w = solve(trimatl(L), g);
    v = solve(trimatu(U), w);

    for (int i = 0; i < m_N-2; i++) 
        m_v[i+1] = v(i);
}

double Poisson::get_maxerror()
{
    // Using the last term as the largest error
    // (could also have checked all the terms)

    return log10(abs(1-m_v[m_N-2]/m_u[m_N-2]));

    // Freeing memory
    Poisson::Free();
}

void Poisson::Print()
{
    // Print function for convenience
    for (int i = 1; i < m_N-1; i++)
        cout << m_v[i] << ", " << m_u[i] << " " << log10(abs(1-m_v[i]/m_u[i])) << endl;

    // Freeing memory
    Poisson::Free();
}

void Poisson::Write(const string& tag)
{
    // Write function to plot in python
    string name = "N";
    name += to_string((int)m_exp);
    name += tag;
    name += ".txt";

    ofstream outFile;
    outFile.open(name);

    for (int i = 0; i <= m_N - 1; i++) 
        outFile << setw(15) << setprecision(8) << scientific << m_x[i] << " " << m_v[i] << " " << m_u[i] << " " << log10(abs(1-m_v[i]/m_u[i])) << endl;

    outFile.close();

    // Freeing memory
    Poisson::Free();
}

void Poisson::Free()
{
    delete[] m_x;
    delete[] m_v;
    delete[] m_u;
}

double** Poisson::Dmatrix(int row, int col)
{
    double** M = new double* [row];
    for (int i = 0; i < row; i++) 
        M[i] = new double[col];

    return M;
}

void Poisson::delete_Dmatrix(double** M, int row, int col)
{
    for (int i = 0; i < row; i++)
       delete[] M[i];

    delete[] M;
}