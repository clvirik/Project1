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
    m_N = round(pow(10.0, exp)) + 1;
    m_u = new double[m_N];
    m_v = new double[m_N];
    m_x = new double[m_N];

    m_u[0] = 0, m_u[m_N-1] = 0;
    m_v[0] = 0, m_v[m_N-1] = 0;
    m_x[0] = 0, m_x[m_N-1] = 1;

    m_h = 1.0/(m_N-1.0);

    for (int i = 1; i < m_N-1; i++) {
        m_x[i] = m_x[i-1] + m_h;
        m_u[i] = Poisson::exact(m_x[i]);
    }

}

double Poisson::f(double x)
{
    return m_h*m_h*100.0*exp(-10.0*x);
}

double Poisson::exact(double x)
{
    return 1.0 - (1.0 - exp(-10.0))*x - exp(-10.0*x);
}


void Poisson::Integrate_proper()
{
    double* a, * c;
    a = new double[m_N-3];
    c = new double[m_N-3];

    double** B;
    B = Poisson::Dmatrix(4, m_N-2);

    for (int i = 1; i < m_N-1; i++) {
        B[0][i-1] = 2;
        B[2][i-1] = Poisson::f(m_x[i]);
    }

    B[3][0] = B[2][0];
    B[1][0] = 1.0/B[0][0];

    for (int i = 0; i < m_N-3; i++) {
        a[i] = -1;
        c[i] = -1;
    }

    for (int i = 1; i < m_N-2; i++) {
        B[1][i] = 1.0/(B[0][i] - a[i-1]*B[1][i-1]*c[i-1]);
        B[3][i] = B[2][i] - a[i-1]*B[1][i-1]*B[3][i-1];
    }

    m_v[m_N-2] = B[3][m_N-3]*B[1][m_N-3];

    for (int i = m_N-3; i > 0; i--) {
        m_v[i] = (B[3][i-1] - m_v[i+1]*c[i-1])*B[1][i-1];
    }

    delete[] a;
    delete[] c;
    Poisson::delete_Dmatrix(B, 4, m_N-2);
}

void Poisson::Integrate_efficient()
{
    double d_tilda;
    double* B;
    B = new double[m_N-2];
    B[0] = Poisson::f(m_x[1]);

    for (int i = 1; i < m_N-2; i++) {
        d_tilda = i/(i+1.0);
        B[i] = Poisson::f(m_x[i+1]) + d_tilda*B[i-1];
    }

    m_v[m_N-2] = B[m_N-3]*(m_N-2.0)/(m_N-1.0);

    for (int i = m_N-3; i > 0; i--) {
        d_tilda = i/(i+1.0);
        m_v[i] = (B[i-1] + m_v[i+1])*d_tilda;
    }

    delete[] B;
}

void Poisson::LU_decomp()
{
    mat A(m_N-2, m_N-2);
    vec g(m_N-2);
    
    for (int i = 0; i < m_N-2; i++) {
        g(i) = Poisson::f(m_x[i+1]);

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
    return log10(abs(1-m_v[m_N-2]/m_u[m_N-2]));

    Poisson::Free();
}

void Poisson::Print()
{
    for (int i = 1; i < m_N-1; i++)
        cout << m_v[i] << ", " << m_u[i] << " " << log10(abs(1-m_v[i]/m_u[i])) << endl;

    Poisson::Free();
}

void Poisson::Write(const string& tag)
{
    string name = "N";
    name += to_string((int)m_exp);
    name += tag;
    name += ".txt";

    ofstream outFile;
    outFile.open(name);

    for (int i = 0; i <= m_N - 1; i++) 
        outFile << setw(15) << setprecision(8) << scientific << m_x[i] << " " << m_v[i] << " " << m_u[i] << " " << log10(abs(1-m_v[i]/m_u[i])) << endl;

    outFile.close();

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