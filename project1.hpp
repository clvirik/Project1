#pragma once
#include <string>

class Poisson {
private:
	int m_N;
	double m_h, m_exp;
	double *m_u, *m_v, *m_x;

	double f(double x);
	double exact(double x);
	double** Dmatrix(int row, int col);
	void delete_Dmatrix(double** M, int row, int col);

public:

	void Initialize(double exp);
	void Integrate_proper();
	void Integrate_efficient();
	void LU_decomp();
	double get_maxerror();
	void Print();
	void Write(const std::string& tag = std::string());
	void Free();
};