#ifndef _CPP_CHEM_H_
#define _CPP_CHEM_H_

#include <vector>
#include <stdexcept>

using std::vector;


class ReactionDiffusion
{
private:
  int * coeff_reac;
  int * coeff_prod;
  int * coeff_totl;
  void _fill_local_r(double*, double*);
  double * _get_p(int, double*);
  double * _get_n(int, double*);
  

public:
  int n, N, nr, mode, nfeval, njeval;
  vector<vector<int> > stoich_reac;
  vector<vector<int> > stoich_prod;
  vector<double> k;
  vector<double> D;

  ReactionDiffusion(int, int, vector<vector<int> >, vector<vector<int> >, 
		    vector<double>, vector<double>, int);
  ~ReactionDiffusion();
  void f(double, double*, double*);
  int banded_jaclu_solve(double, double *, double *, double, int);
  int dense_jaclu_solve(double, double *, double *, double, int);
  void dense_jac_rmaj(double, double*, double*, int, double, int);
  void dense_jac_cmaj(double, double*, double*, int, double, int);
  void banded_jac_cmaj(double, double*, double*, int, double, int);
  void banded_packed_jac_cmaj(double, double*, double*, int, double, int);

};
#endif
