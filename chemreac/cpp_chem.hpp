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
  void _get_dx(int, double *, double *, double *);
  

public:
  int n, N, nr, mode, nfeval, njeval;
  vector<vector<int> > stoich_reac;
  vector<vector<int> > stoich_actv;
  vector<vector<int> > stoich_prod;
  vector<double> k;
  vector<double> D;
  vector<double> x; // separation

  ReactionDiffusion(int, int,
		    vector<vector<int> >, 
		    vector<vector<int> >, 
		    vector<vector<int> >, 
		    vector<double>, 
		    vector<double>, 
		    vector<double>,
		    int);
  ~ReactionDiffusion();
  void f(double, double*, double*);
  void dense_jac_rmaj(double, double*, double*, int);
  void dense_jac_cmaj(double, double*, double*, int);
  void banded_jac_cmaj(double, double*, double*, int);
  void banded_packed_jac_cmaj(double, double*, double*, int);

};
#endif
