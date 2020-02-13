#include "mfem.hpp"
#include <iostream>
using std::cout;
using std::endl;
using std::flush;

#include <vector>
using std::vector;

using namespace mfem;


double Initial(const Vector &x);
double wrapper(const vector<double> &x, vector<double> &grad, void *other);

class Area : public Coefficient
{
private:
  GridFunction &u;
  Vector grad;
public:
  Area( GridFunction &_u ) : u(_u) { }
  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
  virtual ~Area() { }
};

class InverArea : public Coefficient
{
private:
  GridFunction &u;
  Vector grad;
public:
  InverArea( GridFunction &_u ) : u(_u) { }
  virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
  virtual ~InverArea() { }
};

class Surface
{
private:
  FiniteElementSpace *fespace;
  Array<int> &ess_bdr;
  GridFunction *y;
  Area *pp_coeff;
  InverArea *ipp_coeff;
  LinearForm *vol;
  BilinearForm *a;
  int iter;
  bool visualization;
  socketstream sout;
  Mesh *mesh;
public:
  GridFunction *u_gf;
  Surface(FiniteElementSpace *f, Array<int> &es) : 
    fespace(f), ess_bdr(es)
  {
    iter = 0;
    visualization = true;
    mesh = fespace->GetMesh();
  }
  
  ~Surface();
  void Builder();
  double myfunc(const vector<double> &x, vector<double> &grad);
};
