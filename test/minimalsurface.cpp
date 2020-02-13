#include "minimalsurface.hpp"

// valores en la frontera para el problema
double Initial(const Vector &x)
{
   return cos(sin(x(0)*x(1)));
}

// wrapper para pasar la función a nlopt
double wrapper(const vector<double> &x, vector<double> &grad, void *other) 
{
	return reinterpret_cast<Surface*>(other)->myfunc(x, grad);
}

// Método Eval para el coeficiente (calcula $\sqrt{1+(\nabla u)^2}$)
double Area::Eval(ElementTransformation &T, const IntegrationPoint &ip)
{
   u.GetGradient(T, grad);  // grad = grad(u)
   double sig = 1.0 + grad*grad;
   return sqrt(sig);
}

// Método Eval para el coeficiente (calcula $\frac{1}{\sqrt{1+(\nabla u)^2}}$)
double InverArea::Eval(ElementTransformation &T, const IntegrationPoint &ip)
{
	u.GetGradient(T, grad);  // grad = grad(u)
	double sig = 1.0 + grad*grad;
	return 1./sqrt(sig);
}


void Surface::Builder()
{
  FunctionCoefficient coeff(Initial);
	
  y = new GridFunction(fespace);
  *y = 1.;

  u_gf = new GridFunction(fespace);
  u_gf->ProjectCoefficient(coeff);

  pp_coeff = new Area(*u_gf);
  vol = new LinearForm(fespace);
  vol->AddDomainIntegrator(new DomainLFIntegrator(*pp_coeff));

  ipp_coeff = new InverArea(*u_gf);
  a = new BilinearForm(fespace);
  a->AddDomainIntegrator(new DiffusionIntegrator(*ipp_coeff));
  if (visualization)
    {
      char vishost[] = "localhost";
      int  visport   = 19916;
      sout.open(vishost, visport);
    }

}

double Surface::myfunc(const vector<double> &x, vector<double> &grad)
{	
  GridFunction &u = *u_gf;
  for (int i = 0; i<x.size() ; i++)
    u(i) = x[i];
  vol->Update();
  vol->Assemble();
  if (!grad.empty())
    {
      GridFunction X(fespace);
      a->Update();
      a->Assemble();
      a->Mult(*u_gf,X);
      ConstantCoefficient dirichlet(0.);
      X.ProjectBdrCoefficient(dirichlet,ess_bdr);
      for (int i=0;i<X.Size();i++)
	grad[i] = X(i);
    }
   
    double C=((Vector) *y) * ((Vector) *vol);
    cout << "Iter.: " << ++iter << " Objective: " << C << endl;

  if (visualization && sout.good() && iter%50 == 0)
      {
         sout.precision(8);
         sout << "solution\n" << *mesh << u << flush
         << "window_title 'Minimal surface'" << flush
	      << "plot_caption 'Iter.: " << iter << " Obj: " << C << "'" << endl;
      }
    return C;
}

Surface::~Surface()
{
	delete y;
	delete u_gf;
	delete pp_coeff;
	delete ipp_coeff;
	delete vol;
	delete a;
}
