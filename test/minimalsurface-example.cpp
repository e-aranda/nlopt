// Calculo de superficies mínimas
// Se minimiza el funcional $$J(u) = \int_\Omega \sqrt{1+|\nabla u |^2}\,dx$$
// con $u=g$ en $\partial\Omega$
// Se utiliza nlopt-> MMA como algoritmo de minimización, para lo cual hay 
// que definir una función objetivo en la que se calcula tanto la integral anterior
// como la derivada que es
// \delta J(u;v) = \int_\Omega \frac{1}{\sqrt{1+|\nabla u|^2}}\nabla u\cdot \nabla v\,dx$$
// con $v=0$ en $\partial \Omega$

// La clase Surface es la encargada de construir el funcional y su derivada
// Se inicializan los objetos compartidos en Builder y se pasan a myfunc para calcular
// el valor de la función y su derivada
// El cálculo de la función se hace con una forma lineal en la que se pasa un 
// Coefficient calculado mediante la clase Area (básicamente calcula $\sqrt{1+|\nabla u|^2}$)
// y el valor de la integral es el producto del vector calculado por la forma lineal y
// un vector de 1
// Para calcular la derivada, al no haber una forma lineal que use $\nabla v$ usamos una
// forma bilineal de difussion con coeficiente $\frac{1}{\sqrt{1+|\nabla u|^2}}$ (construido
// mediante InverArea) y luego multiplicamos por el vector $u$. Además hay que anular los
// valores en la frontera (mediante un vector ess_bdr construido en el programa principal)
// que se pasa en el constructor de la clase
// Para nlopt es necesario pasar una función wrapper como función objetivo

#include "minimalsurface.hpp"
#include <nlopt.hpp>

int main(int argc, char *argv[])
{
	unsigned t0,t1;  
	t0 = clock(); 
	const char *mesh_file = "cuadrado.mesh";

	Mesh *mesh = new Mesh(mesh_file, 1, 1);
	int dim = mesh->Dimension();
	// Nodos en la frontera
	Array<int> ess_bdr(mesh->bdr_attributes.Max());
	ess_bdr = 1;

	FiniteElementCollection *fec;
	fec = new H1_FECollection(1,dim);
	FiniteElementSpace *fespace = new FiniteElementSpace(mesh, fec);
    
	// inicializamos objeto
	Surface obj(fespace,ess_bdr);

	unsigned n = fespace->GetTrueVSize(); 

	nlopt::opt opt(nlopt::LD_MMA, n);
	opt.set_min_objective(wrapper, &obj);
	opt.set_ftol_abs(1.e-6);

	obj.Builder();  
	GridFunction &u = *(obj.u_gf);
	vector<double> x(n);
	for (int i = 0; i<n ; i++)
		x[i] = u(i);
	double minf;
	nlopt::result result;
    
	try
	  {
	    result = opt.optimize(x, minf);
	    cout << "Resultado: " << result << "  Valor mínimo: " <<minf << endl;  
	  }
	catch(std::exception &e)
	  {
	    std::cout << e.what() << std::endl;
	  }

	delete fespace;
	delete fec;
	delete mesh;

	t1 = clock();
	double time = (double(t1-t0)/CLOCKS_PER_SEC);
	cout << "Execution Time: " << time << endl;
	return 0;
}
