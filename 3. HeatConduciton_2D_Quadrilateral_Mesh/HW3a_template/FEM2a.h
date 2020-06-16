/*This is a template file for use with 2D finite elements (scalar field).
  The portions of the code you need to fill in are marked with the comment "//EDIT".

  Do not change the name of any existing functions, but feel free
  to create additional functions, variables, and constants.
  It uses the deal.II FEM library.*/

//Include files
//Data structures and solvers
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
//Mesh related classes
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
//Finite element implementation classes
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
//Standard C++ libraries
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace dealii;

template <int dim>
class FEM
{
 public:
  //Class functions
  FEM(unsigned int problem); // Class constructor 
  ~FEM();									 //Class destructor

  //Define your 2D basis functions and derivatives
  double basis_function(unsigned int node, 
			double xi_1,
			double xi_2);
  std::vector<double> basis_gradient(unsigned int node, 
				     double xi_1,
				     double xi_2);

  //Solution steps
  void generate_mesh(std::vector<unsigned int> numberOfElements);
  void define_boundary_conds();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results();

		//Function to calculate the l2 norm of the error in the finite element sol'n vs. the exact solution (problem 2)
		double l2norm_of_error();

  //Class objects
  Triangulation<dim>   triangulation; //mesh
  FESystem<dim>        fe;	      //FE element
  DoFHandler<dim>      dof_handler;   //Connectivity matrices

  //Gaussian quadrature - These will be defined in setup_system()
  unsigned int	      quadRule;    //quadrature rule, i.e. number of quadrature points
  std::vector<double> quad_points; //vector of Gauss quadrature points
  std::vector<double> quad_weight; //vector of the quadrature point weights
    
  //Data structures
  SparsityPattern      	        sparsity_pattern; //Sparse matrix pattern
  SparseMatrix<double>    	K;                //Global stiffness (sparse) matrix
  Vector<double>                D, F;             //Global vectors - Solution vector (D) and Global force vector (F)
  Table<2,double>	        nodeLocation;	  //Table of the coordinates of nodes by global dof number
  std::map<unsigned int,double> boundary_values;  //Map of dirichlet boundary conditions 
	double prob;

  //solution name array
  std::vector<std::string> nodal_solution_names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
};

// Class constructor for a scalar field
template <int dim>
FEM<dim>::FEM (unsigned int problem)
:
fe(FE_Q<dim>(QIterated<1>(QTrapez<1>(),1)), 1),
  dof_handler (triangulation)
{
	prob = problem;

  //Nodal Solution names - this is for writing the output file
  nodal_solution_names.push_back("D");
  nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
}

//Class destructor
template <int dim>
FEM<dim>::~FEM (){
  dof_handler.clear ();
}

//Define basis functions
template <int dim>
double FEM<dim>::basis_function(unsigned int node, double xi_1, double xi_2){
  /*"node" specifies which node the basis function corresponds to, 
    "xi" is the point (in the bi-unit domain) where the function is being evaluated.
    You need to calculate the value of the specified basis function and order at the given quadrature pt.*/

  double value = 0.; //Store the value of the basis function in this variable

  //EDIT

  return value;
}

//Define basis function gradient
template <int dim>
std::vector<double> FEM<dim>::basis_gradient(unsigned int node, double xi_1, double xi_2){
  /*"node" specifies which node the basis function corresponds to, 
    "xi" is the point (in the bi-unit domain) where the function is being evaluated.
    You need to calculate the value of the derivative of the specified basis function and order at the given quadrature pt.
    Note that this is the derivative with respect to xi (not x)*/

  std::vector<double> values(dim,0.0); //Store the value of the gradient of the basis function in this variable

  //EDIT

  return values;
}

//Define the problem domain and generate the mesh
template <int dim>
void FEM<dim>::generate_mesh(std::vector<unsigned int> numberOfElements){

  //Define the limits of your domain
  double x_min = , //EDIT - define the left limit of the domain, etc.
    x_max = , //EDIT
    y_min = , //EDIT
    y_max = ; //EDIT

  Point<dim,double> min(x_min,y_min),
    max(x_max,y_max);
  GridGenerator::subdivided_hyper_rectangle (triangulation, numberOfElements, min, max);
}

//Specify the Dirichlet boundary conditions
template <int dim>
void FEM<dim>::define_boundary_conds(){
	
  /*Note: this will be very similiar to the define_boundary_conds function
    in the HW2 template. You will loop over all nodes and use "nodeLocations"
    to check if the node is on the boundary with a Dirichlet condition. If it is,
    then add the node number and the specified value (temperature in this problem)
    to the boundary values map, something like this:

    boundary_values[globalNodeIndex] = dirichletTemperatureValue

    Note that "nodeLocation" is now a Table instead of just a vector. The row index is
    the global node number; the column index refers to the x or y component (0 or 1 for 2D).
    e.g. nodeLocation[7][1] is the y coordinate of global node 7

		Problem 1 and problem 2 have different Dirichlet boundary conditions.*/

  const unsigned int totalNodes = dof_handler.n_dofs(); //Total number of nodes

	//Identify dirichlet boundary nodes and specify their values.
	if(prob == 1){
  	//EDIT - Define the Dirichlet boundary conditions.

	}
	else if(prob == 2){
  	//EDIT - Define the Dirichlet boundary conditions.
	}
}

//Setup data structures (sparse matrix, vectors)
template <int dim>
void FEM<dim>::setup_system(){

  //Let deal.II organize degrees of freedom
  dof_handler.distribute_dofs (fe);

  //Fill in the Table "nodeLocations" with the x and y coordinates of each node by its global index
  MappingQ1<dim,dim> mapping;
  std::vector< Point<dim,double> > dof_coords(dof_handler.n_dofs());
  nodeLocation.reinit(dof_handler.n_dofs(),dim);
  DoFTools::map_dofs_to_support_points<dim,dim>(mapping,dof_handler,dof_coords);
  for(unsigned int i=0; i<dof_coords.size(); i++){
    for(unsigned int j=0; j<dim; j++){
      nodeLocation[i][j] = dof_coords[i][j];
    }
  }

  //Specify boundary condtions (call the function)
  define_boundary_conds();

  //Define the size of the global matrices and vectors
  sparsity_pattern.reinit (dof_handler.n_dofs(), dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();
  K.reinit (sparsity_pattern);
  F.reinit (dof_handler.n_dofs());
  D.reinit (dof_handler.n_dofs());

  //Define quadrature rule - again, you decide what quad rule is needed
  quadRule = 2; //EDIT - Number of quadrature points along one dimension
  quad_points.resize(quadRule); quad_weight.resize(quadRule);

  quad_points[0] = -sqrt(1./3.); //EDIT
  quad_points[1] = sqrt(1./3.); //EDIT

  quad_weight[0] = 1.; //EDIT
  quad_weight[1] = 1.; //EDIT


	//EDIT you need "xi_at_node" like firt coding assignment, but in 2D. You can define a function, a tensor etc. to look up the value.

  //Just some notes...
  std::cout << "   Number of active elems:       " << triangulation.n_active_cells() << std::endl;
  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;   
}

//Form elmental vectors and matrices and assemble to the global vector (F) and matrix (K)
template <int dim>
void FEM<dim>::assemble_system(){

  K=0; F=0;

  const unsigned int   	    dofs_per_elem = fe.dofs_per_cell; //This gives you number of degrees of freedom per element
  FullMatrix<double> 	    Klocal (dofs_per_elem, dofs_per_elem);
  Vector<double>      	    Flocal (dofs_per_elem);
  std::vector<unsigned int> local_dof_indices (dofs_per_elem);

  //loop over elements  
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active(), 
    endc = dof_handler.end();
  for (;elem!=endc; ++elem){

    /*Retrieve the effective "connectivity matrix" for this element
      "local_dof_indices" relates local dofs to global dofs,
      i.e. local_dof_indices[i] gives the global dof number for local dof i.*/
    elem->get_dof_indices (local_dof_indices);

    //Loop over local DOFs and quadrature points to populate Flocal and Klocal.
    FullMatrix<double> Jacobian(dim,dim);
    double detJ, f = 0.;

    
    Flocal = 0.;
		
		//Loop over quadrature  points and local DOFs to populate Jacobian and Flocal 
		// you may use detJ = Jacobian.determinant() inside you loop.
		// you can check FullMatrix in deal.ii for syntax.
		
		

    //Loop over local DOFs and quadrature points to populate Klocal
		//now you need inverse of Jacobian
		//you can look at HW2, and it should be similar here
		
    FullMatrix<double> invJacob(dim,dim), kappa(dim,dim);

    //"kappa" is the conductivity tensor
    kappa = 0.;
    kappa[0][0] = 385.;
    kappa[1][1] = 385.;

    Klocal = 0.;


    //Assemble local K and F into global K and F


  }

  //Apply Dirichlet boundary conditions
  MatrixTools::apply_boundary_values (boundary_values, K, D, F, false);
}

//Solve for D in KD=F
template <int dim>
void FEM<dim>::solve(){

  //Solve for D
  SparseDirectUMFPACK  A;
  A.initialize(K);
  A.vmult (D, F); //D=K^{-1}*F

}

//Output results
template <int dim>
void FEM<dim>::output_results (){

  //Write results to VTK file
  std::ofstream output1("solution.vtk");
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);

  //Add nodal DOF data
  data_out.add_data_vector(D, nodal_solution_names, DataOut<dim>::type_dof_data,
			   nodal_data_component_interpretation);
  data_out.build_patches();
  data_out.write_vtk(output1);
  output1.close();
}

template <int dim>
double FEM<dim>::l2norm_of_error(){
	
	double l2norm = 0.;

	//Find the l2 norm of the error between the finite element sol'n and the exact sol'n
	//(For problem 2 only)
	const unsigned int   			dofs_per_elem = fe.dofs_per_cell; //This gives you dofs per element
	std::vector<unsigned int> local_dof_indices (dofs_per_elem);
	double u_exact, u_h, x, y;
	FullMatrix<double> Jacobian(dim,dim);
	double detJ;

	//loop over elements  
	typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active (), 
																								 endc = dof_handler.end();
	for (;elem!=endc; ++elem){

		//Retrieve the effective "connectivity matrix" for this element
		elem->get_dof_indices (local_dof_indices);
		for(unsigned int q1=0; q1<quadRule; q1++){
			for(unsigned int q2=0; q2<quadRule; q2++){
				Jacobian = 0.;
				/*EDIT - you'll need to find the determinant of the Jacobian to perform the
				numerical integration (see, for example, Flocal or Klocal)*/
				x = 0.; y = 0.; u_h = 0.;
				for(unsigned int B=0; B<dofs_per_elem; B++){
					x += nodeLocation[local_dof_indices[B]][0]*basis_function(B,quad_points[q1],quad_points[q2]);
					//EDIT - You'll also need y and u_h (the finite element solution) at the quadrature points.
				}
		    //EDIT - Find the l2-norm of the error through numerical integration - for problem 2 only.
		    /*This includes evaluating the exact solution at the quadrature points*/
			}
		}
	}

	return sqrt(l2norm);
}