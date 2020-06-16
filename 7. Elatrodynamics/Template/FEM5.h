/*This is a skeleton code file for use with the Finite Element Method for Problems in Physics.
  It uses the deal.II FEM library, dealii.org*/

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

//Define the order of the basis functions (Lagrange polynomials)
//and the order of the quadrature rule globally
const unsigned int order = 1;
const unsigned int quadRule = 2;

template <int dim>
class FEM
{
 public:
  //Class functions
  FEM (double Beta, double Gamma); // Class constructor 
  ~FEM(); //Class destructor

  //Solution steps
  void generate_mesh(std::vector<unsigned int> numberOfElements);
  void define_boundary_conds();
  void setup_system();
  void assemble_system();
  void apply_initial_conditions();
  void solve_trans();
  void output_trans_results(unsigned int index);

  //Function to calculate components of the elasticity tensor
  double C(unsigned int i,unsigned int j,unsigned int k,unsigned int l);

  //Class objects
  Triangulation<dim>   triangulation; //mesh
  FESystem<dim>        fe; //FE element
  DoFHandler<dim>      dof_handler; // Connectivity matrices

  QGauss<dim>  	 	 		 quadrature_formula; //Quadrature

  //Data structures
  SparsityPattern      					sparsity_pattern; //Sparse matrix pattern
  SparseMatrix<double> 					M, K, system_matrix; //Global stiffness matrix - Sparse matrix - used in the solver
  Vector<double>       					D_steady, D_trans, V_trans, A_trans, F, RHS; //Global vectors - Solution vector (D) and Global force vector (F)

  Table<2,double>								dofLocation;			//Table of the coordinates of nodes by global dof number
  std::map<unsigned int,double> boundary_values_of_D; //Map of dirichlet boundary conditions for D
  std::map<unsigned int,double> boundary_values_of_V; //Map of dirichlet boundary conditions for V
  std::map<unsigned int,double> boundary_values_of_A; //Map of dirichlet boundary conditions for Acc
  double							beta, gamma;
    
  //solution name array
  std::vector<std::string> nodal_solution_names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
};

// Class constructor for a scalar field
template <int dim>
FEM<dim>::FEM (double Beta, double Gamma)
:
//fe (FE_Q<dim>(order), dim),
fe(FE_Q<dim>(QIterated<1>(QTrapez<1>(),order)), dim),
  dof_handler (triangulation),
  quadrature_formula(quadRule)
{
  beta = Beta;
  gamma = Gamma;

  for (unsigned int i=0; i<dim; ++i){
    nodal_solution_names.push_back("u");
    nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  }
}

//Class destructor
template <int dim>
FEM<dim>::~FEM (){dof_handler.clear ();}

//Function to calculate the components of the 4th order elasticity tensor
template <int dim>
double FEM<dim>::C(unsigned int i,unsigned int j,unsigned int k,unsigned int l){

  double E=;  //EDIT
  double nu=;  //EDIT
  double lambda=(E*nu)/((1.+nu)*(1.-2.*nu)),
    mu=E/(2.*(1.+nu));

  return lambda*(i==j)*(k==l) + mu*((i==k)*(j==l) + (i==l)*(j==k));

}

//Define the problem domain and generate the mesh
template <int dim>
void FEM<dim>::generate_mesh(std::vector<unsigned int> numberOfElements){

  //Define the limits of your domain
  double x_min = , //EDIT - define the left limit of the domain, etc.
    x_max = , //EDIT
    y_min = , //EDIT
    y_max = , //EDIT
    z_min = , //EDIT
    z_max = ; //EDIT

  Point<dim,double> min(x_min,y_min,z_min),
    max(x_max,y_max,z_max);
  GridGenerator::subdivided_hyper_rectangle (triangulation, numberOfElements, min, max);
}

//Specify the Dirichlet boundary conditions
template <int dim>
void FEM<dim>::define_boundary_conds(){


  // Define the Dirichlet boundary conditions.
  // You can use the dofLocation matrix, as before.
  // See HW4 template to see how to get nodaldof, if necessary.
  // Remember to apply boundary conditions on the acceleration

  //EDIT
  
}

//Setup data structures (sparse matrix, vectors)
template <int dim>
void FEM<dim>::setup_system(){

  //Let deal.II organize degrees of freedom
  dof_handler.distribute_dofs (fe);

  //Get a vector of global degree-of-freedom x-coordinates
  MappingQ1<dim,dim> mapping;
  std::vector< Point<dim,double> > dof_coords(dof_handler.n_dofs());
  dofLocation.reinit(dof_handler.n_dofs(),dim);
  DoFTools::map_dofs_to_support_points<dim,dim>(mapping,dof_handler,dof_coords);
  for(unsigned int i=0; i<dof_coords.size(); i++){
    for(unsigned int j=0; j<dim; j++){
      dofLocation[i][j] = dof_coords[i][j];
    }
  }

  //Specify boundary condtions (call the function)
  define_boundary_conds();

  //Define the size of the global matrices and vectors
  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();
  K.reinit (sparsity_pattern);
  M.reinit (sparsity_pattern);
  system_matrix.reinit (sparsity_pattern);
  D_steady.reinit(dof_handler.n_dofs());
  D_trans.reinit(dof_handler.n_dofs());
  V_trans.reinit(dof_handler.n_dofs());
  A_trans.reinit(dof_handler.n_dofs());
  RHS.reinit(dof_handler.n_dofs());
  F.reinit(dof_handler.n_dofs());

  //Just some notes...
  std::cout << "   Number of active elems:       " << triangulation.n_active_cells() << std::endl;
  std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;   
}

//Form elmental vectors and matrices and assemble to the global vector (F) and matrix (K)
template <int dim>
void FEM<dim>::assemble_system(){

  M=0; K=0; F=0;

  FEValues<dim> fe_values(fe,
			  quadrature_formula,
			  update_values | 
			  update_gradients | 
			  update_JxW_values);

  const unsigned int 				nodes_per_elem = GeometryInfo<dim>::vertices_per_cell;
  const unsigned int   			dofs_per_elem = fe.dofs_per_cell; //This gives you dofs per element
  unsigned int 							num_quad_pts = quadrature_formula.size(); //Total number of quad points in the element
  FullMatrix<double> 				Mlocal (dofs_per_elem, dofs_per_elem);
  FullMatrix<double> 				Klocal (dofs_per_elem, dofs_per_elem);
  Vector<double>      			Flocal (dofs_per_elem);
  std::vector<unsigned int> local_dof_indices (dofs_per_elem); //This relates local dof numbering to global dof numbering
  double										rho = 7600;

  //loop over elements  
  typename DoFHandler<dim>::active_cell_iterator elem = dof_handler.begin_active (),
    endc = dof_handler.end();
  for (;elem!=endc; ++elem){

    //Retrieve the effective "connectivity matrix" for this element
    elem->get_dof_indices (local_dof_indices);

    fe_values.reinit(elem); //Retrieve values from current element
    elem->get_dof_indices (local_dof_indices);

    //Loop over local DOFs, quadrature points, etc. to populate Mlocal and Klocal
    Mlocal = 0.;
    Klocal = 0;
    
    //EDIT

    //assemble K and M from Klocal and Mlocal
    
    //EDIT
  }
}

//Apply initial conditions for the transient problem
template <int dim>
void FEM<dim>::apply_initial_conditions(){

  /*Loop over global nodes. Use dofLocation to determine the position of the node
    and add the correct initial displacement value in D_trans and V_trans.*/

  //EDIT
  
  //Find A_0 = M^{-1}*(F_0 - K*D_0) (no damping)
  
  //EDIT

  //Output initial state
  output_trans_results(0);
}

//Solve for D_transient
template <int dim>
void FEM<dim>::solve_trans(){

  //Define delta_t
  const double delta_t = 1e-6;

  const unsigned int totalDOFs = dof_handler.n_dofs(); //Total number of nodes
  Vector<double> D_tilde(totalDOFs), V_tilde(totalDOFs);

  //Loop over time steps. For each time step, update D_transient from D_n to D_{n+1} using the V method
  for(unsigned int t_step=1; t_step<1001; t_step++){
	
    //Find D_tilde and V_tilde. Remember, at this point D_trans = D_n, V_trans = V_n, and A_trans = A_n

    //EDIT
    
    //Construct system_matrix and RHS to solve for A_trans (which will then be A_{n+1})

    //EDIT
    
    //Apply boundary conditions on A_trans before solving the matrix/vector system

    //EDIT

    //Solve for A_trans (A_{n+1}) in system_matrix*A_trans= RHS

    //EDIT

    //Update D_trans to D_{n+1} using D_tilde, V_tilde, and A_trans (A_{n+1})

    //EDIT
    
    //Output the results every 100 seconds

    //EDIT
  }
}

//Output transient results for a given time step
template <int dim>
void FEM<dim>::output_trans_results (unsigned int index){
  //This adds an index to your filename so that you can distinguish between time steps

  //Write results to VTK file
  char filename[100];
  snprintf(filename, 100, "solution_%d.vtk", index);
  std::ofstream output1 (filename);
  DataOut<dim> data_out; data_out.attach_dof_handler (dof_handler);

  //Add nodal DOF data
  data_out.add_data_vector (D_trans, nodal_solution_names, DataOut<dim>::type_dof_data, nodal_data_component_interpretation);
  data_out.build_patches ();
  data_out.write_vtk (output1);
  output1.close();
}
