/*This is a skeleton code file for use with the Finite Element Method for Problems in Physics.
  It uses the deal.II FEM library, dealii.org*/

//Include files
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "FEM5.h"

using namespace dealii;

//The main program, using the FEM class
int main (){
  try{
    deallog.depth_console (0);

    const int dimension = 3;
    double beta = 0.25; //Specify beta and gamma (Newmark family) 
    double gamma = 0.5;
    FEM<dimension> problemObject(beta,gamma);
    
    //NOTE: This is where you define the number of elements in the mesh
    std::vector<unsigned int> num_of_elems(dimension);
    num_of_elems[0] = 5;
    num_of_elems[1] = 5;
    num_of_elems[2] = 5; //For example, a 5x5x5 mesh
    
    problemObject.generate_mesh(num_of_elems);
    problemObject.setup_system();
    problemObject.assemble_system();
    problemObject.apply_initial_conditions();
    problemObject.solve_trans();
        
  }
  catch (std::exception &exc){
    std::cerr << std::endl << std::endl
	      << "----------------------------------------------------"
	      << std::endl;
    std::cerr << "Exception on processing: " << std::endl
	      << exc.what() << std::endl
	      << "Aborting!" << std::endl
	      << "----------------------------------------------------"
	      << std::endl;

    return 1;
  }
  catch (...){
    std::cerr << std::endl << std::endl
	      << "----------------------------------------------------"
	      << std::endl;
    std::cerr << "Unknown exception!" << std::endl
	      << "Aborting!" << std::endl
	      << "----------------------------------------------------"
	      << std::endl;
    return 1;
  }

  return 0;
}
