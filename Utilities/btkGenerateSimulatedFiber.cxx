/*==========================================================================
  
  © Université de Strasbourg - Centre National de la Recherche Scientifique
  
  Date: 10/01/2013 (last updated: 23/05/2013)
  Author(s): Larbi Boubchir (boubchir at unistra dot fr)
  
  This software is governed by the CeCILL-B license under French law and
  abiding by the rules of distribution of free software.  You can  use,
  modify and/ or redistribute the software under the terms of the CeCILL-B
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and,  more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL-B license and that you accept its terms.
  
==========================================================================*/

// library includes
#include <vtkVersion.h>
#include <vnl/vnl_vector.h>
#include <sstream>
#include <cstdlib>
 
// ---------------------------------
// Usage : btkGenerateSimulatedFiber arg1 arg2 arg3 arg4 arg5
//
// 	Input argument:
// 	[arg1]: fiber number
//      [arg2,arg3,arg4]: incrementation step ({arg2,arg3,arg4} will be added, respectively, to the coordinates X Y Z of each point)
//      [arg5]: new fiber number for the generated fiber
//
// 	Output:
// 	Generate a new simulated fiber. The data is stored in a text file namely 'Simulated_Fiber%.txt' where '%' corresponds to [arg5] 
//
// NOTE: each input fiber data is stored into a text file (.txt) containing the coordinates X Y Z of all the points of the fiber.
// each text file is named 'Fiber%.txt' where '%' corresponds to [arg1].


// Example : btkGenerateSimulatedFiber 4 -1 -2 -3 1
//
// In this example, the code add {-1,-2,-3} to the coordinates {X,Y,Z} of each point of the input fiber 'Fiber4.txt' and then generate the  new simulated fiber 'Simulated_Fiber1.txt' 
// 
// ---------------------------------


//
// main
int main(int argc, char* argv[])
{

// Read input arguments
int filenumber = atoi( argv[1] );
double xT = atoi( argv[2] );
double yT = atoi( argv[3] );
double zT = atoi( argv[4] );
int outfilenumber = atoi( argv[5] );

  std::ostringstream oss;
  oss << "Input_Fiber" << filenumber << ".txt"; // The fiber data is stored into a text file namely Fiber%.txt where % corresponds to 'filenumber'
  std::string filename = oss.str();
  std::ifstream fin(filename.c_str());
 
  std::string line;
  vnl_vector<double> xg(201); // 201 is the number of points for the input fiber (201='nb_points')
  vnl_vector<double> yg(201);
  vnl_vector<double> zg(201);
  int nb_points = 0;
  double x0,y0,z0;
  
  while(std::getline(fin, line))
    {
      double x,y,z;
      std::stringstream linestream;
      linestream << line;
      linestream >> x >> y >> z;
  
    if (nb_points == 0)
    {
      xg(nb_points) = x + xT; 
      yg(nb_points) = y + yT;
      zg(nb_points) = z + zT;
      x0 = x; 
      y0 = y;
      z0 = z;
    }
    else
    {
      xg(nb_points) = (x-x0) + xg(nb_points-1);
      yg(nb_points) = (y-y0) + yg(nb_points-1);
      zg(nb_points) = (z-z0) + zg(nb_points-1);;
      x0 = x;
      y0 = y;
      z0 = z;
    }
  
    nb_points++; 
    }
  fin.close();

    std::ofstream outfile;
    std::ostringstream oss2;
    oss2 << "Fiber" << outfilenumber << ".txt"; // save the output data
    std::string filename2 = oss2.str();
 
    outfile.open(filename2.c_str()); 
     
    for(int i=0; i < nb_points; i++)
    { 
      double randnum = 0.1 + (double)rand()/((double)RAND_MAX/(0.5-0.1));
      outfile << xg[i]+ randnum << " " << yg[i]+ randnum << " " << zg[i]+ randnum << std::endl;
    }
      outfile.close();


   return EXIT_SUCCESS;
}
