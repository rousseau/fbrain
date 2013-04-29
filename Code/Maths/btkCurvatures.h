/*==========================================================================

  © Université de Strasbourg - Centre National de la Recherche Scientifique

  Date: 08/04/2013
  Author(s): Aïcha Bentaieb (abentaieb@unistra.fr)

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

/* This algorithm computes different types of curvatures on triangular mesh polydata. The curvature
 * is set to be Gaussian, Barycentric, Mean or Principal direction curvatures.
 * The curvature tensor is computed based on Taubin Curvature tensor ( ESTIMATING THE TENSOR OF
 * CURVATURE OF A SURFACE FROM A POLYHEDRAL APPROXIMATION -Gabriel Taubin)
 **/


#ifndef BTKCURVATURES_H
#define BTKCURVATURES_H

#include "vtkPolyDataAlgorithm.h"

#define CURVATURE_TENSOR 0
#define CURVATURE_GAUSS 1
#define CURVATURE_MEAN 2
#define CURVATURE_BAR 3

class VTK_GRAPHICS_EXPORT btkCurvatures : public vtkPolyDataAlgorithm
{
public:
    vtkTypeRevisionMacro(btkCurvatures,vtkPolyDataAlgorithm);
    void PrintSelf(ostream& os, vtkIndent indent);

    static btkCurvatures *New();

    vtkSetMacro(CurvatureType,int);
    vtkGetMacro(CurvatureType,int);
    void SetCurvatureTypeToTensor()
    {
        this->SetCurvatureType(CURVATURE_TENSOR);
    }
    void SetCurvatureTypeToGaussian()
    {
        this->SetCurvatureType(CURVATURE_GAUSS);
    }
    void SetCurvatureTypeToMean()
    {
        this->SetCurvatureType(CURVATURE_MEAN);
    }
    void SetCurvatureTypeToBar()
    {
        this->SetCurvatureType(CURVATURE_BAR);
    }

    // Description:
    // Tolerance for floating point check as we compute the Jacobi of a diagonalized matrix
    vtkSetMacro(Tolerance,double);
    vtkGetMacro(Tolerance,double);


private:
    btkCurvatures();

    // Usual data generation method
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

    // main computation
    void GetCurvatureTensor(vtkPolyData *input);
    void GetGaussCurvature(vtkPolyData *input);
    void GetMeanCurvature(vtkPolyData *input);
    void GetBarCurvature(vtkPolyData *input);

    // variables
    int CurvatureType;
    double Tolerance;


};


#endif // BTKCURVATURES_H

