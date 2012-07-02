
#include "btkHighResolutionReconstructionFilter.h"

namespace btk
{
HighResolutionReconstructionFilter::HighResolutionReconstructionFilter()
{
    btkCoutMacro("HighResolutionReconstructionFilter : Constructor");
    //TODO: Default parameters
    m_InterpolationOrderPSF = 5;
    m_InterpolationOrderIBP = 1;

}
//-----------------------------------------------------------------------------------------------------------
HighResolutionReconstructionFilter::~HighResolutionReconstructionFilter()
{
    btkCoutMacro("HighResolutionReconstructionFilter : Destructor");
}
//-----------------------------------------------------------------------------------------------------------
}


