#include "btkHighResolutionIBPFilter.h"

namespace btk
{

HighResolutionIBPFilter::HighResolutionIBPFilter()
{
    btkCoutMacro(HighResolutionIBPFilter : Constructor );

    m_ImageHR = NULL;
}

HighResolutionIBPFilter::~HighResolutionIBPFilter()
{
    btkCoutMacro(HighResolutionIBPFilter : Destructor );
}

void HighResolutionIBPFilter::Update()
{
    btkCoutMacro(HighResolutionIBPFilter : Update Method );

}

}
