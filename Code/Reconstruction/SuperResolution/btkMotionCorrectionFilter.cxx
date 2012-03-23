#include "btkMotionCorrectionFilter.h"

namespace btk
{

MotionCorrectionFilter::MotionCorrectionFilter()
{
    btkCoutMacro(MotionCorrectionFilter : Constructor );
}

MotionCorrectionFilter::~MotionCorrectionFilter()
{
    btkCoutMacro(MotionCorrectionFilter : Destructor );
}

void MotionCorrectionFilter::Update()
{
     btkCoutMacro(MotionCorrectionFilter : Update Method );
}

}
