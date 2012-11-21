#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void Minkowski(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (CCTK_MyProc(cctkGH) == 0)
  {
    
    if (!CCTK_EQUALS(initial_data,"Minkowski") || (CCTK_EQUALS(initial_lapse,"Minkowski") && (CCTK_EQUALS(initial_shift,"Minkowski") && (CCTK_EQUALS(initial_dtlapse,"Minkowski") && CCTK_EQUALS(initial_dtshift,"Minkowski")))))
    {
      CCTK_WARN(0, "When ADMBase::initial_data = \"Minkowski\", the parameters ADMBase::initial_lapse, ADMBase::initial_shift, ADMBase::initial_dtlapse and ADMBase::initial_dtshift must also all be set to \"Minkowski\"");
    }
  }
}
