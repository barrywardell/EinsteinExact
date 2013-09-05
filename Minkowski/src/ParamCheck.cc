#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void Minkowski_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (CCTK_MyProc(cctkGH) == 0)
  {
    
    if (((CCTK_EQUALS(initial_data,"Minkowski")) || (CCTK_EQUALS(initial_lapse,"Minkowski")) || (CCTK_EQUALS(initial_shift,"Minkowski")) || (CCTK_EQUALS(initial_dtlapse,"Minkowski")) || (CCTK_EQUALS(initial_dtshift,"Minkowski"))) && ((!CCTK_EQUALS(initial_data,"Minkowski")) || (!CCTK_EQUALS(initial_lapse,"Minkowski")) || (!CCTK_EQUALS(initial_shift,"Minkowski")) || (!CCTK_EQUALS(initial_dtlapse,"Minkowski")) || (!CCTK_EQUALS(initial_dtshift,"Minkowski"))))
    {
      CCTK_WARN(CCTK_WARN_ABORT, "If one of the parameters ADMBase::initial_data, ADMBase::initial_lapse, ADMBase::initial_shift, ADMBase::initial_dtlapse, and ADMBase::initial_dtshift are set to \"Minkowski\", then all must be set to this value");
    }
  }
}
