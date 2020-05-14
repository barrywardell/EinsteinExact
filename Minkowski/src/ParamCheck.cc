#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
extern "C" void Minkowski_ParamCheck(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_Minkowski_ParamCheck
  DECLARE_CCTK_ARGUMENTS_CHECKED(Minkowski_ParamCheck);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (CCTK_MyProc(cctkGH) == 0)
  {
    if (((CCTK_EQUALS(initial_data,"Minkowski")) || (CCTK_EQUALS(initial_lapse,"Minkowski")) || (CCTK_EQUALS(initial_shift,"Minkowski")) || (CCTK_EQUALS(initial_dtlapse,"Minkowski")) || (CCTK_EQUALS(initial_dtshift,"Minkowski"))) && ((!CCTK_EQUALS(initial_data,"Minkowski")) || (!CCTK_EQUALS(initial_lapse,"Minkowski")) || (!CCTK_EQUALS(initial_shift,"Minkowski")) || (!CCTK_EQUALS(initial_dtlapse,"Minkowski")) || (!CCTK_EQUALS(initial_dtshift,"Minkowski"))))
    {
      CCTK_ERROR("If one of the parameters ADMBase::initial_data, ADMBase::initial_lapse, ADMBase::initial_shift, ADMBase::initial_dtlapse, and ADMBase::initial_dtshift are set to \"Minkowski\", then all must be set to this value");
    }
  }
}
