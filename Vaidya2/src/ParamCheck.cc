#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
extern "C" void Vaidya2_ParamCheck(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_Vaidya2_ParamCheck
  DECLARE_CCTK_ARGUMENTS_CHECKED(Vaidya2_ParamCheck);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (CCTK_MyProc(cctkGH) == 0)
  {
    if (((CCTK_EQUALS(initial_data,"Vaidya")) || (CCTK_EQUALS(initial_lapse,"Vaidya")) || (CCTK_EQUALS(initial_shift,"Vaidya")) || (CCTK_EQUALS(initial_dtlapse,"Vaidya")) || (CCTK_EQUALS(initial_dtshift,"Vaidya"))) && ((!CCTK_EQUALS(initial_data,"Vaidya")) || (!CCTK_EQUALS(initial_lapse,"Vaidya")) || (!CCTK_EQUALS(initial_shift,"Vaidya")) || (!CCTK_EQUALS(initial_dtlapse,"Vaidya")) || (!CCTK_EQUALS(initial_dtshift,"Vaidya"))))
    {
      CCTK_ERROR("If one of the parameters ADMBase::initial_data, ADMBase::initial_lapse, ADMBase::initial_shift, ADMBase::initial_dtlapse, and ADMBase::initial_dtshift are set to \"Vaidya\", then all must be set to this value");
    }
  }
}
