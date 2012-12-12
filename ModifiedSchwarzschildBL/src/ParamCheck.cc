#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ModifiedSchwarzschildBL_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (CCTK_MyProc(cctkGH) == 0)
  {
    
    if (((CCTK_EQUALS(initial_data,"ModifiedSchwarzschildBL") || CCTK_EQUALS(initial_lapse,"ModifiedSchwarzschildBL") || CCTK_EQUALS(initial_shift,"ModifiedSchwarzschildBL") || CCTK_EQUALS(initial_dtlapse,"ModifiedSchwarzschildBL") || CCTK_EQUALS(initial_dtshift,"ModifiedSchwarzschildBL")) && (!CCTK_EQUALS(initial_data,"ModifiedSchwarzschildBL") || !CCTK_EQUALS(initial_lapse,"ModifiedSchwarzschildBL") || !CCTK_EQUALS(initial_shift,"ModifiedSchwarzschildBL") || !CCTK_EQUALS(initial_dtlapse,"ModifiedSchwarzschildBL") || !CCTK_EQUALS(initial_dtshift,"ModifiedSchwarzschildBL"))))
    {
      CCTK_WARN(0, "The parameters ADMBase::initial_data ADMBase::initial_lapse, ADMBase::initial_shift, ADMBase::initial_dtlapse and ADMBase::initial_dtshift must all be set to \"ModifiedSchwarzschildBL\"");
    }
  }
}