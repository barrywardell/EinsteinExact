#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
extern "C" void KerrSchild_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS_KerrSchild_ParamCheck;
  DECLARE_CCTK_PARAMETERS;
  
  if (CCTK_MyProc(cctkGH) == 0)
  {
    if (((CCTK_EQUALS(initial_data,"KerrSchild")) || (CCTK_EQUALS(initial_lapse,"KerrSchild")) || (CCTK_EQUALS(initial_shift,"KerrSchild")) || (CCTK_EQUALS(initial_dtlapse,"KerrSchild")) || (CCTK_EQUALS(initial_dtshift,"KerrSchild"))) && ((!CCTK_EQUALS(initial_data,"KerrSchild")) || (!CCTK_EQUALS(initial_lapse,"KerrSchild")) || (!CCTK_EQUALS(initial_shift,"KerrSchild")) || (!CCTK_EQUALS(initial_dtlapse,"KerrSchild")) || (!CCTK_EQUALS(initial_dtshift,"KerrSchild"))))
    {
      CCTK_ERROR("If one of the parameters ADMBase::initial_data, ADMBase::initial_lapse, ADMBase::initial_shift, ADMBase::initial_dtlapse, and ADMBase::initial_dtshift are set to \"KerrSchild\", then all must be set to this value");
    }
  }
}
