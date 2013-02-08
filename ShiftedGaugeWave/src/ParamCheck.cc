#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ShiftedGaugeWave_ParamCheck(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  
  if (CCTK_MyProc(cctkGH) == 0)
  {
    
    if (((CCTK_EQUALS(initial_data,"ShiftedGaugeWave") || CCTK_EQUALS(initial_lapse,"ShiftedGaugeWave") || CCTK_EQUALS(initial_shift,"ShiftedGaugeWave") || CCTK_EQUALS(initial_dtlapse,"ShiftedGaugeWave") || CCTK_EQUALS(initial_dtshift,"ShiftedGaugeWave")) && (!CCTK_EQUALS(initial_data,"ShiftedGaugeWave") || !CCTK_EQUALS(initial_lapse,"ShiftedGaugeWave") || !CCTK_EQUALS(initial_shift,"ShiftedGaugeWave") || !CCTK_EQUALS(initial_dtlapse,"ShiftedGaugeWave") || !CCTK_EQUALS(initial_dtshift,"ShiftedGaugeWave"))))
    {
      CCTK_WARN(0, "If one of the parameters ADMBase::initial_data, ADMBase::initial_lapse, ADMBase::initial_shift, ADMBase::initial_dtlapse, and ADMBase::initial_dtshift are set to \"ShiftedGaugeWave\", then all must be set to this value");
    }
  }
}
