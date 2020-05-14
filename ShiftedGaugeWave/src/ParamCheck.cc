#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
extern "C" void ShiftedGaugeWave_ParamCheck(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ShiftedGaugeWave_ParamCheck
  DECLARE_CCTK_ARGUMENTS_CHECKED(ShiftedGaugeWave_ParamCheck);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (CCTK_MyProc(cctkGH) == 0)
  {
    if (((CCTK_EQUALS(initial_data,"ShiftedGaugeWave")) || (CCTK_EQUALS(initial_lapse,"ShiftedGaugeWave")) || (CCTK_EQUALS(initial_shift,"ShiftedGaugeWave")) || (CCTK_EQUALS(initial_dtlapse,"ShiftedGaugeWave")) || (CCTK_EQUALS(initial_dtshift,"ShiftedGaugeWave"))) && ((!CCTK_EQUALS(initial_data,"ShiftedGaugeWave")) || (!CCTK_EQUALS(initial_lapse,"ShiftedGaugeWave")) || (!CCTK_EQUALS(initial_shift,"ShiftedGaugeWave")) || (!CCTK_EQUALS(initial_dtlapse,"ShiftedGaugeWave")) || (!CCTK_EQUALS(initial_dtshift,"ShiftedGaugeWave"))))
    {
      CCTK_ERROR("If one of the parameters ADMBase::initial_data, ADMBase::initial_lapse, ADMBase::initial_shift, ADMBase::initial_dtlapse, and ADMBase::initial_dtshift are set to \"ShiftedGaugeWave\", then all must be set to this value");
    }
  }
}
