/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ShiftedGaugeWave_Startup(void)
{
  const char * banner = "ShiftedGaugeWave";
  CCTK_RegisterBanner(banner);
  return 0;
}
