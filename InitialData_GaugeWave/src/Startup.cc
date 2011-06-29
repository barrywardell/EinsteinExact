/*  File produced by Kranc */

#include "cctk.h"

extern "C" int InitialData_GaugeWave_Startup(void)
{
  const char * banner = "InitialData_GaugeWave";
  CCTK_RegisterBanner(banner);
  return 0;
}
