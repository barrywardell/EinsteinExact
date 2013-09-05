/*  File produced by Kranc */

#include "cctk.h"

extern "C" int GaugeWave_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "GaugeWave";
  CCTK_RegisterBanner(banner);
  return 0;
}
