/*  File produced by Kranc */

#include "cctk.h"

extern "C" int Vaidya2_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "Vaidya2";
  CCTK_RegisterBanner(banner);
  return 0;
}
