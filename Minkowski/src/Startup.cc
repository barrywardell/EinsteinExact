/*  File produced by Kranc */

#include "cctk.h"

extern "C" int Minkowski_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "Minkowski";
  CCTK_RegisterBanner(banner);
  return 0;
}
