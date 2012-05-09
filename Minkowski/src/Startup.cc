/*  File produced by Kranc */

#include "cctk.h"

extern "C" int Minkowski_Startup(void)
{
  const char * banner = "Minkowski";
  CCTK_RegisterBanner(banner);
  return 0;
}
