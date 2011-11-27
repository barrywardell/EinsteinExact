/*  File produced by Kranc */

#include "cctk.h"

extern "C" int Vaidya2_Startup(void)
{
  const char * banner = "Vaidya2";
  CCTK_RegisterBanner(banner);
  return 0;
}
