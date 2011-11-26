/*  File produced by Kranc */

#include "cctk.h"

extern "C" int Vaidya_Startup(void)
{
  const char * banner = "Vaidya";
  CCTK_RegisterBanner(banner);
  return 0;
}
