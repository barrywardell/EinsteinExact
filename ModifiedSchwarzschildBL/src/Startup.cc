/*  File produced by Kranc */

#include "cctk.h"

extern "C" int ModifiedSchwarzschildBL_Startup(void)
{
  const char* banner CCTK_ATTRIBUTE_UNUSED = "ModifiedSchwarzschildBL";
  CCTK_RegisterBanner(banner);
  return 0;
}
