/*  File produced by Kranc */

#include "cctk.h"

extern "C" int KerrSchild_Startup(void)
{
  const char * banner = "KerrSchild";
  CCTK_RegisterBanner(banner);
  return 0;
}
