/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void ModifiedSchwarzschildBL_RegisterVars(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_ModifiedSchwarzschildBL_RegisterVars
  DECLARE_CCTK_ARGUMENTS_CHECKED(ModifiedSchwarzschildBL_RegisterVars);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  /* Register all the evolved Array functions with MoL */
  return;
}
