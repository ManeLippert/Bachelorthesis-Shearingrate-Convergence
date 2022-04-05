/**********************************************
  Cross section tables header
**********************************************/

#ifndef __xsect_h
#define __xsect_h
#include "pdp1.h"

struct crossSectStruct
{
  int n;/* number of data pairs */
  SCALAR thresholdE; /* threshold energy */
  SCALAR* E; /* energy array */
  SCALAR* sigma; /* cross section array */
};

typedef struct crossSectStruct CrossSect;

int readXSectionTables(CrossSect** xSect, CrossSect* totalMom, char* fileName);
SCALAR interpolateTableLinear(CrossSect *xSect, SCALAR E);

#endif
