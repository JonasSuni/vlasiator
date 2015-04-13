/*
This file is part of Vlasiator.

Copyright 2012 Finnish Meteorological Institute












*/

#ifndef BACKGROUNDFIELD_H
#define BACKGROUNDFIELD_H

#include "fieldfunction.hpp"
#include "../definitions.h"
void setBackgroundField(
   FieldFunction& bgFunction,
   Real* cellParams,
   Real* faceDerivatives,
   Real* volumeDerivatives,
   bool append=false
);

void setBackgroundFieldToZero(
   Real* cellParams,
   Real* faceDerivatives,
   Real* volumeDerivatives
);

#endif

