//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file coordinates.cpp
//! \brief

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "coordinates.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"

//----------------------------------------------------------------------------------------
// constructor, initializes coordinates data

Coordinates::Coordinates(ParameterInput *pin, MeshBlockPack *ppack) :
    pmy_pack(ppack) {
  // Check for relativistic dynamics
  is_special_relativistic = pin->GetOrAddBoolean("coord","special_rel",false);
}
