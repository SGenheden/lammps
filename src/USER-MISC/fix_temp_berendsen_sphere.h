/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(temp/berendsen/sphere,FixTempBerendsenSphere)

#else

#ifndef LMP_FIX_TEMP_BERENDSEN_SPHERE_H
#define LMP_FIX_TEMP_BERENDSEN_SPHERE_H

#include "fix_temp_berendsen.h"

namespace LAMMPS_NS {

class FixTempBerendsenSphere : public FixTempBerendsen {
 public:
  FixTempBerendsenSphere(class LAMMPS *, int, char **);
  ~FixTempBerendsenSphere();
  void init();
  void end_of_step();
  int modify_param(int, char **);
 private:

  char *id_temp_rot;
  class Compute *temperature_rot;
};

}

#endif
#endif

