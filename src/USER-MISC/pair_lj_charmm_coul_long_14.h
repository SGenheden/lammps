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

#ifdef PAIR_CLASS

PairStyle(lj/charmm/coul/long/14,PairLJCharmmCoulLong14)

#else

#ifndef LMP_PAIR_LJ_CHARMM_COUL_LONG_14_H
#define LMP_PAIR_LJ_CHARMM_COUL_LONG_14_H

#include "pair_lj_charmm_coul_long.h"

namespace LAMMPS_NS {

class PairLJCharmmCoulLong14 : public PairLJCharmmCoulLong {
 public:
  PairLJCharmmCoulLong14(class LAMMPS *);
  virtual ~PairLJCharmmCoulLong14();
  void* extract(const char *str, int &dim);
  void init_style();
  virtual void compute(int, int);
  void compute_inner();
  void compute_middle();
  virtual void compute_outer(int, int);
 protected :
   double** elescale;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style lj/charmm/coul/long requires atom attribute q

The atom style defined does not have these attributes.

E: Pair inner cutoff >= Pair outer cutoff

The specified cutoffs for the pair style are inconsistent.

E: Pair cutoff < Respa interior cutoff

One or more pairwise cutoffs are too short to use with the specified
rRESPA cutoffs.

E: Pair inner cutoff < Respa interior cutoff

One or more pairwise cutoffs are too short to use with the specified
rRESPA cutoffs.

E: Pair style requires a KSpace style

No kspace style is defined.

*/
