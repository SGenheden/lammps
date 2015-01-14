/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "fix_temp_berendsen_sphere.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "group.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};
enum{CONSTANT,EQUAL};

/* ---------------------------------------------------------------------- */

FixTempBerendsenSphere::FixTempBerendsenSphere(LAMMPS *lmp, int narg, char **arg) :
  FixTempBerendsen(lmp, narg, arg)
{

  // create a new compute temp style
  // id = fix-ID + temp, compute group = fix group

  int n = strlen(id) + 9;
  id_temp_rot = new char[n];
  strcpy(id_temp_rot,id);
  strcat(id_temp_rot,"_temprot");

  char **newarg = new char*[5];
  newarg[0] = id_temp_rot;
  newarg[1] = group->names[igroup];
  newarg[2] = (char *) "temp/sphere";
  newarg[3] = (char *) "dof";
  newarg[4] = (char *) "rotate";
  modify->add_compute(5,newarg);
  delete [] newarg;

}

/* ---------------------------------------------------------------------- */

FixTempBerendsenSphere::~FixTempBerendsenSphere()
{

  delete [] tstr;

  // delete temperature if fix created it

  if (tflag) {
    modify->delete_compute(id_temp_rot); 
//    modify->delete_compute(id_temp);
  }
  //delete [] id_temp;
  delete [] id_temp_rot;
}

/* ---------------------------------------------------------------------- */

void FixTempBerendsenSphere::init()
{

  FixTempBerendsen::init();

  // check rotational temperature compute

  int icompute = modify->find_compute(id_temp_rot);
  if (icompute < 0)
    error->all(FLERR,"Temperature ID for fix temp/berendsen does not exist");
  temperature_rot = modify->compute[icompute];

}

/* ---------------------------------------------------------------------- */

void FixTempBerendsenSphere::end_of_step()
{
  double t_current = temperature->compute_scalar();
  double t_rot_current = temperature_rot->compute_scalar();
  if (t_current == 0.0 || t_rot_current == 0.0)
    error->all(FLERR,
               "Computed temperature for fix temp/berendsen cannot be 0.0");

  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  // set current t_target and t_rot_target
  // if variable temp, evaluate variable, wrap with clear/add

  if (tstyle == CONSTANT)
    t_target = t_start + delta * (t_stop-t_start);
  else {
    modify->clearstep_compute();
    t_target = input->variable->compute_equal(tvar);
    if (t_target < 0.0)
      error->one(FLERR,
                 "Fix temp/berendsen/sphere variable returned negative temperature");
    modify->addstep_compute(update->ntimestep + nevery);
  }

  // rescale velocities by lamda
  // for BIAS:
  //   temperature is current, so do not need to re-compute
  //   OK to not test returned v = 0, since lamda is multiplied by v

  double lamda = sqrt(1.0 + update->dt/t_period*(t_target/t_current - 1.0));
  double efactor = 0.5 * force->boltz * temperature->dof;
  energy += t_current * (1.0-lamda*lamda) * efactor;

  double lamda_rot = sqrt(1.0 + update->dt/t_period*(t_target/t_rot_current - 1.0));
  double efactor_rot = 0.5 * force->boltz * temperature_rot->dof;
  energy += t_rot_current * (1.0-lamda*lamda) * efactor_rot;

  double **v = atom->v;
  double **omega = atom->omega;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (which == NOBIAS) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        v[i][0] *= lamda;
        v[i][1] *= lamda;
        v[i][2] *= lamda;
        omega[i][0] *= lamda_rot;
        omega[i][1] *= lamda_rot;
        omega[i][2] *= lamda_rot;
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        temperature->remove_bias(i,v[i]);
        v[i][0] *= lamda;
        v[i][1] *= lamda;
        v[i][2] *= lamda;
        omega[i][0] *= lamda_rot;
        omega[i][1] *= lamda_rot;
        omega[i][2] *= lamda_rot;
        temperature->restore_bias(i,v[i]);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixTempBerendsenSphere::modify_param(int narg, char **arg)
{

  int parent_return = FixTempBerendsen::modify_param(narg,arg);

  if (strcmp(arg[0],"temp_rot") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    if (tflag) {
      modify->delete_compute(id_temp_rot);
      tflag = 0;
    }
    delete [] id_temp_rot;
    int n = strlen(arg[1]) + 1;
    id_temp_rot = new char[n];
    strcpy(id_temp_rot,arg[1]);

    int icompute = modify->find_compute(id_temp_rot);
    if (icompute < 0)
      error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature_rot = modify->compute[icompute];

    if (temperature_rot->tempflag == 0)
      error->all(FLERR,
                 "Fix_modify temperature ID does not compute temperature");
    if (temperature_rot->igroup != igroup && comm->me == 0)
      error->warning(FLERR,"Group for fix_modify temp != fix group");
    return 2;
  }
  return parent_return;
}

