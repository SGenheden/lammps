/***************************************************************************
                               base_ellipsoid.h
                             -------------------
                               W. Michael Brown

  Base class for acceleration of ellipsoid potentials

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : Thu May 5 2011
    email                : brownw@ornl.gov
 ***************************************************************************/

#ifndef BASE_ELLIPSOID_H
#define BASE_ELLIPSOID_H

#include "pair_gpu_device.h"
#include "pair_gpu_balance.h"
#include "mpi.h"

#ifdef USE_OPENCL
#include "geryon/ocl_texture.h"
#else
#include "geryon/nvd_texture.h"
#endif

namespace LAMMPS_AL {

template <class numtyp, class acctyp>
class BaseEllipsoid {
 public:
  BaseEllipsoid();
  virtual ~BaseEllipsoid();

  /// Clear any previous data and set up for a new LAMMPS run
  /** \param max_nbors initial number of rows in the neighbor matrix
    * \param cell_size cutoff + skin
    * \param gpu_split fraction of particles handled by device
    * 
    * Returns:
    * -  0 if successfull
    * - -1 if fix gpu not found
    * - -3 if there is an out of memory error
    * - -4 if the GPU library was not compiled for GPU
    * - -5 Double precision is not supported on card **/
  int init_base(const int nlocal, const int nall, const int max_nbors,
                const int maxspecial, const double cell_size,
                const double gpu_split, FILE *screen, const int ntypes,
                int **h_form, const char *nbor_program, 
                const char *ellipsoid_program, const char *lj_program);

  /// Estimate the overhead for GPU context changes and CPU driver
  void estimate_gpu_overhead();


  /// Check if there is enough storage for atom arrays and realloc if not
  /** \param success set to false if insufficient memory **/
  inline void resize_atom(const int nall, bool &success) {
    atom->resize(nall, success);
  }

  /// Check if there is enough storage for neighbors and realloc if not
  /** \param nlocal number of particles whose nbors must be stored on device
    * \param host_inum number of particles whose nbors need to copied to host
    * \param current maximum number of neighbors
    * \param olist_size size of list of particles from CPU neighboring
    * \note host_inum is 0 if the host is performing neighboring
    * \note if GPU is neighboring nlocal+host_inum=total number local particles
    * \note if CPU is neighboring olist_size=total number of local particles 
    * \note if GPU is neighboring olist_size=0 **/
  inline void resize_local(const int nlocal, const int host_inum,
                           const int max_nbors, const int olist_size,
                           bool &success) {
    ans->resize(nlocal, success);
    if (_multiple_forms) ans->dev_ans.zero();

    if (olist_size>static_cast<int>(host_olist.numel())) {
      host_olist.clear();
      int new_size=static_cast<int>(static_cast<double>(olist_size)*1.10);
      success=success && (host_olist.alloc(new_size,*ucl_device)==UCL_SUCCESS);
    }

    nbor->resize(nlocal,host_inum,max_nbors,success);
    double bytes=ans->gpu_bytes()+nbor->gpu_bytes();
    if (bytes>_max_bytes)
      _max_bytes=bytes;
  }

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear_base();
  
  /// Output any timing information
  void output_times();

  /// Returns memory usage on device per atom
  int bytes_per_atom(const int max_nbors) const;

  /// Total host memory used by library for pair style
  double host_memory_usage_base() const;

  /// Accumulate timers
  inline void acc_timers() {
    if (device->time_device()) {
      if (nbor_time_avail) {
        nbor->time_nbor.add_to_total();
        nbor->time_nbor.add_to_total();
        nbor_time_avail=false;
      }
      time_nbor1.add_to_total();
      time_ellipsoid.add_to_total();
      if (_multiple_forms) {
        time_nbor2.add_to_total();
        time_ellipsoid2.add_to_total();
        time_lj.add_to_total();
      }
      atom->acc_timers();
      ans->acc_timers();
    }
  }
  
  /// Zero timers
  inline void zero_timers() {
    nbor_time_avail=false;
    time_nbor1.zero();
    time_ellipsoid.zero();
    if (_multiple_forms) {
      time_nbor2.zero();
      time_ellipsoid2.zero();
      time_lj.zero();
    }
    atom->zero_timers();
    ans->zero_timers();
  }

  /// Pack neighbors to limit thread divergence for lj-lj and ellipse 
  void pack_nbors(const int GX, const int BX, const int start, const int inum,
                  const int form_low, const int form_high, 
                  const bool shared_types, int ntypes);

  /// Copy neighbor list from host
  void reset_nbors(const int nall, const int inum, const int osize, int *ilist,
                   int *numj, int *type, int **firstneigh, bool &success);

  /// Build neighbor list on device
  void build_nbor_list(const int inum, const int host_inum,
                       const int nall, double **host_x, int *host_type,
                       double *sublo, double *subhi, int *tag, int **nspecial,
                       int **special, bool &success);

  /// Pair loop with host neighboring
  int* compute(const int f_ago, const int inum_full, const int nall,
               double **host_x, int *host_type, int *ilist, int *numj,
               int **firstneigh, const bool eflag, const bool vflag,
               const bool eatom, const bool vatom, int &host_start,
               const double cpu_time, bool &success, double **quat);

  /// Pair loop with device neighboring
  int** compute(const int ago, const int inum_full, const int nall,
                double **host_x, int *host_type, double *sublo,
                double *subhi, int *tag, int **nspecial,
                int **special, const bool eflag, const bool vflag, 
                const bool eatom, const bool vatom, int &host_start, 
                int **ilist, int **numj, const double cpu_time, bool &success,
                double **host_quat);

  /// Build neighbor list on accelerator
  void build_nbor_list(const int inum, const int host_inum, const int nall, 
                       double **host_x, int *host_type, double *sublo,
                       double *subhi, bool &success);
                       
  // -------------------------- DEVICE DATA ------------------------- 

  /// Device Properties and Atom and Neighbor storage
  PairGPUDevice<numtyp,acctyp> *device;

  /// Geryon device
  UCL_Device *ucl_device;

  /// Device Timers
  UCL_Timer time_nbor1, time_ellipsoid, time_nbor2, time_ellipsoid2, time_lj;

  /// Host device load balancer
  PairGPUBalance<numtyp,acctyp> hd_balancer;

  /// LAMMPS pointer for screen output
  FILE *screen;

  // --------------------------- ATOM DATA --------------------------

  /// Atom Data
  PairGPUAtom<numtyp,acctyp> *atom;

  // --------------------------- TYPE DATA -------------------------- 

  /// cut_form.x = cutsq, cut_form.y = form
  UCL_D_Vec<numtyp2> cut_form;

  // ------------------------ FORCE/ENERGY DATA -----------------------

  PairGPUAns<numtyp,acctyp> *ans;

  // --------------------------- NBOR DATA ----------------------------

  /// Neighbor data
  PairGPUNbor *nbor;
  /// ilist with particles sorted by type
  UCL_H_Vec<int> host_olist;
  /// True if we need to accumulate time for neighboring
  bool nbor_time_avail;

  // ------------------------- DEVICE KERNELS -------------------------
  UCL_Program *nbor_program, *ellipsoid_program, *lj_program;
  UCL_Kernel k_nbor_fast, k_nbor;
  UCL_Kernel k_ellipsoid, k_sphere_ellipsoid, k_lj_fast, k_lj;
  inline int block_size() { return _block_size; }

  // --------------------------- TEXTURES -----------------------------
  UCL_Texture pos_tex;
  UCL_Texture q_tex;

 protected:
  bool _compiled;
  int _block_size, _threads_per_atom;
  double  _max_bytes, _max_an_bytes;
  double _gpu_overhead, _driver_overhead;
  UCL_D_Vec<int> *_nbor_data;

  // True if we want to use fast GB-sphere or sphere-sphere calculations 
  bool _multiple_forms;
  int **_host_form;
  int _last_ellipse, _max_last_ellipse;

  void compile_kernels(UCL_Device &dev, const char *nbor_string,
                       const char *ellipsoid_string, const char *lj_string);

  virtual void loop(const bool _eflag, const bool _vflag) = 0;
};

}

#endif

