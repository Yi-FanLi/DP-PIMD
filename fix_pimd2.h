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

FixStyle(pimd2,FixPIMD2)

#else

#ifndef FIX_PIMD2_H
#define FIX_PIMD2_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPIMD2 : public Fix {
 public:
  FixPIMD2(class LAMMPS *, int, char **);
//  virtual ~FixPIMD2();

  int setmask();

  void init();
  void setup(int);
  void post_force(int);
  void initial_integrate(int);
  void final_integrate();
  void end_of_step();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int,int,int);
  int pack_exchange(int,double*);
  int unpack_exchange(int,double*);
  int pack_restart(int,double*);
  void unpack_restart(int,int);
  int maxsize_restart();
  int size_restart(int);
  double compute_vector(int);

  int pack_forward_comm(int, int*, double *, int, int*);
  void unpack_forward_comm(int, int, double *);

  int method;
  int np;
  double inverse_np;
  double temp;
  int thermostat;

  /* ring-polymer model */

  double omega_np, fbond, spring_energy, sp;
  int x_last, x_next;

  void spring_force();

  /* fictitious mass */

  double fmass, *mass;

  /* inter-partition communication */

  int max_nsend;
  tagint* tag_send;
  double *buf_send;

  int max_nlocal;
  double *buf_recv, **buf_beads;

  int size_plan;
  int *plan_send, *plan_recv;
  double **comm_ptr;

  void comm_init();
  void comm_exec(double **);

  double **coords;
  int nsend, nrecv;
  tagint* tags_send;
  double *coords_send, *coords_recv;
  
  void comm_coords();

  /* centroid-virial kinetic energy estimator computation */
  double *xc;
  void compute_xc();
  double xf, t_vir;
  
  void compute_t_vir();

  /* primitive kinetic energy estimator computation */
  double total_spring_energy;
  double t_prim;

  void compute_t_prim();

  /* normal-mode operations */

  double *lam, **M_x2xp, **M_xp2x, **M_f2fp, **M_fp2f;
  int *mode_index;

  void nmpimd_init();
  void nmpimd_fill(double**);
  void nmpimd_transform(double**, double**, double*);

  /* Nose-hoover chain integration */

  int nhc_offset_one_1, nhc_offset_one_2;
  int nhc_size_one_1, nhc_size_one_2;
  int nhc_nchain;
  bool nhc_ready;
  double nhc_temp, dtv, dtf, t_sys;

  double **nhc_eta;        /* coordinates of NH chains for ring-polymer beads */
  double **nhc_eta_dot;    /* velocities of NH chains                         */
  double **nhc_eta_dotdot; /* acceleration of NH chains                       */
  double **nhc_eta_mass;   /* mass of NH chains                               */

  void nhc_init();
  void nhc_update_v();
  void nhc_update_x();

  /* Langevin thermostat BAOAB integration */

  bool baoab_ready;
  double gamma, c1, c2;
  double baoab_temp;

  class RanMars *random;
  int seed;

  void baoab_init();
  void baoab_update_v();
  void baoab_update_x();
  void random_v();
 
  double r1, r2, r3;
  
  /* harmonic oscillator model system */
  int harmonicflag;
  double omega;

  /* potential energy and total energy of the extended system */
  double hope, tote, totke;
  
  void compute_hope();
  void compute_tote();

  /* thermodynamic integration */
  int tiflag;
  int timethod;
  double lambda, dfdl;
  double **x_scaled;
  
  void compute_xscaled();
  void compute_dfdl();

};


}

#endif
#endif
