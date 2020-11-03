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

/* ----------------------------------------------------------------------
   Package      FixPIMD2
   Purpose      Quantum Path Integral Algorithm for Quantum Chemistry
   Copyright    Voth Group @ University of Chicago
   Authors      Chris Knight & Yuxing Peng (yuxing at uchicago.edu)
                Yifan Li (mail_liyifan@163.com)

   Updated      Nov-02-2020
   Version      2.0
------------------------------------------------------------------------- */

#include "fix_pimd2.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include "universe.h"
#include "comm.h"
#include "force.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "math_const.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{PIMD,NMPIMD,CMD};
enum{nhc, baoab};
enum{MSTI, SCTI};

/* ---------------------------------------------------------------------- */

FixPIMD2::FixPIMD2(LAMMPS *lmp, int narg, char **arg) : 
  Fix(lmp, narg, arg),
  random(NULL)
{
  method       = PIMD;
  thermostat   = nhc;
  fmass        = 1.0;
  temp         = 298.15;
  nhc_temp     = 298.15;
  nhc_nchain   = 2;
  baoab_temp   = 298.15;
  sp           = 1.0;
  harmonicflag = 0;
  omega        = 0.0;
  tiflag       = 0;
  timethod     = MSTI;
  lambda       = 0.0;

  for(int i=3; i<narg-1; i+=2)
  {
    if(strcmp(arg[i],"method")==0)
    {
      if(strcmp(arg[i+1],"pimd")==0) method=PIMD;
      else if(strcmp(arg[i+1],"nmpimd")==0) method=NMPIMD;
      else if(strcmp(arg[i+1],"cmd")==0) method=CMD;
      else error->universe_all(FLERR,"Unknown method parameter for fix pimd");
    }
    else if(strcmp(arg[i],"fmass")==0)
    {
      fmass = atof(arg[i+1]);
      if(fmass<0.0 || fmass>1.0) error->universe_all(FLERR,"Invalid fmass value for fix pimd");
    }
    else if(strcmp(arg[i],"sp")==0)
    {
      sp = atof(arg[i+1]);
      if(fmass<0.0) error->universe_all(FLERR,"Invalid sp value for fix pimd");
    }
    else if(strcmp(arg[i],"temp")==0)
    {
      temp = atof(arg[i+1]);
      if(temp<0.0) error->universe_all(FLERR,"Invalid temp value for fix pimd");
    } 
    else if(strcmp(arg[i], "thermostat")==0)
    {
      if(strcmp(arg[i+1],"nhc")==0)
      {
         thermostat = nhc;
      //   nhc_temp = temp;
         nhc_nchain = atoi(arg[i+2]);
         i++;
         if(nhc_nchain<2) error->universe_all(FLERR,"Invalid nhc value for fix pimd");
    
      }
      else if(strcmp(arg[i+1],"baoab")==0) 
      {
        thermostat = baoab;
       //  baoab_temp = temp;
        seed = atoi(arg[i+2]);
        i++;
      }
      else error->universe_all(FLERR,"Unknown thermostat parameter for fix pimd");
    }
    else if(strcmp(arg[i], "ti")==0)
    {
      tiflag = 1;
      if(strcmp(arg[i+1], "MSTI")==0)  timethod = MSTI;
      else if(strcmp(arg[i+1], "SCTI")==0)  timethod = SCTI;
      else error->universe_all(FLERR, "Unknown method parameter for thermodynamic integration");
      lambda = atof(arg[i+2]);
      i++;
    }
 
//    else if(strcmp(arg[i],"nhc")==0)
//    {
//      nhc_nchain = atoi(arg[i+1]);
//      if(nhc_nchain<2) error->universe_all(FLERR,"Invalid nhc value for fix pimd");
//    }
    else if(strcmp(arg[i], "model")==0)
    {
      harmonicflag = 1;
      omega = atof(arg[i+1]);
      if(omega<0) error->universe_all(FLERR,"Invalid model frequency value for fix pimd");
    }
    else error->universe_all(arg[i],i+1,"Unknown keyword for fix pimd");
  }

  // initialize Marsaglia RNG with processor-unique seed

  if(thermostat==baoab)
  {
    baoab_temp = temp;
    random = new RanMars(lmp, seed + universe->me);
  }
  
  if(thermostat==nhc)
  {
    nhc_temp = temp;
  }
  /* Initiation */

  max_nsend = 0;
  tag_send = NULL;
  buf_send = NULL;

  max_nlocal = 0;
  buf_recv = NULL;
  buf_beads = NULL;

  coords_send = coords_recv = NULL;
  nsend = nrecv = 0;
  tags_send = NULL;
  coords = NULL;
  size_plan = 0;
  plan_send = plan_recv = NULL;

  xc = NULL;
  xf = 0.0;
  t_vir = 0.0;
  total_spring_energy = 0.0;
  t_prim = 0.0;

  hope = 0.0;
  tote = 0.0;
  totke = 0.0;

  dfdl = 0.0;
  x_scaled = NULL;

  M_x2xp = M_xp2x = M_f2fp = M_fp2f = NULL;
  lam = NULL;
  mode_index = NULL;

  mass = NULL;

  array_atom = NULL;
  nhc_eta = NULL;
  nhc_eta_dot = NULL;
  nhc_eta_dotdot = NULL;
  nhc_eta_mass = NULL;

  size_peratom_cols = 12 * nhc_nchain + 3;

  nhc_offset_one_1 = 3 * nhc_nchain;
  nhc_offset_one_2 = 3 * nhc_nchain +3;
  nhc_size_one_1 = sizeof(double) * nhc_offset_one_1;
  nhc_size_one_2 = sizeof(double) * nhc_offset_one_2;

  gamma = 0.0;
  c1 = 0.0;
  c2 = 0.0;

  restart_peratom = 1;
  peratom_flag    = 1;
  peratom_freq    = 1;

  global_freq = 1;
  thermo_energy = 1;
  vector_flag = 1;
  size_vector = 7;
  extvector   = 1;
  comm_forward = 3;

  atom->add_callback(0); // Call LAMMPS to allocate memory for per-atom array
  atom->add_callback(1); // Call LAMMPS to re-assign restart-data for per-atom array

  grow_arrays(atom->nmax);

  // some initilizations

  baoab_ready = false;
  nhc_ready = false;

  r1 = 0.0;
  r2 = 0.0;
}

/* ---------------------------------------------------------------------- */
/*
FixPIMD2::~FixPIMD2()
{
  if(thermostat==baoab)
  {
    delete random;
  }
}
*/
/* ---------------------------------------------------------------------- */

int FixPIMD2::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::end_of_step()
{
  compute_hope();
  compute_tote();
  if(update->ntimestep % 10000 == 0)
  {
  if(universe->me==0) fprintf(stdout, "This is the end of step %d.\n", update->ntimestep);
  }
}
/* ---------------------------------------------------------------------- */

void FixPIMD2::init()
{
  if (atom->map_style == 0)
    error->all(FLERR,"Fix pimd requires an atom map, see atom_modify");

  if(universe->me==0 && screen) fprintf(screen,"Fix pimd initializing Path-Integral ...\n");
  // fprintf(stdout, "Fix pimd initilizing Path-Integral ...\n");

  // prepare the constants

  np = universe->nworlds;
  inverse_np = 1.0 / np;

  /* The first solution for the force constant, using SI units

  const double Boltzmann = 1.3806488E-23;    // SI unit: J/K
  const double Plank     = 6.6260755E-34;    // SI unit: m^2 kg / s

  double hbar = Plank / ( 2.0 * MY_PI ) * sp;
  double beta = 1.0 / ( Boltzmann * input.nh_temp);

  // - P / ( beta^2 * hbar^2)   SI unit: s^-2
  double _fbond = -1.0 / (beta*beta*hbar*hbar) * input.nbeads;

  // convert the units: s^-2 -> (kcal/mol) / (g/mol) / (A^2)
  fbond = _fbond * 4.184E+26;

  */

  /* The current solution, using LAMMPS internal real units */

  const double Boltzmann = force->boltz;
  const double Plank     = force->hplanck;

  double hbar   = Plank / ( 2.0 * MY_PI );
  double beta   = 1.0 / (Boltzmann * temp);
  double _fbond = 1.0 * np / (beta*beta*hbar*hbar) ;

  omega_np = sqrt(np) / (hbar * beta) * sqrt(force->mvv2e);
  fbond = - _fbond * force->mvv2e;

  if(universe->me==0)
    printf("Fix pimd -P/(beta^2 * hbar^2) = %20.7lE (kcal/mol/A^2)\n\n", fbond);

  if(thermostat==nhc) dtv = update->dt;
  else if(thermostat==baoab) dtv = 0.5 * update->dt;
  else
  {
    error->universe_all(FLERR,"Unknown thermostat parameter for fix pimd");
  }
  dtf = 0.5 * update->dt * force->ftm2v;

  comm_init();

  mass = new double [atom->ntypes+1];

  if(method==CMD || method==NMPIMD) nmpimd_init();
  else for(int i=1; i<=atom->ntypes; i++) mass[i] = atom->mass[i] / np * fmass;

  if(thermostat==baoab)
  {
    if(!baoab_ready)
    {
      baoab_init();
    }
    // fprintf(stdout, "baoab thermostat initialized!\n");
  }
  else if(thermostat==nhc)
  {
    if(!nhc_ready)
    {
      nhc_init();
    }
    // fprintf(stdout, "nhc thermostat initialized!\n");
  }
  else error->universe_all(FLERR,"Unknown thermostat parameter for fix pimd");

  // if(universe->me==0) fprintf(screen, "Fix pimd successfully initialized!\n");
  // fprintf(stdout, "Fix pimd successfully initialized!\n");
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::setup(int vflag)
{
  if(universe->me==0 && screen) fprintf(screen,"Setting up Path-Integral ...\n");
  // fprintf(stdout, "just before post_force()\n");
  // fprintf(stdout, "temperature: temp = %2.6f, baoab_temp = %2.6f.\n", temp, baoab_temp);
  post_force(vflag);
  // fprintf(stdout, "PIMD successfully set up!\n");
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::initial_integrate(int /*vflag*/)
{
  if(thermostat==baoab)
  {
    baoab_update_v();
    baoab_update_x();
    random_v();
    baoab_update_x();
  }
  else if(thermostat==nhc)
  {
    nhc_update_v();
    nhc_update_x();
  }
  else
  {
    error->universe_all(FLERR,"Unknown thermostat parameter for fix pimd");
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::final_integrate()
{
  if(thermostat==baoab)
  {
    baoab_update_v();
  }
  else if(thermostat==nhc)
  {
    nhc_update_v();
  }
  else
  {
    error->universe_all(FLERR,"Unknown thermostat parameter for fix pimd");
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::post_force(int /*flag*/)
{
  // fprintf(stdout, "Coming into post_force!\n");
  if(harmonicflag)
  {
    double* _mass = atom->mass;
    int *type = atom->type;
    for(int i=0; i<atom->nlocal; i++)
    {
      for(int j=0; j<3; j++)
      {
        atom->f[i][j] = -_mass[type[i]] * omega * omega * atom->x[i][j] * force->mvv2e;
      }
    }
  }
  // fprintf(stdout, "Successfully assigned harmonic forces!\n");
  for(int i=0; i<atom->nlocal; i++) 
  {
    for(int j=0; j<3; j++) 
    {
      atom->f[i][j] /= np;
    }
  }
  // fprintf(stdout, "me = %d, atom->x: %2.6f, %2.6f, %2.6f\n", universe->me, atom->x[0][0], atom->x[0][1], atom->x[0][2]);
  // fprintf(stdout, "me = %d, before ti, atom->f: %2.6f, %2.6f, %2.6f\n", universe->me, atom->f[0][0], atom->f[0][1], atom->f[0][2]);
  compute_t_vir();
  // fprintf(stdout, "me = %d, mass = %2.2f, omega = %2.6f, mvv2e = %2.6f, t_vir = %2.6e.\n", universe->me, atom->mass[1], omega, force->mvv2e, t_vir);
  // fprintf(stdout, "Force /np is OK!^_^\n");
  if(tiflag)
  {
    comm_coords();
    compute_xc();
    if(harmonicflag)
    {
      double *_mass = atom->mass;
      int *type = atom->type;
      if(timethod==MSTI)
      {
        for(int i=0; i<atom->nlocal; i++)
        {
          for(int j=0; j<3; j++)
          {
            atom->f[i][j] *= lambda;
            atom->f[i][j] -= (1.0 - lambda) / np * _mass[type[i]] * omega * omega * xc[3*i+j] * force->mvv2e;
          }
        }
      }
      else if(timethod==SCTI)
      {
        compute_xscaled();
        for(int i=0; i<atom->nlocal; i++)
        {
          for(int j=0; j<2; j++)
          {
            atom->f[i][j] = -1.0 * lambda / np * _mass[type[i]] * omega * omega * x_scaled[universe->iworld][3*i+j] * force->mvv2e;
            for(int k=0; k<np; k++)
            {
              atom->f[i][j] -= (1.0 - lambda) / np / np * _mass[type[i]] * omega * omega * x_scaled[k][3*i+j] * force->mvv2e;
            }
          }
        }
      }
    }
    else
    {
      error->universe_all(FLERR, "Only harmonic oscillator model system is supported in thermodynamic integration at this moment");
    }
    // fprintf(stdout, "me = %d, after ti, atom->f: %2.6f, %2.6f, %2.6f\n", universe->me, atom->f[0][0], atom->f[0][1], atom->f[0][2]);

    compute_dfdl();
  }
  // comm_coords();
  // fprintf(stdout, "Communicating coordinates has been finished!\n");
  // fprintf(stdout, "Successfully computed t_vir!\n");
  comm_exec(atom->x);
  // fprintf(stdout, "Successfully communicated x!\n");
  spring_force();
  // fprintf(stdout, "Successfully executed spring_force()!\n");
  MPI_Allreduce(&spring_energy, &total_spring_energy, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  total_spring_energy *= 0.25;
  compute_t_prim();
  // fprintf(stdout, "Successfully computed t_prim!\n");
  if(method==CMD || method==NMPIMD)
  {
    /* forward comm for the force on ghost atoms */

    nmpimd_fill(atom->f);

    /* inter-partition comm */

    comm_exec(atom->f);

    /* normal-mode transform */

    nmpimd_transform(buf_beads, atom->f, M_f2fp[universe->iworld]);
  }
  // fprintf(stdout, "This is the end of post_force()!\n");
  // if(universe->me == 0) 
  // fprintf(stdout, "me = %d, coordinate at the end of pf: %2.6e, %2.6e, %2.6e.\n", universe->me, atom->x[0][0], atom->x[0][1], atom->x[0][2]);

}

/* ----------------------------------------------------------------------
   Nose-Hoover Chains
------------------------------------------------------------------------- */

void FixPIMD2::nhc_init()
{
  double tau = 1.0 / omega_np;
  double KT  = force->boltz * nhc_temp;

  double mass0 = KT * tau * tau;
  int max = 3 * atom->nlocal;

  for(int i=0; i<max; i++)
  {
    for(int ichain=0; ichain<nhc_nchain; ichain++)
    {
      nhc_eta[i][ichain]        = 0.0;
      nhc_eta_dot[i][ichain]    = 0.0;
      nhc_eta_dot[i][ichain]    = 0.0;
      nhc_eta_dotdot[i][ichain] = 0.0;
      nhc_eta_mass[i][ichain]   = mass0;
      if((method==CMD || method==NMPIMD) && universe->iworld==0) ; else nhc_eta_mass[i][ichain]  *= fmass;
    }

    nhc_eta_dot[i][nhc_nchain]    = 0.0;

    for(int ichain=1; ichain<nhc_nchain; ichain++)
      nhc_eta_dotdot[i][ichain] = (nhc_eta_mass[i][ichain-1] * nhc_eta_dot[i][ichain-1]
        * nhc_eta_dot[i][ichain-1] * force->mvv2e - KT) / nhc_eta_mass[i][ichain];
  }

  // Zero NH acceleration for CMD

  if(method==CMD && universe->iworld==0) for(int i=0; i<max; i++)
    for(int ichain=0; ichain<nhc_nchain; ichain++) nhc_eta_dotdot[i][ichain] = 0.0;

  nhc_ready = true;
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::nhc_update_x()
{
  int n = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;

  if(method==CMD || method==NMPIMD)
  {
    nmpimd_fill(atom->v);
    comm_exec(atom->v);

    /* borrow the space of atom->f to store v in cartisian */

    v = atom->f;
    nmpimd_transform(buf_beads, v, M_xp2x[universe->iworld]);
  }

  for(int i=0; i<n; i++)
  {
    x[i][0] += dtv * v[i][0];
    x[i][1] += dtv * v[i][1];
    x[i][2] += dtv * v[i][2];
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::nhc_update_v()
{
  int n = atom->nlocal;
  int *type = atom->type;
  double **v = atom->v;
  double **f = atom->f;

  for(int i=0; i<n; i++)
  {
    double dtfm = dtf / mass[type[i]];
    v[i][0] += dtfm * f[i][0];
    v[i][1] += dtfm * f[i][1];
    v[i][2] += dtfm * f[i][2];
  }

  t_sys = 0.0;
  if(method==CMD && universe->iworld==0) return;

  double expfac;
  int nmax = 3 * atom->nlocal;
  double KT = force->boltz * nhc_temp;
  double kecurrent, t_current;

  double dthalf = 0.5   * update->dt;
  double dt4    = 0.25  * update->dt;
  double dt8    = 0.125 * update->dt;

  for(int i=0; i<nmax; i++)
  {
    int iatm = i/3;
    int idim = i%3;

    double *vv = v[iatm];

    kecurrent = mass[type[iatm]] * vv[idim]* vv[idim] * force->mvv2e;
    t_current = kecurrent / force->boltz;

    double *eta = nhc_eta[i];
    double *eta_dot = nhc_eta_dot[i];
    double *eta_dotdot = nhc_eta_dotdot[i];

    eta_dotdot[0] = (kecurrent - KT) / nhc_eta_mass[i][0];

    for(int ichain=nhc_nchain-1; ichain>0; ichain--)
    {
      expfac = exp(-dt8 * eta_dot[ichain+1]);
      eta_dot[ichain] *= expfac;
      eta_dot[ichain] += eta_dotdot[ichain] * dt4;
      eta_dot[ichain] *= expfac;
    }

    expfac = exp(-dt8 * eta_dot[1]);
    eta_dot[0] *= expfac;
    eta_dot[0] += eta_dotdot[0] * dt4;
    eta_dot[0] *= expfac;

    // Update particle velocities half-step

    double factor_eta = exp(-dthalf * eta_dot[0]);
    vv[idim] *= factor_eta;

    t_current *= (factor_eta * factor_eta);
    kecurrent = force->boltz * t_current;
    eta_dotdot[0] = (kecurrent - KT) / nhc_eta_mass[i][0];

    for(int ichain=0; ichain<nhc_nchain; ichain++)
      eta[ichain] += dthalf * eta_dot[ichain];

    eta_dot[0] *= expfac;
    eta_dot[0] += eta_dotdot[0] * dt4;
    eta_dot[0] *= expfac;

    for(int ichain=1; ichain<nhc_nchain; ichain++)
    {
      expfac = exp(-dt8 * eta_dot[ichain+1]);
      eta_dot[ichain] *= expfac;
      eta_dotdot[ichain] = (nhc_eta_mass[i][ichain-1] * eta_dot[ichain-1] * eta_dot[ichain-1]
                           - KT) / nhc_eta_mass[i][ichain];
      eta_dot[ichain] += eta_dotdot[ichain] * dt4;
      eta_dot[ichain] *= expfac;
    }

    t_sys += t_current;
  }

  t_sys /= nmax;

}

/* ----------------------------------------------------------------------
   Langevin thermostat, BAOAB integrator
------------------------------------------------------------------------- */

void FixPIMD2::baoab_init()
{
  double KT = force->boltz * baoab_temp;
/*
  gamma = (double*) memory->srealloc(gamma, sizeof(double)*np*3, "FixPIMD2::gamma");
  c1 = (double*) memory->srealloc(c1, sizeof(double)*np*3, "FixPIMD2::c1");
  c2 = (double*) memory->srealloc(c2, sizeof(double)*np*3, "FixPIMD2::c2");
  for(int i=0; i<3*np; i++)
  {
    gamma[i] = omega_np;
    c1[i] = exp(-gamma[i] * update->dt);
    c2[i] = sqrt(1.0 - c1[i] * c1[i]);
  }
*/
  double beta = 1.0 / force->boltz / temp;
  double hbar = force->hplanck / (2.0 * MY_PI);
  gamma = sqrt(np) / beta / hbar;
  // gamma = omega_np;
  c1 = exp(-gamma * update->dt);
  c2 = sqrt(1.0 - c1 * c1);
  // fprintf(stdout, "gamma = %2.6e, c1 = %2.6e, c2 = %2.6e.\n", gamma, c1, c2);
  baoab_ready = true;
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::baoab_update_v()
{
  int n = atom->nlocal;
  int *type = atom->type;
  double **v = atom->v;
  double **f = atom->f;

  for(int i=0; i<n; i++)
  {
    double dtfm = dtf / mass[type[i]];
    v[i][0] += dtfm * f[i][0];
    v[i][1] += dtfm * f[i][1];
    v[i][2] += dtfm * f[i][2];
  }

}

/* ---------------------------------------------------------------------- */

void FixPIMD2::baoab_update_x(){
  int n = atom->nlocal;
  double **x = atom->x;
  double **v = atom->v;

  if(method==NMPIMD)
  {
    nmpimd_fill(atom->v);
    comm_exec(atom->v);
    v = atom->f;
    nmpimd_transform(buf_beads, v, M_xp2x[universe->iworld]);
  }

  for(int i=0; i<n; i++)
  {
    x[i][0] += dtv * v[i][0];
    x[i][1] += dtv * v[i][1];
    x[i][2] += dtv * v[i][2];
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::random_v()
{
  int nlocal = atom->nlocal;
  int *type = atom->type;
  double beta = 1.0 / force->boltz / baoab_temp * force->mvv2e;
   // if(universe->me == 0) fprintf(stdout, "r1 = %2.6f, r2 = %2.6f, r3 = %2.6f.\n", r1, r2, r3);
  // if(universe->me == 0) fprintf(stdout, "velocity before rv: %2.6e, %2.6e, %2.6e.\n", atom->v[0][0], atom->v[0][1], atom->v[0][2]);
  for(int i=0; i<nlocal; i++)
  {
   /* atom->v[i][0] = c1[3*universe->iworld] * atom->v[i][0] + c2[3*universe->iworld] * sqrt(1.0 / mass[type[i]] / beta) * random->gaussian(); 
    atom->v[i][1] = c1[3*universe->iworld+1] * atom->v[i][1] + c2[3*universe->iworld+1] * sqrt(1.0 / mass[type[i]] / beta) * random->gaussian();
    atom->v[i][2] = c1[3*universe->iworld+2] * atom->v[i][2] + c2[3*universe->iworld+2] * sqrt(1.0 / mass[type[i]] / beta) * random->gaussian();*/
    r1 = random->gaussian();
    r2 = random->gaussian();
    r3 = random->gaussian();
    atom->v[i][0] = c1 * atom->v[i][0] + c2 * sqrt(1.0 / mass[type[i]] / beta) * r1; 
    atom->v[i][1] = c1 * atom->v[i][1] + c2 * sqrt(1.0 / mass[type[i]] / beta) * r2;
    atom->v[i][2] = c1 * atom->v[i][2] + c2 * sqrt(1.0 / mass[type[i]] / beta) * r3;
  }
  // if(universe->me == 0) fprintf(stdout, "velocity after rv: %2.6e, %2.6e, %2.6e.\n", atom->v[0][0], atom->v[0][1], atom->v[0][2]); 
}

/* ----------------------------------------------------------------------
   Normal Mode PIMD
------------------------------------------------------------------------- */

void FixPIMD2::nmpimd_init()
{
  memory->create(M_x2xp, np, np, "fix_feynman:M_x2xp");
  memory->create(M_xp2x, np, np, "fix_feynman:M_xp2x");
  memory->create(M_f2fp, np, np, "fix_feynman:M_f2fp");
  memory->create(M_fp2f, np, np, "fix_feynman:M_fp2f");

  lam = (double*) memory->smalloc(sizeof(double)*np, "FixPIMD2::lam");

  // Set up  eigenvalues

  lam[0] = 0.0;
  if(np%2==0) lam[np-1] = 4.0 * np;

  for(int i=2; i<=np/2; i++)
  {
    lam[2*i-3] = lam[2*i-2] = 2.0 * np * (1.0 - 1.0 *cos(2.0*MY_PI*(i-1)/np));
  }

  // Set up eigenvectors for non-degenerated modes

  for(int i=0; i<np; i++)
  {
    M_x2xp[0][i] = 1.0 / np;
    if(np%2==0) M_x2xp[np-1][i] = 1.0 / np * pow(-1.0, i);
  }

  // Set up eigenvectors for degenerated modes

  for(int i=0; i<(np-1)/2; i++) for(int j=0; j<np; j++)
  {
    M_x2xp[2*i+1][j] =   sqrt(2.0) * cos ( 2.0 * MY_PI * (i+1) * j / np) / np;
    M_x2xp[2*i+2][j] = - sqrt(2.0) * sin ( 2.0 * MY_PI * (i+1) * j / np) / np;
  }

  // Set up Ut

  for(int i=0; i<np; i++)
    for(int j=0; j<np; j++)
    {
      M_xp2x[i][j] = M_x2xp[j][i] * np;
      M_f2fp[i][j] = M_x2xp[i][j] * np;
      M_fp2f[i][j] = M_xp2x[i][j];
    }

  // Set up masses

  int iworld = universe->iworld;

  for(int i=1; i<=atom->ntypes; i++)
  {
    mass[i] = atom->mass[i];

    if(iworld)
    {
      mass[i] *= lam[iworld];
      mass[i] *= fmass;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::nmpimd_fill(double **ptr)
{
  comm_ptr = ptr;
  comm->forward_comm_fix(this);
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::nmpimd_transform(double** src, double** des, double *vector)
{
  int n = atom->nlocal;
  int m = 0;

  for(int i=0; i<n; i++) for(int d=0; d<3; d++)
  {
    des[i][d] = 0.0;
    for(int j=0; j<np; j++) { des[i][d] += (src[j][m] * vector[j]); }
    m++;
  }
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::spring_force()
{
  spring_energy = 0.0;

  double **x = atom->x;
  double **f = atom->f;
  double* _mass = atom->mass;
  int* type = atom->type;
  int nlocal = atom->nlocal;

  double* xlast = buf_beads[x_last];
  double* xnext = buf_beads[x_next];

  for(int i=0; i<nlocal; i++)
  {
    double delx1 = xlast[0] - x[i][0];
    double dely1 = xlast[1] - x[i][1];
    double delz1 = xlast[2] - x[i][2];
    xlast += 3;
    domain->minimum_image(delx1, dely1, delz1);

    double delx2 = xnext[0] - x[i][0];
    double dely2 = xnext[1] - x[i][1];
    double delz2 = xnext[2] - x[i][2];
    xnext += 3;
    domain->minimum_image(delx2, dely2, delz2);

    double ff = fbond * _mass[type[i]];

    double dx = delx1+delx2;
    double dy = dely1+dely2;
    double dz = delz1+delz2;

    f[i][0] -= (dx) * ff;
    f[i][1] -= (dy) * ff;
    f[i][2] -= (dz) * ff;

    // spring_energy += (dx*dx+dy*dy+dz*dz);
    spring_energy += -ff * (delx1*delx1+dely1*dely1+delz1*delz1+delx2*delx2+dely2*dely2+delz2*delz2);
  }
}

/* ----------------------------------------------------------------------
   Comm operations
------------------------------------------------------------------------- */

void FixPIMD2::comm_init()
{
  if(size_plan)
  {
    delete [] plan_send;
    delete [] plan_recv;
  }

  if(method == PIMD)
  {
    size_plan = 2;
    plan_send = new int [2];
    plan_recv = new int [2];
    mode_index = new int [2];

    int rank_last = universe->me - comm->nprocs;
    int rank_next = universe->me + comm->nprocs;
    if(rank_last<0) rank_last += universe->nprocs;
    if(rank_next>=universe->nprocs) rank_next -= universe->nprocs;

    plan_send[0] = rank_next; plan_send[1] = rank_last;
    plan_recv[0] = rank_last; plan_recv[1] = rank_next;

    mode_index[0] = 0; mode_index[1] = 1;
    x_last = 1; x_next = 0;
  }
  else
  {
    size_plan = np - 1;
    plan_send = new int [size_plan];
    plan_recv = new int [size_plan];
    mode_index = new int [size_plan];

    for(int i=0; i<size_plan; i++)
    {
      plan_send[i] = universe->me + comm->nprocs * (i+1);
      if(plan_send[i]>=universe->nprocs) plan_send[i] -= universe->nprocs;

      plan_recv[i] = universe->me - comm->nprocs * (i+1);
      if(plan_recv[i]<0) plan_recv[i] += universe->nprocs;

      mode_index[i]=(universe->iworld+i+1)%(universe->nworlds);
    }

    x_next = (universe->iworld+1+universe->nworlds)%(universe->nworlds);
    x_last = (universe->iworld-1+universe->nworlds)%(universe->nworlds);
  }

  if(buf_beads)
  {
    for(int i=0; i<np; i++) if(buf_beads[i]) delete [] buf_beads[i];
    delete [] buf_beads;
  }

  buf_beads = new double* [np];
  for(int i=0; i<np; i++) buf_beads[i] = NULL;
  
  if(coords)
  {
    for(int i=0; i<np; i++) if(coords[i]) delete [] coords[i];
    delete [] coords;
  }
  
  coords = new double* [np];
  for(int i=0; i<np; i++) coords[i] = NULL;
  
  if(x_scaled)
  {
    for(int i=0; i<np; i++) if(x_scaled[i]) delete [] x_scaled[i];
    delete [] x_scaled;
  }

  x_scaled = new double* [np];
  for(int i=0; i<np; i++) x_scaled[i] = NULL;
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::comm_exec(double **ptr)
{
  int nlocal = atom->nlocal;

  if(nlocal > max_nlocal)
  {
    max_nlocal = nlocal+200;
    int size = sizeof(double) * max_nlocal * 3;
    buf_recv = (double*) memory->srealloc(buf_recv, size, "FixPIMD2:x_recv");

    for(int i=0; i<np; i++)
      buf_beads[i] = (double*) memory->srealloc(buf_beads[i], size, "FixPIMD2:x_beads[i]");
  }

  // copy local positions

  memcpy(buf_beads[universe->iworld], &(ptr[0][0]), sizeof(double)*nlocal*3);

  // go over comm plans

  for(int iplan = 0; iplan<size_plan; iplan++)
  {
    // sendrecv nlocal

    int nsend;

    MPI_Sendrecv( &(nlocal), 1, MPI_INT, plan_send[iplan], 0,
                  &(nsend),  1, MPI_INT, plan_recv[iplan], 0, universe->uworld, MPI_STATUS_IGNORE);

    // allocate arrays

    if(nsend > max_nsend)
    {
      max_nsend = nsend+200;
      tag_send = (tagint*) memory->srealloc(tag_send, sizeof(tagint)*max_nsend, "FixPIMD2:tag_send");
      buf_send = (double*) memory->srealloc(buf_send, sizeof(double)*max_nsend*3, "FixPIMD2:x_send");
    }

    // send tags

    MPI_Sendrecv( atom->tag, nlocal, MPI_LMP_TAGINT, plan_send[iplan], 0,
                  tag_send,  nsend,  MPI_LMP_TAGINT, plan_recv[iplan], 0, universe->uworld, MPI_STATUS_IGNORE);

    // wrap positions

    double *wrap_ptr = buf_send;
    int ncpy = sizeof(double)*3;

    for(int i=0; i<nsend; i++)
    {
      int index = atom->map(tag_send[i]);

      if(index<0)
      {
        char error_line[256];

        sprintf(error_line, "Atom " TAGINT_FORMAT " is missing at world [%d] "
                "rank [%d] required by  rank [%d] (" TAGINT_FORMAT ", "
                TAGINT_FORMAT ", " TAGINT_FORMAT ").\n", tag_send[i],
                universe->iworld, comm->me, plan_recv[iplan],
                atom->tag[0], atom->tag[1], atom->tag[2]);

        error->universe_one(FLERR,error_line);
      }

      memcpy(wrap_ptr, ptr[index], ncpy);
      wrap_ptr += 3;
    }

    // sendrecv x

    MPI_Sendrecv( buf_send, nsend*3,  MPI_DOUBLE, plan_recv[iplan], 0,
                  buf_recv, nlocal*3, MPI_DOUBLE, plan_send[iplan], 0, universe->uworld, MPI_STATUS_IGNORE);

    // copy x

    memcpy(buf_beads[mode_index[iplan]], buf_recv, sizeof(double)*nlocal*3);
  }
/*  fprintf(stdout, "me = %d, nlocal = %d, tags: ", universe->me, nlocal);
  for(int i=0; i<nlocal; i++) fprintf(stdout, "%d ", atom->tag[i]);
  fprintf(stdout, "\nme = %d, buf_beads: \n", universe->me);
  for(int i=0; i<np; i++)
  {
    fprintf(stdout, "me = %d, ", universe->me);
    for(int j=0; j<nlocal; j++)
    {
      fprintf(stdout, "%2.6f, %2.6f, %2.6f,     ", buf_beads[i][3*j], buf_beads[i][3*j+1], buf_beads[i][3*j+2]);
    }
    fprintf(stdout, "\n");
  } 
  fprintf(stdout, "\n");*/
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::comm_coords()
{
  int nlocal = atom->nlocal;

  // fprintf(stdout, "me = %d, nlocal = %d, tags: \n", universe->me, nlocal);
  // for(int i=0; i<nlocal; i++)
  // {
  //   fprintf(stdout, "%d ", atom->tag[i]);
  // }
  // fprintf(stdout, "\n");
  // fprintf(stdout, "atom->x: \n");
  // for(int i=0; i<nlocal; i++)
  // {
  //   fprintf(stdout, "%2.6e, %2.6e, %2.6e, \n", atom->x[i][0], atom->x[i][1], atom->x[i][2]);
  // }
  // fprintf(stdout, "\n\n");

  // assign memory for arrays
  int size_coords = sizeof(double) * nlocal * 3;
  int size_tags;// = sizeof(tagint) * nlocal;
  coords_recv = (double*) memory->srealloc(coords_recv, size_coords, "FixPIMD2:coords_recv");
  // tags_recv = (tagint*) memory->srealloc(tags_recv, size_tags, "FixPIMD2:tags_recv");
  for(int i=0; i<np; i++)
  {
    coords[i] = (double*) memory->srealloc(coords[i], size_coords, "FixPIMD2:coords[i]");
  }
  
  // copy local positions and tags
//  memcpy(coords_send, &(atom->x[0][0]), size_coords);
//  memcpy(tags_send, &(atom->tag[0]), size_tags);
  memcpy(coords[universe->iworld], &(atom->x[0][0]), size_coords);
  // fprintf(stdout, "me = %d, before traversing.\n", universe->me); 
  // traversing over all the other worlds
  for(int dworld=1; dworld<=np-1; dworld++)
  {
//    for(int iproc=0; iproc<comm->nprocs; iproc++)
      // send the tags and coords to the process proc_send
      // receive the tags and coords from the process proc_recv
//      int proc_send = ((universe->iworld + dworld) % universe->nworlds) * comm->nprocs + iproc;
//      int proc_recv = ((universe->iworld - dworld + universe->nworlds) % universe->nworlds) * comm->nprocs + iproc;
    int proc_send = (universe->me + dworld * comm->nprocs) % universe->nprocs; 
    int proc_recv = (universe->me - dworld * comm->nprocs + universe->nprocs) % universe->nprocs;
    int world_recv = (int)(proc_recv / comm->nprocs);
      // fprintf(stdout, "me = %d, proc_send = %d, proc_recv = %d, world_recv = %d.\n", universe->me, proc_send, proc_recv, world_recv);
    
    // determine the number of atoms to be sent to and received from the other worlds
    MPI_Sendrecv(&(nlocal), 1, MPI_INT, proc_send, 0, 
                 &(nsend), 1, MPI_INT, proc_recv, 0, 
                 universe->uworld, MPI_STATUS_IGNORE);
    nrecv = nlocal;

    size_coords = sizeof(double) * nsend * 3;
    size_tags = sizeof(tagint) * nsend;

    coords_send = (double*) memory->srealloc(coords_send, size_coords, "FixPIMD2:coords_send");
    tags_send = (tagint*) memory->srealloc(tags_send, size_tags, "FixPIMD2:tags_send");

    MPI_Sendrecv(atom->tag, nlocal, MPI_LMP_TAGINT, proc_send, 0,
                 tags_send, nsend, MPI_LMP_TAGINT, proc_recv, 0,
                 universe->uworld, MPI_STATUS_IGNORE);
    // fprintf(stdout, "me = %d, dowrld = %d, before warpping.\n", universe->me, dworld);
    // wrap positions
    double *wrap_ptr = coords_send;
    int ncpy = sizeof(double)*3;

    for(int i=0; i<nsend; i++)
    {
      int index = atom->map(tags_send[i]);
      if(index < 0)
      {
        char error_line[256];

        sprintf(error_line, "Atom " TAGINT_FORMAT " is missing at world [%d] "
                "rank [%d] required by  rank [%d] (" TAGINT_FORMAT ", "
                TAGINT_FORMAT ", " TAGINT_FORMAT ").\n", tags_send[i],
                universe->iworld, comm->me, proc_recv,
                atom->tag[0], atom->tag[1], atom->tag[2]);

        error->universe_one(FLERR,error_line);
      }    
      memcpy(wrap_ptr, atom->x[index], ncpy);
      wrap_ptr += 3;
    }
    // fprintf(stdout, "me = %d, dworld = %d, after wrapping, before Sendrecv.\n", universe->me, dworld); 
    MPI_Sendrecv(coords_send, nsend*3, MPI_DOUBLE, proc_recv, 0,
                coords_recv, nrecv*3, MPI_DOUBLE, proc_send, 0,
                universe->uworld, MPI_STATUS_IGNORE);
    // fprintf(stdout, "me = %d, dworld = %d, after Sendrecv.\n", universe->me, dworld);

    memcpy(coords[world_recv], coords_recv, sizeof(double)*nlocal*3);          
  }
  // fprintf(stdout, "me = %d, comm_coords completed!\n", universe->me);
  // test comm_coords()
/*  fprintf(stdout, "me = %d, tags:\n", universe->me);
  for(int i=0; i<nlocal; i++)
  {
    fprintf(stdout, "%d ", atom->tag[i]);
  }
  fprintf(stdout, "\n");
  fprintf(stdout, "me = %d, coords:\n", universe->me);
  for(int i=0; i<np; i++)
  {
    for(int j=0; j<nlocal; j++)
    {
      fprintf(stdout, "%2.6e, %2.6e, %2.6e,     ", coords[i][3*j], coords[i][3*j+1], coords[i][3*j+2]);
    }
    fprintf(stdout, "\n");
  }
*/
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::compute_xc()
{
  int nlocal = atom->nlocal;
  xc = (double*) memory->srealloc(xc, sizeof(double) * nlocal * 3, "FixPIMD2:xc");
  for(int i=0; i<nlocal; i++)
  {
    xc[3*i] = xc[3*i+1] = xc[3*i+2] = 0.0;
    for(int j=0; j<np; j++)
    {
      xc[3*i] += coords[j][3*i];
      xc[3*i+1] += coords[j][3*i+1];
      xc[3*i+2] += coords[j][3*i+2];
    }
    xc[3*i] /= np;
    xc[3*i+1] /= np;
    xc[3*i+2] /= np;
  } 
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::compute_xscaled()
{
  int nlocal = atom->nlocal;
  for(int i=0; i<np; i++)
  {
    x_scaled[i] = (double*) memory->srealloc(x_scaled[i], sizeof(double) * nlocal * 3, "FixPIMD2:x_scaled[i]");
  }
  for(int i=0; i<np; i++)
  {
    for(int j=0; j<nlocal; j++)
    {
    x_scaled[i][3*j] = lambda * coords[i][3*j] + (1.0 - lambda) * xc[3*j];
    x_scaled[i][3*j+1] = lambda * coords[i][3*j+1] + (1.0 - lambda) * xc[3*j+1];
    x_scaled[i][3*j+2] = lambda * coords[i][3*j+2] + (1.0 - lambda) * xc[3*j+2];
    }
  }
  
  // test x_scaled
  /* for(int i=0; i<np; i++)
  {
    for(int j=0; j<nlocal; j++)
    {
      fprintf(stdout, "me = %d, i = %d, xscaled: %2.6f, %2.6f, %2.6f.\n", universe->me, i, x_scaled[i][3*j], x_scaled[i][3*j+1], x_scaled[i][3*j+2]);
    }
  }
  fprintf(stdout, "\n");*/
}

/* ---------------------------------------------------------------------- */

int FixPIMD2::pack_forward_comm(int n, int *list, double *buf,
                             int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = comm_ptr[j][0];
    buf[m++] = comm_ptr[j][1];
    buf[m++] = comm_ptr[j][2];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    comm_ptr[i][0] = buf[m++];
    comm_ptr[i][1] = buf[m++];
    comm_ptr[i][2] = buf[m++];
  }
}

/* ----------------------------------------------------------------------
   Memory operations
------------------------------------------------------------------------- */

double FixPIMD2::memory_usage()
{
  double bytes = 0;
  bytes = atom->nmax * size_peratom_cols * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::grow_arrays(int nmax)
{
  if (nmax==0) return;
  int count = nmax*3;

  memory->grow(array_atom, nmax, size_peratom_cols, "FixPIMD2::array_atom");
  memory->grow(nhc_eta,        count, nhc_nchain,   "FixPIMD2::nh_eta");
  memory->grow(nhc_eta_dot,    count, nhc_nchain+1, "FixPIMD2::nh_eta_dot");
  memory->grow(nhc_eta_dotdot, count, nhc_nchain,   "FixPIMD2::nh_eta_dotdot");
  memory->grow(nhc_eta_mass,   count, nhc_nchain,   "FixPIMD2::nh_eta_mass");
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::copy_arrays(int i, int j, int /*delflag*/)
{
  int i_pos = i*3;
  int j_pos = j*3;

  memcpy(nhc_eta       [j_pos], nhc_eta       [i_pos], nhc_size_one_1);
  memcpy(nhc_eta_dot   [j_pos], nhc_eta_dot   [i_pos], nhc_size_one_2);
  memcpy(nhc_eta_dotdot[j_pos], nhc_eta_dotdot[i_pos], nhc_size_one_1);
  memcpy(nhc_eta_mass  [j_pos], nhc_eta_mass  [i_pos], nhc_size_one_1);
}

/* ---------------------------------------------------------------------- */

int FixPIMD2::pack_exchange(int i, double *buf)
{
  int offset=0;
  int pos = i * 3;

  memcpy(buf+offset, nhc_eta[pos],        nhc_size_one_1); offset += nhc_offset_one_1;
  memcpy(buf+offset, nhc_eta_dot[pos],    nhc_size_one_2); offset += nhc_offset_one_2;
  memcpy(buf+offset, nhc_eta_dotdot[pos], nhc_size_one_1); offset += nhc_offset_one_1;
  memcpy(buf+offset, nhc_eta_mass[pos],   nhc_size_one_1); offset += nhc_offset_one_1;

  return size_peratom_cols;
}

/* ---------------------------------------------------------------------- */

int FixPIMD2::unpack_exchange(int nlocal, double *buf)
{
  int offset=0;
  int pos = nlocal*3;

  memcpy(nhc_eta[pos],        buf+offset, nhc_size_one_1); offset += nhc_offset_one_1;
  memcpy(nhc_eta_dot[pos],    buf+offset, nhc_size_one_2); offset += nhc_offset_one_2;
  memcpy(nhc_eta_dotdot[pos], buf+offset, nhc_size_one_1); offset += nhc_offset_one_1;
  memcpy(nhc_eta_mass[pos],   buf+offset, nhc_size_one_1); offset += nhc_offset_one_1;

  return size_peratom_cols;
}

/* ---------------------------------------------------------------------- */

int FixPIMD2::pack_restart(int i, double *buf)
{
  int offset=0;
  int pos = i * 3;
  // pack buf[0] this way because other fixes unpack it
  buf[offset++] = size_peratom_cols+1;

  memcpy(buf+offset, nhc_eta[pos],        nhc_size_one_1); offset += nhc_offset_one_1;
  memcpy(buf+offset, nhc_eta_dot[pos],    nhc_size_one_2); offset += nhc_offset_one_2;
  memcpy(buf+offset, nhc_eta_dotdot[pos], nhc_size_one_1); offset += nhc_offset_one_1;
  memcpy(buf+offset, nhc_eta_mass[pos],   nhc_size_one_1); offset += nhc_offset_one_1;

  return size_peratom_cols+1;
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values
  // unpack the Nth first values this way because other fixes pack them

  int m = 0;
  for (int i=0; i<nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  int pos = nlocal * 3;

  memcpy(nhc_eta[pos],        extra[nlocal]+m, nhc_size_one_1); m += nhc_offset_one_1;
  memcpy(nhc_eta_dot[pos],    extra[nlocal]+m, nhc_size_one_2); m += nhc_offset_one_2;
  memcpy(nhc_eta_dotdot[pos], extra[nlocal]+m, nhc_size_one_1); m += nhc_offset_one_1;
  memcpy(nhc_eta_mass[pos],   extra[nlocal]+m, nhc_size_one_1); m += nhc_offset_one_1;

  nhc_ready = true;
}

/* ---------------------------------------------------------------------- */

int FixPIMD2::maxsize_restart()
{
  return size_peratom_cols+1;
}

/* ---------------------------------------------------------------------- */

int FixPIMD2::size_restart(int /*nlocal*/)
{
  return size_peratom_cols+1;
}

/* ----------------------------------------------------------------------
   Compute centroid-virial kinetic energy estimator
------------------------------------------------------------------------- */

void FixPIMD2::compute_t_vir()
{
  int nlocal = atom->nlocal;
/*  xc = (double*) memory->srealloc(xc, sizeof(double) * nlocal * 3, "FixPIMD2:xc");
  for(int i=0; i<nlocal; i++)
  {
    xc[3*i] = xc[3*i+1] = xc[3*i+2] = 0.0;
    for(int j=0; j<np; j++)
    {
      xc[3*i] += coords[j][3*i];
      xc[3*i+1] += coords[j][3*i+1];
      xc[3*i+2] += coords[j][3*i+2];
    }
    xc[3*i] /= np;
    xc[3*i+1] /= np;
    xc[3*i+2] /= np;
  }
*/
  xf = 0.0;
  t_vir = 0.0;
  for(int i=0; i<nlocal; i++)
  {
    for(int j=0; j<3; j++)
    {
 //     xf -= (atom->x[i][j] - xc[3*i+j]) * atom->f[i][j];
      xf -= atom->x[i][j] * atom->f[i][j];
    }
  }

  MPI_Allreduce(&xf,&t_vir,1,MPI_DOUBLE,MPI_SUM,universe->uworld);
  // fprintf(stdout, "me = %d, before const, xf = %2.6e, t_vir = %2.6e.\n", universe->me, xf, t_vir); 
  t_vir *= 0.5;
  // t_vir += 1.5 * atom->natoms * force->boltz * temp;
  // fprintf(stdout, "me = %d, atom->x: %2.6f, %2.6f, %2.6f\n", universe->me, atom->x[0][0], atom->x[0][1], atom->x[0][2]);
  // fprintf(stdout, "me = %d, atom->f: %2.6f, %2.6f, %2.6f\n", universe->me, atom->f[0][0], atom->f[0][1], atom->f[0][2]);
  // fprintf(stdout, "me = %d, after const, xf = %2.6e, t_vir = %2.6e.\n", universe->me, xf, t_vir); 
}

/* ----------------------------------------------------------------------
   Compute primitive kinetic energy estimator
------------------------------------------------------------------------- */

void FixPIMD2::compute_t_prim()
{
  // fprintf(stdout, "in compute_t_prim, me = %d, N = %d, np = %d, force->boltz = %2.8f, temp = %2.8f, total_spring_energy = %2.8e.\n", universe->me, atom->natoms, np, force->boltz, temp, total_spring_energy);
  t_prim = 1.5 * atom->natoms * np * force->boltz * temp - total_spring_energy;
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::compute_hope()
{
  int nlocal = atom->nlocal;
  double *_mass = atom->mass;
  int *type = atom->type;
  double hoe = 0.0;
  if(harmonicflag)
  {
    for(int i=0; i<nlocal; i++)
    {
      double ff = _mass[type[i]] * omega * omega * force->mvv2e;
      for(int j=0; j<3; j++)
      {
        hoe += 0.5 * ff * atom->x[i][j] * atom->x[i][j]; 
      }
    }
  }
  hoe /= np;
  MPI_Allreduce(&hoe, &hope, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  hope += total_spring_energy;
}

/* ---------------------------------------------------------------------- */

void FixPIMD2::compute_tote()
{
  double kine = 0.0;
  int nlocal = atom->nlocal;
  // double *_mass = atom->mass;
  int *type = atom->type;
  for(int i=0; i<nlocal; i++)
  {
    for(int j=0; j<3; j++)
    {
      kine += 0.5 * mass[type[i]] * atom->v[i][j] * atom->v[i][j];
    }
  }
  MPI_Allreduce(&kine, &totke, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
  totke *= force->mvv2e;
  tote = totke + hope;
}

/* ----------------------------------------------------------------------
   Compute free energy derivative estimator for thermodynamic integration
------------------------------------------------------------------------- */

void FixPIMD2::compute_dfdl()
{
  if(timethod==MSTI)
  {
    int nlocal = atom->nlocal;
    double *_mass = atom->mass;
    int *type = atom->type;
    double pote = 0.0;
    for(int i=0; i<nlocal; i++)
    {
      for(int j=0; j<3; j++)
      {
        pote += 0.5 * _mass[type[i]] * omega * omega * atom->x[i][j] * atom->x[i][j];
      }
    }
    pote /= np;
    MPI_Allreduce(&pote, &dfdl, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
    for(int i=0; i<nlocal; i++)
    {
      for(int j=0; j<3; j++)
      {
        dfdl -= 0.5 * _mass[type[i]] * omega * omega * xc[3*i+j] * xc[3*i+j];
      }
    }
    dfdl *= force->mvv2e;
  }
  else if(timethod==SCTI)
  {
    int nlocal = atom->nlocal;
    double *_mass = atom->mass;
    int *type = atom->type;
    double pote = 0.0;
    for(int i=0; i<nlocal; i++)
    {
      for(int j=0; j<3; j++)
      {
        pote += (atom->x[i][j] - xc[3*i+j]) * _mass[type[i]] * omega * omega * x_scaled[universe->iworld][3*i+j];
      }
    }
    pote /= np;
    MPI_Allreduce(&pote, &dfdl, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
    dfdl *= force->mvv2e;
  }
}

/* ---------------------------------------------------------------------- */

double FixPIMD2::compute_vector(int n)
{
  if(n==0) { return total_spring_energy; }
//  if(n==1) { return t_sys; }
  if(n==1) { return totke; }
  if(n==2) { return hope; }
  if(n==3) { return tote; }
  if(n==4) { return t_prim; }
  if(n==5) { return t_vir; }
//  if(n==6) { return atom->v[0][0]; }
//  if(n==7) { return r1; }
  if(n==6) { return dfdl; }
  return 0.0;
}
