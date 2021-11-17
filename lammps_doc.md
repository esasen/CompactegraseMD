# Interesting LAMMPS features

##### Create Molecule Template: [molecule command](https://docs.lammps.org/molecule.html)

##### Rigid Bodies: [fix rigid](https://docs.lammps.org/fix_rigid.html)

##### Restart [read_restart command](https://docs.lammps.org/read_restart.html)

[fix pour](https://docs.lammps.org/fix_pour.html)

- Insert finite-size particles or molecules into the simulation box every few timesteps within a specified region until N particles or molecules have been inserted

[fix deposit](https://docs.lammps.org/fix_deposit.html)

- Insert a single atom or molecule into the simulation domain every M timesteps until N atoms or molecules have been inserted.

## To read





# Features

#### Grand Canonical Monte Carlo: [fix gcmc](https://docs.lammps.org/fix_gcmc.html) 

mol keyword:

​		 set: 		translations and rotations only on full molecules

​	 	not set:  translations on constituent atoms

rigid body insertion:  set rigid with value id of fix rigid/small (checkout what that means)



- temperature of imaginary reservoir and fix nvt should be the same
-  Some fixes have an associated potential energy. Examples of such fixes include: [efield](https://docs.lammps.org/fix_efield.html), [gravity](https://docs.lammps.org/fix_gravity.html), [addforce](https://docs.lammps.org/fix_addforce.html), [langevin](https://docs.lammps.org/fix_langevin.html), [restrain](https://docs.lammps.org/fix_restrain.html), [temp/berendsen](https://docs.lammps.org/fix_temp_berendsen.html), [temp/rescale](https://docs.lammps.org/fix_temp_rescale.html), and [wall fixes](https://docs.lammps.org/fix_wall.html). For that energy to be included in the total potential energy of the system (the quantity used when performing GCMC exchange and MC moves), you MUST enable the [fix_modify](https://docs.lammps.org/fix_modify.html) *energy* option for that fix.  The doc pages for individual [fix](https://docs.lammps.org/fix.html) commands specify if this should be done.
- Use of this fix typically will cause the number of atoms to fluctuate, therefore, you will want to use the [compute_modify dynamic/dof](https://docs.lammps.org/compute_modify.html) command to insure that the current number of atoms is used as a normalizing factor each time temperature is computed. A simple example of this is:
  - ​	compute_modify thermo_temp dynamic yes ()
- timestep cannot be changed when restarting
- when using shake or rigid only GCMC exchange moves are supported: argument M must be set to zero

