SCALE Job started on Sun 2023-07-02 03:32:39
With Host name owl
working directory /tmp/scale.o.218476
and process id 218476.
Input file named /home/o/git/IMSR/SCALE_IMSR_Model/01-vol_mass_count/EIRENE.inp
and output file named /home/o/git/IMSR/SCALE_IMSR_Model/01-vol_mass_count/EIRENE.out
----------------------------------------------------------
>>> Loading SCALE Standard Composition Library from /opt/scale6.3.0_data/scale.rev40.sclib
*** Adjusted number of batches to 36 to be divisible by the number of nodes (36)
Calculating volume batches 1-36 of 36
Tracing 50013184 rays per core along the Z axis
[owl:218531] 35 more processes have sent help message help-mpi-btl-base.txt / btl:no-nics
[owl:218531] Set MCA parameter "orte_base_help_aggregate" to 0 to see all help / error messages
Reinitializing LFG control, ignoring new seed 2272013 (actual seed 112358) and resetting stream ID
*** Disabling parallel HDF5 calls
::: Building model
>>> Splitting compound 3006000 into at% {3006012: 0.9893, 3006013: 0.0107}
>>> Splitting compound 6000 into at% {6012: 0.9893, 6013: 0.0107}
>>> Creating default boundary mesh from (-215 -215 -30) to (215 215 460) for ORANGE geometry
::: Building physics 'physics_ce_mc'
>>> Loading CE library /opt/scale6.3.1-mpi/data/ce_v7.1_endf.xml
Loading      nuclide       li-6 @ 900K
Loading      nuclide       li-6 @ 1200K
Broadening   nuclide       li-6 @ 923.15K
Loading      nuclide       li-7 @ 900K
Loading      nuclide       li-7 @ 1200K
Broadening   nuclide       li-7 @ 923.15K
Loading      nuclide       be-9 @ 900K
Loading      nuclide       be-9 @ 1200K
Broadening   nuclide       be-9 @ 923.15K
Loading      nuclide          c @ 900K
Loading      nuclide          c @ 1200K
Broadening   nuclide          c @ 923.15K
Loading      nuclide c-graphite @ 800K
Loading      nuclide c-graphite @ 1000K
Broadening   nuclide c-graphite @ 923.15K
Loading      nuclide       f-19 @ 900K
Loading      nuclide       f-19 @ 1200K
Broadening   nuclide       f-19 @ 923.15K
Loading      nuclide      si-28 @ 900K
Loading      nuclide      si-28 @ 1200K
Broadening   nuclide      si-28 @ 923.15K
Loading      nuclide      si-29 @ 900K
Loading      nuclide      si-29 @ 1200K
Broadening   nuclide      si-29 @ 923.15K
Loading      nuclide      si-30 @ 900K
Loading      nuclide      si-30 @ 1200K
Broadening   nuclide      si-30 @ 923.15K
Loading      nuclide       p-31 @ 900K
Loading      nuclide       p-31 @ 1200K
Broadening   nuclide       p-31 @ 923.15K
Loading      nuclide      cr-50 @ 900K
Loading      nuclide      cr-50 @ 1200K
Broadening   nuclide      cr-50 @ 923.15K
Loading      nuclide      cr-52 @ 900K
Loading      nuclide      cr-52 @ 1200K
Broadening   nuclide      cr-52 @ 923.15K
Loading      nuclide      cr-53 @ 900K
Loading      nuclide      cr-53 @ 1200K
Broadening   nuclide      cr-53 @ 923.15K
Loading      nuclide      cr-54 @ 900K
Loading      nuclide      cr-54 @ 1200K
Broadening   nuclide      cr-54 @ 923.15K
Loading      nuclide      mn-55 @ 900K
Loading      nuclide      mn-55 @ 1200K
Broadening   nuclide      mn-55 @ 923.15K
Loading      nuclide      fe-54 @ 900K
Loading      nuclide      fe-54 @ 1200K
Broadening   nuclide      fe-54 @ 923.15K
Loading      nuclide      fe-56 @ 900K
Loading      nuclide      fe-56 @ 1200K
Broadening   nuclide      fe-56 @ 923.15K
Loading      nuclide      fe-57 @ 900K
Loading      nuclide      fe-57 @ 1200K
Broadening   nuclide      fe-57 @ 923.15K
Loading      nuclide      fe-58 @ 900K
Loading      nuclide      fe-58 @ 1200K
Broadening   nuclide      fe-58 @ 923.15K
Loading      nuclide      ni-58 @ 900K
Loading      nuclide      ni-58 @ 1200K
Broadening   nuclide      ni-58 @ 923.15K
Loading      nuclide      ni-60 @ 900K
Loading      nuclide      ni-60 @ 1200K
Broadening   nuclide      ni-60 @ 923.15K
Loading      nuclide      ni-61 @ 900K
Loading      nuclide      ni-61 @ 1200K
Broadening   nuclide      ni-61 @ 923.15K
Loading      nuclide      ni-62 @ 900K
Loading      nuclide      ni-62 @ 1200K
Broadening   nuclide      ni-62 @ 923.15K
Loading      nuclide      ni-64 @ 900K
Loading      nuclide      ni-64 @ 1200K
Broadening   nuclide      ni-64 @ 923.15K
Loading      nuclide      mo-92 @ 900K
Loading      nuclide      mo-92 @ 1200K
Broadening   nuclide      mo-92 @ 923.15K
Loading      nuclide      mo-94 @ 900K
Loading      nuclide      mo-94 @ 1200K
Broadening   nuclide      mo-94 @ 923.15K
Loading      nuclide      mo-95 @ 900K
Loading      nuclide      mo-95 @ 1200K
Broadening   nuclide      mo-95 @ 923.15K
Loading      nuclide      mo-96 @ 900K
Loading      nuclide      mo-96 @ 1200K
Broadening   nuclide      mo-96 @ 923.15K
Loading      nuclide      mo-97 @ 900K
Loading      nuclide      mo-97 @ 1200K
Broadening   nuclide      mo-97 @ 923.15K
Loading      nuclide      mo-98 @ 900K
Loading      nuclide      mo-98 @ 1200K
Broadening   nuclide      mo-98 @ 923.15K
Loading      nuclide     mo-100 @ 900K
Loading      nuclide     mo-100 @ 1200K
Broadening   nuclide     mo-100 @ 923.15K
Loading      nuclide     gd-152 @ 900K
Loading      nuclide     gd-152 @ 1200K
Broadening   nuclide     gd-152 @ 923.15K
Loading      nuclide     gd-155 @ 900K
Loading      nuclide     gd-155 @ 1200K
Broadening   nuclide     gd-155 @ 923.15K
Loading      nuclide     gd-156 @ 900K
Loading      nuclide     gd-156 @ 1200K
Broadening   nuclide     gd-156 @ 923.15K
Loading      nuclide     gd-157 @ 900K
Loading      nuclide     gd-157 @ 1200K
Broadening   nuclide     gd-157 @ 923.15K
Loading      nuclide     gd-158 @ 900K
Loading      nuclide     gd-158 @ 1200K
Broadening   nuclide     gd-158 @ 923.15K
Loading      nuclide     gd-160 @ 900K
Loading      nuclide     gd-160 @ 1200K
Broadening   nuclide     gd-160 @ 923.15K
Loading      nuclide      u-234 @ 900K
Loading      nuclide      u-234 @ 1200K
Broadening   nuclide      u-234 @ 923.15K
Loading      nuclide      u-235 @ 900K
Loading      nuclide      u-235 @ 1200K
Broadening   nuclide      u-235 @ 923.15K
Loading      nuclide      u-236 @ 900K
Loading      nuclide      u-236 @ 1200K
Broadening   nuclide      u-236 @ 923.15K
Loading      nuclide      u-238 @ 900K
Loading      nuclide      u-238 @ 1200K
Broadening   nuclide      u-238 @ 923.15K
*** The following nuclide IDs were remapped due to missing neutron data
*** Remapped 4 nuclide IDs: 6012->6000, 6013->6000, 3006012->3006000, 3006013->3006000
::: Building tallies
::: Building sources
Building source 'cuboid_fission'
::: Building Shift solver internals
::: Building Shift sources
::: Initializing Shift solver
::: Building Shift tallies
::: Skipping transport due to user option
::: Run complete


***Warning: This problem will not be run because parm=check was entered in the analytical sequence specification.
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
thread 0HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2)   #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5G.c line 893 in H5Gclose(): decrementing group ID failed
    major: Symbol table
    minor: Unable to decrement reference count
  #001: H5Iint.c line 1155 in H5I_dec_app_ref(): can't decrement ID ref count
    major: Object ID
    minor: Unable to decrement reference count
  #002: H5Iint.c line 1107 in H5I__dec_app_ref(): can't decrement ID ref count
    major: Object ID
    minor: Unable to decrement reference count
  #003: H5Iint.c line 1011 in H5I__dec_ref(): can't locate ID
    major: Object ID
    minor: Unable to find ID information (already closed?)
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5G.c line 893 in H5Gclose(): decrementing group ID failed
    major: Symbol table
    minor: Unable to decrement reference count
  #001: H5Iint.c line 1155 in H5I_dec_app_ref(): can't decrement ID ref count
    major: Object ID
    minor: Unable to decrement reference count
  #002: H5Iint.c line 1107 in H5I__dec_app_ref(): can't decrement ID ref count
    major: Object ID
    minor: Unable to decrement reference count
  #003: H5Iint.c line 1011 in H5I__dec_ref(): can't locate ID
    major: Object ID
    minor: Unable to find ID information (already closed?)
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5F.c line 982 in H5Fflush(): unable to synchronously flush file
    major: File accessibility
    minor: Unable to flush data from cache
  #001: H5F.c line 946 in H5F__flush_api_common(): invalid object identifier
    major: Invalid arguments to routine
    minor: Inappropriate type
  #002: H5VLint.c line 1779 in H5VL_vol_object(): invalid identifier
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5F.c line 982 in H5Fflush(): unable to synchronously flush file
    major: File accessibility
    minor: Unable to flush data from cache
  #001: H5F.c line 946 in H5F__flush_api_common(): invalid object identifier
    major: Invalid arguments to routine
    minor: Inappropriate type
  #002: H5VLint.c line 1779 in H5VL_vol_object(): invalid identifier
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5I.c line 550 in H5Iget_ref(): can't get ID ref count
    major: Object ID
    minor: Can't get value
  #001: H5Iint.c line 1361 in H5I_get_ref(): can't locate ID
    major: Object ID
    minor: Unable to find ID information (already closed?)
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5F.c line 1060 in H5Fclose(): decrementing file ID failed
    major: File accessibility
    minor: Unable to close file
  #001: H5Iint.c line 1155 in H5I_dec_app_ref(): can't decrement ID ref count
    major: Object ID
    minor: Unable to decrement reference count
  #002: H5Iint.c line 1107 in H5I__dec_app_ref(): can't decrement ID ref count
    major: Object ID
    minor: Unable to decrement reference count
  #003: H5Iint.c line 1011 in H5I__dec_ref(): can't locate ID
    major: Object ID
    minor: Unable to find ID information (already closed?)
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5I.c line 550 in H5Iget_ref(): can't get ID ref count
    major: Object ID
    minor: Can't get value
  #001: H5Iint.c line 1361 in H5I_get_ref(): can't locate ID
    major: Object ID
    minor: Unable to find ID information (already closed?)
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5F.c line 1060 in H5Fclose(): decrementing file ID failed
    major: File accessibility
    minor: Unable to close file
  #001: H5Iint.c line 1155 in H5I_dec_app_ref(): can't decrement ID ref count
    major: Object ID
    minor: Unable to decrement reference count
  #002: H5Iint.c line 1107 in H5I__dec_app_ref(): can't decrement ID ref count
    major: Object ID
    minor: Unable to decrement reference count
  #003: H5Iint.c line 1011 in H5I__dec_ref(): can't locate ID
    major: Object ID
    minor: Unable to find ID information (already closed?)
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
HDF5-DIAG: Error detected in HDF5 (1.14.1-2) thread 0:
  #000: H5T.c line 1911 in H5Tclose(): not a datatype
    major: Invalid arguments to routine
    minor: Inappropriate type
