The files mpole_k[1-3]*.seq are madx sequences containing various multipoles.

The corresponding file mpole_k[1-3]*.npy are the results of transporting a
set of particles through the multipole with the script madx script
file madx_track_multipole.madx.

The script run_multipole_sequences.sh runs madx on all the .seq files and
generates the .npy files.

These are used in the test_madx_multipoles.py test that checks that
synergia tracking of madx multipoles matches madx tracking of madx multipoles,
making sure that we have the skew and tilt definitions correct.
