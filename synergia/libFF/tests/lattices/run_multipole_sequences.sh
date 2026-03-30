#!/bin/bash

seqs="mpole_k1 mpole_k1s mpole_k1_tilt mpole_k2 mpole_k2s mpole_k2_tilt mpole_k3
      mpole_k3s mpole_k3_tilt mpole_k4 mpole_k4s mpole_k4_tilt mpole_k5 mpole_k5s
      mpole_k5_tilt mpole_k6 mpole_k6s mpole_k6_tilt cfsbend_dipole cfsbend_dipole_skew cfsbend_sext cfsbend_sext_skew"

for s in $seqs
do
    cp $s.seq mpole1.seq
    madx madx_track_multipole.madx
    mv madx_track_multipole.txtone madx_$s.out
    rm mpole1.seq
    #python convert_madx_tracks_to_np.py madx_track_multipole.txtone
    #mv madx_track_multipole.npy $s.npy
done
