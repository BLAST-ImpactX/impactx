#!/bin/bash

set -e

seqs="quad quad_long drift rfc sext sext_long oct oct_long solenoid
      kicker hkicker vkicker kicker_long hkicker_long vkicker_long"

for s in $seqs
do
    cp l_$s.seq element.seq
    madx madx_track_element.madx
    mv madx_track_element.txtone madx_$s.out
    rm element.seq
    #python convert_madx_tracks_to_np.py madx_track_multipole.txtone
    #mv madx_track_multipole.npy $s.npy
done
