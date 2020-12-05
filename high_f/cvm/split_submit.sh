#!/bin/bash -l
for F in {1..7}; do cp ucvm_la_habra_ext_large_0.conf ucvm_la_habra_ext_large_${F}.conf; sed -i -e "/z0/s/.*/z0=$((384*20*F)).0/" -e "/meshfile/s/.*/meshfile=la_habra_ext_large_cvmsi_20m_${F}.media/" -e "/gridfile/s/.*/gridfile=la_habra_ext_large_cvmsi_20m_${F}.grid/" ucvm_la_habra_ext_large_${F}.conf; done
#sleep 4h; for F in 16 17 18 {22..31}; do qsub extract_lahabra_large_"$F".pbs; sleep 2h; done
