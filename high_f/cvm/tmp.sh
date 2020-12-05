#!/bin/bash -l
sleep 4h; for F in 16 17 18 {22..31}; do qsub extract_lahabra_large_"$F".pbs; sleep 2h; done
