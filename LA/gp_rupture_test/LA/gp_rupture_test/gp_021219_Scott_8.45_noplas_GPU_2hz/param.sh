  -X 6320 -Y 4200 -Z 400 -x 40 -y 12 -G 1
  --TMAX 400.0 --DH 100. --DT 0.005
  --NPC 0 --ND 50 --ARBC 0.95
  --IFAULT 2 --MEDIASTART 2 --READ_STEP 100
  --NSRC 1125520 --NST 40000
  --INVEL mesh --NVAR 3
  --INSRC input_rst/srcpart/tpsrc/tpsrc
  --INSRC_I2 input_rst/srcpart/split_faults/fault
  --NVE 1 
  --WRITE_STEP 200 --NTISKP 10
  --NBGX 1 --NEDX 6320 
  --NBGY 1 --NEDY 4200
  --NBGZ 1 --NEDZ 1
  -c output_ckp/chkp -o output_sfc
  --FOLLOWBATHY 0 --SoCalQ 1
  --FAC 1.0 --Q0 150. --EX 0.6 --FP 0.5