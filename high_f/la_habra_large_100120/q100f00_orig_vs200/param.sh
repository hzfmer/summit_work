  -X 2160 -Y 1656 -Z 160,380,680 -x 48 -y 36 -G 1 
  --TMAX 120.00 --DH 72.0 --DT 0.0004 
  --ND 44 --ARBC 0.95 
  --NSRC 0,15625,0 --NST 12500 --READ_STEP 12500
  --IFAULT 1 --INSRC source
  --MEDIASTART 2 --IDYNA 0
  --INVEL mesh --NVE 1 --NVAR 3 
  --NBGX 1,1,1 --NEDX 9504,3168 --NSKPX 1,8
  --NBGY 1,1,1 --NEDY 7020,2340 --NSKPY 1,8
  --NBGZ 1,1,1 --NEDZ 1,1,1 
  --NTISKP 50 --WRITE_STEP 200
  -c output_ckp/ckp -o output_sfc
  --FAC 1.0 --Q0 150. --EX 0.0 --FP 1.0
