  -X 80 -Y 80 -Z 100 -x 4 -y 4 -G 1
  --TMAX 64.00   --DH 4 --DT 0.00032 
  --ND 20 --ARBC 0.95 --NPC 2
  --NSRC 1 --NST 200000 --READ_STEP 200000
  --IFAULT 6 --INSRC stf.txt  
  --INVEL mesh --MEDIASTART 2 
  --IDYNA 0  --NVE 1 --NVAR 5 
  --NBGX 40 --NEDX 40 --NBGY 41 --NEDY 41 
  --NBGZ 1 --NEDZ 37 --NSKPZ 18 --OUT output_sfc 
  --NTISKP 25 --WRITE_STEP 8000
  --FAC 1.0 --Q0 150. --EX 0 --FP 1.0