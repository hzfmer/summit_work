  -X 400 -Y 400 -Z 400 -x 8 -y 4 -G 1
  --TMAX 150.00   --DH 2.5 --DT 0.0002
  --ND 40 --ARBC 0.95 --NPC 2
  --NSRC 1 --NST 750000 --READ_STEP 750000
  --IFAULT 6 --INSRC stf.txt  
  --INVEL mesh --MEDIASTART 2 
  --IDYNA 0  --NVE 1 --NVAR 5 
  --NBGX 1 --NEDX 351 --NSKPX 50
  --NBGY 1 --NEDY 301 --NSKPY 100 
  --NBGZ 1 --NEDZ 41 --NSKPZ 40 --OUT output_sfc 
  --NTISKP 1 --WRITE_STEP 750000
  --FAC 1.0 --Q0 150. --EX 0.2 --FP 1.0
