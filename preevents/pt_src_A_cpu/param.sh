  -X 1000 -Y 480 -Z 800 -x 40 -y 20 
  --TMAX 5.00   --DH 5.0 --DT 0.0004 
  --ND 20 --ARBC 0.95 
  --NSRC 1 --NST 2501 --READ_STEP 2501 --WRITE_STEP 1250
  --IFAULT 1 --INSRC source 
  --INVEL mesh
  --IDYNA 0  --NVE 1 --NVAR 5 
  --MEDIASTART 2 
  --NBGX2 101 --NEDX2 901 --NSKPX2 50
  --NBGY2 113 --NEDY2 413 --NSKPY2 25
  --NBGZ2 1 --NEDZ2 40 --NSKPZ2 1
  --NBGX 1 --NEDX 1000 --NBGY 1 --NEDY 480 
  --NBGZ 1 --NEDZ 1  --OUT output_sfc 
  --NTISKP 10 --NTISKP2 10 
  --FAC 1.0 --Q0 150. --EX 0.2 --FP 1.0
