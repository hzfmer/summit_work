  -X 864 -Y 720 -Z 184,372 -x 12 -y 8 -G 2 
  --TMAX 120.00 --DH 24.0 --DT 0.0005 
  --ND 50 --ARBC 0.95 
  --NSRC 0,15625 --NST 10000 --READ_STEP 10000
  --IFAULT 1 --INSRC source
  --MEDIASTART 2 --IDYNA 0
  --INVEL mesh --NVE 1 --NVAR 3 
  --NBGX 1,1 --NEDX 2592,864 --NSKPX 1,8
  --NBGY 1,1 --NEDY 2160,720 --NSKPY 1,8
  --NBGZ 1,1 --NEDZ 1,1 
  --NTISKP 40 --WRITE_STEP 200
  -c output_ckp/ckp -o output_sfc
  --RECVFILE receiver.txt
  --FAC 1.0 --Q0 150. --EX 0.0 --FP 1.0
