  -X 1400 -Y 1400 -Z 600  -x 20 -y 20 -G 1 
  --TMAX 30.00   --DH  20.0 --DT 0.001 
  --ND 20 --ARBC 0.95 
  --NSRC 1560 --NST 5000 --READ_STEP 5000 --WRITE_STEP 500
  --IFAULT 1 --INSRC momrate_LaHabra.zf.100m.rev.bin   
  --INVEL mesh  
  --IDYNA 0  --NVE 1 --NVAR 3 
  --MEDIASTART 2 
  --NBGX 1 --NEDX 1400 --NBGY 1 --NEDY 1400
  --NBGZ 1 --NEDZ 1  --OUT output_sfc 
  --NTISKP 10 
  --FAC 1.0 --Q0 150. --EX 0.6 --FP 1.0
