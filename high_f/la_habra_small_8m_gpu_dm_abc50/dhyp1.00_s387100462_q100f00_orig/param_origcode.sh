 -X 1152 -Y 1152 -Z 108,464 -x 12 -y 12 -G 2
 --TMAX 30.0 --DH 24.0 --DT 0.001 --ND 80 --ARBC 0.95
 --IFAULT 1 --INSRC source
 --IDYNA 0 --NSRC 0,15625 --NST 5000 
 --READ_STEP 5000 --WRITE_STEP 300
 --MEDIASTART 2 --NTISKP 20 
 --INVEL mesh --NVE 1 --NVAR 3 
 --NBGX 1,1 --NEDX 3456,1152 --NBGY 1,1 --NEDY 3456,1152
 --NBGZ 1,1 --NEDZ 1,1
 -c output_ckp/ckp -o output_sfc
 --FAC 1.0 --Q0 150. --EX 0.0 --FP 1.0