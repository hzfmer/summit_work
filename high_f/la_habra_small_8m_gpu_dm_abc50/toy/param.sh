 -X 3456 -Y 3456 -Z 1400 -x 24 -y 24 -G 1
 --TMAX 30.0 --DH 8.0 --DT 0.0004 --ND 80 --ARBC 0.95
 --IDYNA 0 --NSRC 0 --NVAR 3 
 --MEDIASTART 2 --NTISKP 50 --WRITE_STEP 500
 --INVEL mesh --NVE 1 
 -c output_ckp/ckp 
 --SOURCEFILE input/source.txt 
 --RECVFILE input/receiver.txt
 --INTOPO topography.bin
 --FAC 1.0 --Q0 150. --EX 0.0 --FP 1.0
