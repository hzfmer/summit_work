 -X 3456 -Y 3456 -Z 108 -x 12 -y 6 -G 1
 --TMAX 10.0 --DH 8.0 --DT 0.00023 --ND 80 --ARBC 0.95
 --IDYNA 0 --NSRC 0 --NVAR 3 
 --MEDIASTART 2 --NTISKP 20 --WRITE_STEP 200
 --INVEL input/mesh --NVE 1 
 --NBGX 1 --NEDX 3456 --NBGY 1 --NEDY 3456
 --NBGZ 1 --NEDZ 1
 -c output_ckp/ckp -o output_sfc
 --SOURCEFILE input/source.txt --INTOPO topography.bin
 --FAC 1.0 --Q0 150. --EX 0.6 --FP 1.0
