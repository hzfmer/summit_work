 -X 240 -Y 150 -Z 128,64 -x 3 -y 1 -G 2
 --TMAX 8.0 --DH 45.0 --DT 0.001 --ND 10 --ARBC 0.95
 --IDYNA 0 --NSRC 0,0
 --MEDIASTART 2 --NTISKP 20 --WRITE_STEP 200
 --INVEL mesh --NVE 1 --NVAR 3
 --NBGX 1,1 --NEDX 720,240 --NSKPX 4,2
 --NBGY 1,1 --NEDY 450,150 --NSKPY 2,2
 --NBGZ 1,1 --NEDZ 128,64 --NSKPZ 2,2
 -c output_ckp/ckp -o output_sfc_3_1
 --SOURCEFILE input/source.txt
 --INTOPO topography.bin
 --FAC 1.0 --Q0 150. --EX 0.0 --FP 1.0