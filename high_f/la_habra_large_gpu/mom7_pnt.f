      parameter(ntall=4300,dt=0.001,pi=3.1415927)
      real src
      dimension src(ntall),src1(ntall)
      open(8,file='source.asc')
      open(9,file=
     +'momrate_src_yz_av_20m_4hz_0',
     +access='direct',recl=4*(3+6*ntall))


      str=(239.-39.9)*pi/180.
      dip=70.*pi/180.
      rakes=38.*pi/180.
      amom=5.764e16
      n1=749
      n2=1116
      n3=250
      
      amom1=0
      nrec=0
      amom1=amom1+amom
      read(8,*) (src(i),i=1,ntall)

      src1(1)=0.
      do 60 i=2,ntall
	      src1(i)=dt*src(i)+src1(i-1)
   60 continue
      ff=amom/src1(ntall)
      print*,ff,src1(ntall),n1,n2,n3
      do l=1,ntall
	      src(l)=src(l)*ff 
      enddo

        ayy_s= -(sin(dip)*cos(rakes)*sin(2.*str)+
     +  sin(2.*dip)*sin(rakes)*sin(str)*sin(str))

        axy_s= sin(dip)*cos(rakes)*cos(2.*str)+
     +  0.5*(sin(2.*dip)*sin(rakes)*sin(2.*str))

        ayz_s= (cos(dip)*cos(rakes)*cos(str)+
     +  cos(2.*dip)*sin(rakes)*sin(str))

        axx_s= sin(dip)*cos(rakes)*sin(2.*str)-
     +  sin(2.*dip)*sin(rakes)*cos(str)*cos(str)

        axz_s=  (cos(dip)*cos(rakes)*sin(str)-
     +  cos(2.*dip)*sin(rakes)*cos(str))

        azz_s= sin(2.*dip)*sin(rakes)

          xx_s=axx_s
          yy_s=ayy_s
          zz_s=azz_s
          xz_s=axz_s
          yz_s=ayz_s
          xy_s=axy_s

      nrec=nrec+1
      write(9,rec=nrec) 
     +n1,n2,n3,(0.25*xx_s*src(l),0.25*yy_s*src(l),0.25*zz_s*src(l),
     +0.25*xz_s*src(l),1.00*yz_s*src(l),0.25*xy_s*src(l),l=1,ntall)

      nrec=nrec+1
      write(9,rec=nrec) 
     +n1,n2+1,n3,(0.25*xx_s*src(l),.25*yy_s*src(l),.25*zz_s*src(l),
     +0.25*xz_s*src(l),0.*yz_s*src(l),0.*xy_s*src(l),l=1,ntall)

      nrec=nrec+1
      write(9,rec=nrec) 
     +n1+1,n2,n3,(0.*xx_s*src(l),0.*yy_s*src(l),0.*zz_s*src(l),
     +0.25*xz_s*src(l),0.*yz_s*src(l),0.25*xy_s*src(l),l=1,ntall)

      nrec=nrec+1
      write(9,rec=nrec) 
     +n1,n2,n3-1,(0.25*xx_s*src(l),0.25*yy_s*src(l),0.25*zz_s*src(l),
     +0.*xz_s*src(l),0.*yz_s*src(l),0.25*xy_s*src(l),l=1,ntall)

      nrec=nrec+1
      write(9,rec=nrec) 
     +n1,n2+1,n3-1,(0.25*xx_s*src(l),0.25*yy_s*src(l),0.25*zz_s*src(l),
     +0.*xz_s*src(l),0.*yz_s*src(l),0.*xy_s*src(l),l=1,ntall)

      nrec=nrec+1
      write(9,rec=nrec) 
     +n1+1,n2,n3-1,(0.*xx_s*src(l),0.*yy_s*src(l),0.*zz_s*src(l),
     +0.*xz_s*src(l),0.*yz_s*src(l),0.25*xy_s*src(l),l=1,ntall)

      nrec=nrec+1
      write(9,rec=nrec) 
     +n1+1,n2+1,n3,(0.*xx_s*src(l),0.*yy_s*src(l),0.*zz_s*src(l),
     +0.25*xz_s*src(l),0.*yz_s*src(l),0.*xy_s*src(l),l=1,ntall)

      print*,amom1
      print*,nrec
      end
