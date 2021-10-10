c c c c c c
c Herman-Skillman Code 1961-62
c

c   re-done version nov 9th 88; Modified 1990-1993; 
c   Most Recent Modifications: 1995, November 1st
c 
c   Herman-Skillman Code 1961-62
c
c   hartree-fock-slater self-consistent atomic field program
c          
c   ***********************************************************
c   ***                                                     ***
c   ***      modified inputs and outputs                    ***
c   ***       ru2  is read from file 13  (441 numbers)      ***
c   ***       ru3                    14  (441 numbers)      ***
c   ***       ru2   written on  file 15  (441 numbers)      ***
c   ***      all other input is  format()                   ***
c   ***      except  the input potential (110 numbers)      ***
c   ***                                                     ***
c   ***********************************************************
c      originally  written by sherwood skillman                         
c      rca laboratories, princeton, new jersey, spring 1961             
c      modified by frank herman, summer 1961                            
c      further modified by richard kortum,  lockheed research           
c       laboratories, palo alto, california,  summer 1962               
c     iteration number,measure of self-consistency,and atomic number z  
c     are always printed on-line.                                       
c     output tape b-5 contains self-consistent atomic potential,        
c      energy eigenvalues, and radial wave functions                    
c      successive solutions separated by end of file                    
c     z = atomic number.  ion = ionicity.  zzz = ion+1                  
c     if key = 0   normalized numerical atomic potential is read in     
c      for every 4th mesh point from 1 to 437.                          
c     if key = 1     numerical atomic potential is to be read in.       
c     if key = 2  extrapolate for starting potential                    
c     ipratt is the number of consecutive iterations to use the pratt   
c      improvement scheme following each application of the arithmetic  
c      average scheme.                                                  
c     maxit = maximum no. of iterations unless maxit is read as 0 ,     
c      in which case maxit is set to 20.                                
c     nocopy =0, b-5 is copied to a-3. nocopy = 1, b-5 is not copied.
c     kut = 1, the cutoff potential is -2(ion+1)/r.              
c     kut = 0, no cutoff potential is imposed.       
c
c     F287 change
c
      Character * 80 outname, outfil, letters, wnlz
      Character * 80 allnam(20), anlz
c
      common LimMesh,xnum(521),rsvale(521),rsatom(521),v(521),vout(521)          
      common snlo(521)  ,r(521)     ,ru(521)    ,ruexch(521),xi(521)   
      common xj(521)    ,snl(22,521),ruinl1(521),rufnl1(521),ruinl2(521)
      common rufnl2(521),ru2(521)   ,ru3(521)   ,nnlz(24)   ,wwnl(24)   
      common nkkk(24)   ,ee(24)     ,a(4,5)                             
      dimension x(521),rscore(521),denm(521),qq(521)                    
      equivalence (xnum,rscore),(rsvale,denm),(rsatom,qq)               
 1001 format(f8.5,9f7.5)                                                
 1003 format(f4.0,3i4)                                                  
 1005 format('  www= ',f4.0,' zzz= ',f4.0,'   z= ',f4.0,
     1' ncores= ',i4,'   nvales= ',i4,   
     X '   ncspvs= ',i4/' control cards incorrect.') 
 1007 format(i4,f4.0,f8.4)                                              
 1008 format(i7,f7.0, 1pe16.7,2(i6,0pf9.3))                             
 1010 format(1e15.7,14e14.7)                                          
 1012 format(i4,2f8.6,5i4)                                              
 1015 format('  key = ',i4) 
c1016 format('Iter ',iter,'  z  del i(d) x(d)  i cut  ')
 1017 format('#  IDENT: ',i4,' E = ',3(1X,f14.6),i4,e14.7,8x,'z') 
 1018 format(f4.0,3i4,56x,'z',i3,i4)                                   
 2018 format('#',f4.0,3i4,56x,'z',i3,i4)                                   
 1019 format(1pe15.7,4e14.7,' z',i3,i4) 
 1020 format (1pe15.7,57x,'z',i3,i4)                                     
 1021 format (72A1)                                                            
 1058 format (A80)
!Shuka 2021.10.7
      istring = 0      
!end Shuka      
      nfiles =0                                                        
 1    continue                                                         
c     read heading card.                                               
      read(*,1058) outname
      read  (*,1021)                                                   
      write (*,1021)                                                   
c     read control cards and input potentials. calculate trial potential
      read  (*, *)key,tol,thresh,mesh,ipratt,maxit,nocopy,kut         
      read(*,*) LimMesh
      ncards=90                                                         
      write (*,1015)key                                               
      if(maxit) 2,2,3                                                  
 2    maxit=20                                                         
 3    continue                                                          
 4    nblock = (mesh)/40                                               
c     construct x mesh and r mesh                                      
      i=1                                                              
      x(i)=0.0                                                         
      r(i)=0.0                                                         
      deltax=0.0025                                                     
      do 6 j=1,nblock                                                 
      do 5 k=1,40                                                      
      i=i+1                                                            
 5    x(i)=x(i-1)+deltax                                               
 6    deltax=deltax+deltax                                            
      if(key-1)8,9,7                                                   
 7    read  (13,1010)(ru2(m),m=1,441)                                   
      read  (14,1010)(ru3(m),m=1,441)                                   
      ze2=-ru2(1)*.5                                                   
      ze3=-ru3(1)*.5                                                   
      go to 10                                                         
c     read in atomic potential                                        
 8    read  (5,1001)(ru2(m),m=1,437,4)                                 
      go to 10                                                        
 9    read  (14,1010)(ru3(m),m=1,441)                                 
      ze3=-ru3(1)*0.5                                                 
 10   read  (5, *)z,ncores,nvales,ion                               
      if (z) 141,1  ,11                                               
 11   nfiles =nfiles+1                                                
      iz=z                                                            
      ncspvs=ncores+nvales                                            
      c=0.88534138 /z**(1. /3. )                                      
      twoion=ion+ion                                                  
c     100 lines
      zzz=ion+1                                                        
c
c     F287
c
      twozzz = zzz+zzz                                                 
      do 12 i=2,mesh                                                  
 12   r(i)=c*x(i)                                                      
 13   read  (5,*)(nnlz(i),wwnl(i),ee(i),i=1,ncspvs)                  
      www=0.0                                                          
      do 14 i=1,ncspvs                                                 
 14   www=www+wwnl(i)                                                  
      if( abs(z+1. -www-zzz)-0.001 ) 16,15,15                          
 15   continue 
c       write (*,1005)www,zzz,z,ncores,nvales,ncspvs                     
      stop                                                             
 16   continue                                                         
      if (key-1) 21,26,17                                              
c     construct atomic potenial                                        
 17   if( abs(ze3-ze2-z+ze3)-0.001 ) 19,19,18                          
 18   write (*,1022)z,ze2,ze3                                          
 1022 format (//'starting potentials and z=',f4.0,' in error ',2f4.0)   
      stop                                           
 19   do 20 i=1,441                    
      ru(i) = ru3(i)+ru3(i)-ru2(i)             
 20   continue                           
      go to 31                                          
 21   twoz=z+z                                                    
      do 22 i=1,437,4                                              
 22   ru(i)=-ru2(i)*twoz                                            
      ru(441)=ru(437)                                                
      ru(445)=ru(437)                                                 
      m=9                                                            
      do25 i=1,437,4                                                  
      m=m-1                                                           
      if(m) 23,24,24                                                  
 23   ru(i+1)=(22. *ru(i)+11. *ru(i+4)-ru(i+8))/32.                    
      ru(i+2)=(10. *ru(i)+15. *ru(i+4)-ru(i+8))/24.                   
      ru(i+3)=( 6. *ru(i)+27. *ru(i+4)-ru(i+8))/32.                   
      m=9                                                             
      go to 25                                                        
 24   ru(i+1)=(21. *ru(i)+14. *ru(i+4)-3. *ru(i+8))/32.               
      ru(i+2)=( 3. *ru(i)+ 6. *ru(i+4)-    ru(i+8))/ 8.                
      ru(i+3)=( 5. *ru(i)+30. *ru(i+4)-3. *ru(i+8))/32.               
 25   continue                                                         
      go to 31                                                         
 26   if( abs(ze3-z)-0.001 )27,27,29                                   
 27   do 28 i=1,441                                                    
      ru(i)=ru3(i)                                                     
 28   continue                                                         
      go to 31                                                         
 29   zoz=z/ze3                                                        
      do 30 i=1,441                                                    
      ru(i)=ru3(i)*zoz                                                 
 30   continue                                                         
 31   v(1)=-9.9e35                                                     
      m=min0(441,mesh)                                                 
      if(kut) 32,37,32                                                 
 32   do 33 i= 2,m                                                     
 33   v(i)= ru(i)/r(i)                                                 
      if(mesh-m)34,34,36                                               
 34   do 35 i=442, mesh                                                
 35   v(i)=-twoion/r(i)                                             
 36   limit=m                                                        
      icut= mesh                                                       
      ic=mesh                                                          
      go to 47                                                         
 37   continue                                                         
      icut=0                                                          
      do 42 i=2,m                                                      
      if (icut) 38,38,40                                               
 38   if( twozzz+ru(i)) 41,41,39                                       
 39   icut =i                                                          
 40   v(i)=-twozzz/r(i)                                                
      go to 42                                                        
 41   v(i)=ru(i)/r(i)                                                 
 42   continue                                                        
      if(icut) 43,43,44                                               
 43   icut=m                                                           
 44   limit=icut                                                       
      if(mesh-m) 47,47,45                                             
 45   continue                                                        
      do 46 i=442,mesh                                                
 46   v(i)=-twozzz/r(i)                                              
 47   continue                                                       
      delta =1000000.                                                
      niter=0                                                          
      nonmon =3                                                       
      iprsw=0                                                          
c     write (*,1016)                                                   
c     start iteration                                                  
 48   mcards=90                                                        
      if(maxit-niter) 49,51,51                                         
 49   continue                                                         
      write (*,1023)                                                  
      do 50 i=1,mesh,5                                                 
c      write (*,1024)i,x (i),ru3(i),ruinl1(i),rufnl1(i),    ruinl2(i),
c     1rufnl2(i),ru(i)                                             
 50   continue                                                        
 1023 format ('1i,x,ru3,ruinl1,rufnl1,ruinl2,rufnl2,ru '  )           
 1024 format (i8,f10.4,1p6e16.7)                                       
      go to 10                                                         
 51   do 52 i=1,mesh                                                   
      rscore(i)=0.0                                                    
 52   rsvale(i)=0.0                                                     
c     solve schroedinger equation for each state in turn. also calculate
c           core and valence charge density.                            
      do 62 m=1,ncspvs                                                  
      e=ee(m)                                                          
      nn=nnlz(m)/100                                                   
      lam =nnlz(m)/10 -10*nn                                           
      xl= lam                                                          
      call scheq(z,e,lam,nn,kkk,mesh,c,thresh)                         
      if(m-ncores)53,53,55                                            
 53   do 54 i=1,kkk                                                   
 54   rscore(i)=rscore(i)+wwnl(m)*snlo(i)**2                          
      go to 57                                                        
 55   do 56 i=1,kkk                                                   
 56   rsvale(i)=rsvale(i)+wwnl(m)*snlo(i)**2                          
 57   do 58 i=1,kkk                                                   
 58   snl(m,i)=snlo(i)                                                 
      k4=kkk+1                                                         
      if(k4-mesh)59,59,61                                              
 59   do 60 i=k4,mesh                                                  
 60   snl(m,i)=0                                                       
c        wave functions from kkk+1 to end of mesh are set to zero      
 61   nkkk(m)=kkk                                                      
      mcards=mcards+2+((kkk-1)/40)*8                                   
 62   ee(m)=e                                                          
c     calculate total charge density and atomic exchange potential     
      do 64 i=1,mesh                                                   
 63   rsatom(i) =rscore(i)+rsvale(i)                                   
 64   ruexch(i)=-6. *((3. *r(i)*rsatom(i))/315.82734 )**(1. /3. )      
c     calculate atomic coulomb potential                               
      a1=0.0                                                           
      asum=0.0                                                         
      b1=0.0                                                           
      bsum=0.0                                                         
      h=0.0025 *c                                                      
      i=1                                                              
      xi(1)=0.0                                                        
      xj(1)=0.0                                                        
      do 66 j=1,nblock                                                 
      do 65 k=1,40                                                    
      i=i+1                                                            
      a2=rsatom(i)*.5                                                  
      a1=a1+a2                                                        
      b2=a2/r(i)                                                      
      b1=b1+b2                                                        
      xi(i)=asum+a1*h                                                 
      xj(i)=bsum+b1*h                                                 
      a1=a1+a2                                                        
 65   b1=b1+b2                                                        
      asum=xi(i)                                                      
      bsum=xj(i)                                                      
      a1=a2                                                           
      b1=b2                                                           
 66   h=h+h                                                           
 67   continue                                                        
 68   do 69 i=1,mesh                                                  
      xi(i)=-2. *z+2. *(xi(i)+r(i)*(xj(mesh)-xj(i)))                  
      xj(i)=xi(i)+ruexch(i)                                           
 69   continue                                                        
 70   do 71 i=1,mesh                                                  
      ruinl1(i)=ruinl2(i)                                             
      rufnl1(i)=rufnl2(i)                                             
      ruinl2(i)=ru(i)                                                 
 71   rufnl2(i)=xj(i)                                                  
 72   niter=niter+1                                                    
      pdelta =delta                                                    
      delta =0.0                                                       
      do 76 i=1,limit                                                  
 73   snlo(i)=ru(i)-xj(i)                                              
 74   xi(i)= abs(snlo(i))                                              
      if (xi(i)-delta)76,76,75                                         
 75   delta=xi(i)                                                      
      idelta=i                                                         
 76   continue                                                         
      write (*,1008)niter,z,delta,idelta,x(idelta),icut,x(icut)        
c     test self-consistency of atomic potential.                        
 77   if(delta-tol) 123,78,78                                           
c     if scf criterion not satisfied, calculate next trial potential.   
 78   if(iprsw) 79,79,87                                                
 79   do 80 i=2,limit                                                   
 80   ru(i)=.5 *(ru(i)+xj(i))                                           
      if(mesh.le.limit) go to 86                                        
 81   ruzm=xj(mesh)                                                     
      if(ruzm.eq.xj(limit)) go to 82
      go to 84 
 82   do 83i=limit,mesh                                                 
      ratio=(i-limit)/(mesh-limit)                                      
 83   ru(i)=.5 *(1. -ratio)*ru(i)+.5 *(1. +ratio)*xj(i)                 
      limit=mesh                                                        
      go to 86                                                          
 84   ratio=(ruzm-ru(limit))/(ruzm-xj(limit))                           
      do 85 i=limit,mesh                                                
 85   ru(i)=ruzm-ratio*(ruzm-xj(i))                                     
      limit =mesh                                                       
 86   iprsw=ipratt                                                      
      go to 110                                                         
 87   continue                                                          
 88   if(nonmon) 79,79,89  
 89   if(pdelta-delta)90,90,91                                          
c     if delta is not monotonic decreasing four times, bypass           
c       pratt improvement scheme                                        
 90   nonmon =nonmon-1                                                  
      if(nonmon)79,79,91                                                
 91   alph=0.5                                                          
c     pratt improvement scheme                                          
      do 101 i=2,icut                                                   
      xnum(i)=ruinl1(i)*rufnl2(i)-ruinl2(i)*rufnl1(i)                   
      denm(i)=rufnl2(i)-rufnl1(i)-ruinl2(i)+ruinl1(i)                   
 92   if( abs(denm(i)/ruinl2(i))-0.0001 ) 93,93,95                      
 93   continue                                                          
 94   alph=0.5                                                          
      go to 100                                                         
 95   alph=(xnum(i)/denm(i)-rufnl2(i))/snlo(i)                          
      if(alph) 96,99,97                                                 
 96   alph=0.0                                                          
      go to 100                                                         
 97   if(.5 -alph) 98,99,99                                             
 98   alph=0.5                                                          
 99   continue                                                          
 100  xi(i)=alph                                                        
 101  continue                                                          
      iprsw=iprsw-1                                                     
      if(kut) 104,102,104                                               
 102  continue                                                          
      ic=icut+20                                                        
      ic1=icut+1                                                        
      adel=xi(icut)*.05                                                 
      do 103 i=ic1,ic                                                   
      xi(i)=xi(i-1)-adel                                                
 103  xj(i)=xi(i)                                                       
 104  continue                                                          
      xj(1)=0.5                                                         
      xj(2)=xi(2)                                                       
      asum=xi(2)+xi(3)+xi(4)+xi(5)                                      
      do 105 i=3,icut                                                   
      xj(i)=asum*0.2                                                    
 105  asum=asum-xi(i-2)+xi(i+3)                                         
      if(kut) 108,106,108                                               
 106  continue                                                          
      ic1=ic+1                                                          
      do 107 i=ic1,mesh                                                 
      xj(i)=0.0                                                         
 107  ru(i)=rufnl2(i)                                                   
 108  continue                                                          
      do 109 i=2,ic                                                     
 109  ru(i)= rufnl2(i)+xj(i)*snlo(i)                                    
 110  continue                                                          
      if(kut) 111,113,111                                               
 111  icut=mesh                                                         
      limit=mesh                                                        
      do 112 i=2,mesh                                                   
      vlast= v(i)                                                       
      v(i)=ru(i)/r(i)                                                   
 112  xi(i)=v(i)-vlast                                                  
      go to 119                                                         
 113  continue                                                          
      icut =0                                                           
      do 118 i=2,mesh                                                   
      vlast=v(i)                                                        
      if(icut) 114,114,116                                              
 114  if(twozzz+ru(i)) 117,117,115                                      
 115  icut=i                                                            
 116  v(i)=-twozzz/r(i)                                                 
      go to 118                                                         
  117 v(i)=ru(i)/r(i)                                                   
 118  xi(i)=v(i)-vlast                                                  
 119  continue                                                          
      xi(1)=0.0                                                         
c     next trial eigenvalues predicted by perturbation theory           
      ncards=90                                                         
      do 122 m=1,ncspvs                                                 
      k=(nkkk(m)-1)/40                                                  
      h=0.0025 *c                                                       
      asum=0.0                                                          
      a1=0.0                                                            
      i=1                                                               
      do 121 j=1,k                                                      
      do 120 l=1,40                                                     
      i=i+1                                                             
      a2=xi(i)*snl(m,i)**2                                              
 120  a1=a1+a2*h                                                        
      asum=asum+a1-(a2*.5 )*h                                           
      h=h+h                                                             
 121  a1=(a2*.5 )*h                                                     
      ee(m)=ee(m)+asum                                                  
 122  ncards=ncards+8*k+2                                               
      go to 48                                                          
c     results transferred from internal memory to output tape(s).       
c     output is on tape b-5 for offline punching (bcd), and is          
c       copied to normal output tape 6 (a3)                             
 123  continue                                                          
 124  continue                                                          
      nc=1                                                              
      write (* ,1018)z,ncores,nvales,ion,iz,nc                          
      do 127 i=1,441                                                    
      if(twoion+ruinl2(i)) 127,127,125                                  
 125  do 126 m=i,441                                                    
 126  ruinl2(m)=-twoion                                                 
      go to 128                                                         
 127  continue                                                          
c   line 403
 128  continue                                                          
      do 129 min=1,440,5                                                
      max= min+4                                                        
      nc=nc+1                                                           
c      write (* ,1019)(ruinl2(m),m=min,max),iz,nc                        
 129  continue                                                          
      nc=nc+1                                                           
c      write ( *,1020)ruinl2(441),iz,nc                                  
c
c     loop over the orbitals
c
      do 132 m=1,ncspvs                                                 
         nlz=nnlz(m)                                                       
         kkk=nkkk(m)                                                       
         xl=nlz/10-10*(nlz/100)                                            
         lp=xl+1.0                                                         
c        compute first term of series. (snl(r)/r**(lam+1) at r=0)          
         do 130 i=1,4                                                      
           a(i,1)=1.0                                                        
           a(i,2)=r(i+1)                                                     
           a(i,3)=r(i+1)*r(i+1)                                              
           a(i,4)=r(i+1)*a(i,3)                                              
 130     a(i,5)=snl (m,i+1)/r(i+1)**lp                                     
         call crosym (4)                                                   
         nc=nc+1                                                           
c        writing IDENT
         write (* ,1017)nlz,ee(m),xl,wwnl(m),kkk,a(1,5)              
c
c        OPENING FILE 9 for wavefunction output
c
         write(letters,'(A4)') outname
         write(wnlz,'(i3)') nlz 
         write(anlz,5332)
5332     format('.h-s')
         outfil = letters(1:4) // wnlz(1:3) // anlz(1:4)
         istring=istring+1
         write(allnam(istring),5331) outfil
5331     format('plot "',A11,'" using 1:2 with lines')
         open(9,file=outfil)
         write (9 ,2018)z,ncores,nvales,ion,iz,nc
         write (9 ,1017)nlz,ee(m),xl,wwnl(m),kkk,a(1,5)              
         k1=kkk-1                                                          
           do 5137 kill=1,LimMesh
5137       write (9, *) r(kill), snl(m,kill)
         close(9)
           do 131 min=1,k1 ,5                                                 
             nc=nc+1                                                           
             max= min+4                                                        
c            write (* ,1019)(snl(m,i),i=min,max),iz,nc                         
 131       continue                                                          
         nc=nc+1                                                           
c        write (* ,1020) snl(m,kkk),iz,nc
c
c                   loop over the orbitals  132  ending
c                                   
 132  continue
c
c 
      write(*,'(A)') 'Finished the loop'                                                          
      ax=2.*z  
      v(1)=1.  
      ruinl2(1)=1. 
      do 8129 m=1,441
      vout(m)=v(m)
      v(m)=-v(m)*r(m)/ax
 8129 ruinl2(m)=-ruinl2(m)/ax
c 
c     Potential output
c
c     the name-construction is copied, not optimized
c
      write(letters,'(A4)') outname
      write(wnlz,'(A4)') '-pot' 
      write(anlz,5332)
!5332  format('.h-s')
      outfil = letters(1:4) // wnlz(1:4) // anlz(1:4)
      istring=istring+1
      write(allnam(istring),5339) outfil

5339  format('plot "',A12,'" using 1:2 with lines')
      open(9,file=outfil)
      write (9 ,1995)
      write (9 ,2018)z,ncores,nvales,ion,iz,nc      
1995  format('# Herman-Skillman potential')                   
c     
c     plot of the energies
c
      do 1332 m=1,ncspvs    
         yy=ee(m)
         xx=0.0
          write(9,*) xx,yy
         xx=20.0
          write(9,*) xx,yy
          write(9,*)
1332   continue                                                                   
       do 5167 kill=1,440
5167   write (9, *) r(kill), vout(kill)
      close(9)
c 
c     Potential output finished
c      
      write(*,6699)
 6699 format('    output to file 12')  
c      write(12,8801) (ruinl2(m),m=1,440,4)
c      write(15,1010) (ru(m),m=1,441) 
c      write(12,8801) (r(m),m=1,440,4) 
c      write(12,8801) (v(m),m=1,440,4) 
      write(6,6699)
 8801 format(f9.5,9f8.5)
      do 135 i=1,ncspvs                                                 
      dlx=0.0025 /3.0                                                   
      rex=0.                                                            
      l=1                                                               
      do 134 j=1,nblock                                                 
      lx=l+40                                                           
      lx2=lx-1                                                          
      l2=l+1                                                            
      w=4.0                                                             
      sum=x(l)*snl(i,l)*snl(i,l)                                        
      do 133 k=l2,lx2                                                   
      sum =sum+w*x(k)*snl(i,k)*snl(i,k)                                 
 133  w=6.0 -w                                                          
      rex=dlx*(sum+x(lx)*snl(i,lx)**2)+rex                              
      l=lx                                                              
      dlx=dlx+dlx                                                       
 134  continue                                                          
      rex=rex*c*c                                                       
 135  continue                                                          
      if(key-1)  136,137,139                                            
 136  key=1                                                             
 137  do 138 i=1,mesh                                                   
 138  ru3(i)=ru(i)                                                      
      ze3=z                                                             
      go to 10                                                          
 139  do 140 i=1,mesh                                                   
      ru2(i)=ru3(i)                                                     
 140  ru3(i)=ru(i)                                                      
      ze2=ze3                                                           
      ze3=z                                                             
      go to 10                                                          
  141 rewind 1 
      do 4433 kkk=1,istring
4433  write(*,'(A80)') allnam(kkk)
      stop 
      end                                                               
      subroutine scheq(zz,en,lambda,nofl,kkk,mess,scf,thresh)           
c     subroutine scheq                                                  
c     compute energy eigenvalue and wave function                       
c      originally  written by sherwood skillman                         
c      rca laboratories, princeton, new jersey, spring 1961             
c      modified by frank herman, summer 1961                            
c      further modified by richard kortum,  lockheed research           
c       laboratories, palo alto, california,  summer 1962               
      common LimMesh,xnum(521),rsvale(521),rsatom(521),v(521),vout(521)                 
      common snlo(521)  ,r(521)     ,ru(521)    ,ruexch(521),xi(521)    
      common xj(521)    ,snl(22,521),ruinl1(521),rufnl1(521),ruinl2(521)
      common rufnl2(521),ru2(521)   ,ru3(521)   ,nnlz(24)   ,wwnl(24)   
      common nkkk(24)   ,ee(24)     ,a(4,5)                             
      dimension p(5),q(5),t(5),d(5)                                     
      dimension x(521),rscore(521),denm(521),qq(521)                    
      equivalence (xnum,rscore),(rsvale,denm),(rsatom,qq)               
c     set up constants and initialize                                   
      z=zz                                                              
      lam=lambda                                                        
      nn=nofl                                                           
      mesh=mess                                                         
      c=scf                                                             
      many = 200                                                        
 1    e=en                                                              
      more =0                                                           
      less =0                                                           
      morev=0                                                           
      lessv=0                                                           
      nprint=0                                                          
      emore=0.                                                          
      eless =0.                                                         
      de=0.                                                             
      lamm=lam-1                                                        
      lamp=lam+1                                                        
      xlp=lamp                                                          
      ndcr=nn-lamp                                                      
      b=lam*lamp                                                        
      oc=r(2)                                                           
      h=oc                                                              
      hsq=h*h                                                           
      b3=(v(3)-v(2))/h-z/hsq                                            
      y=h+h                                                             
      flps=4*lam+6                                                      
      slpt=6*lam+12                                                     
      elpt=8*lam+20                                                     
      a1=-z/xlp                                                         
      ysq=y*y                                                           
      b1=-z-z                                                           
      ab1=a1*b1                                                         
      ab3=a1*b3                                                         
c     raise h and y to lam+1                                            
      htl=h                                                             
      ytl=y                                                             
      if(lam)7,4,2                                                      
 2    do 3 i=1,lam                                                      
      htl=htl*h                                                         
 3    ytl=ytl*y                                                         
 4    h1=hsq                                                            
      bohs=b/hsq                                                        
      boh=b1/h                                                          
      bth=b3*h                                                          
      bq3=bohs+boh+bth                                                  
      bq4=bohs/4. +boh/2. +bth+bth                                      
      epl=8+lam                                                         
      fpl=5+lam                                                         
      xifc=c*.21701389e-4  
c     start outward integration                                         
 5    nprint=nprint+1                                                   
      eps =e-eg                                                         
      eg =e                                                             
      if(many-nprint)  6,9,9                                            
 6    write (*,1001)nn,lam ,z                                           
 1001 format ('   no convergence on',i4,i1,f4.0)                       
      return                                                            
 7    nstop=77                                                          
 8    write (*,1002)nstop                                               
 1002 format('0stop',i4,8hin scheq)                                     
 9    do 10 i=1,mesh                                                    
 10   snlo(i)=0.0                                                       
      if(nprint-1) 7,11,19                                              
 11   continue                                                          
      do 12  i= 4,mesh                                                  
      qq(i) = v(i)+b/(r(i)*r(i))-e                                      
 12   continue                                                          
 13   m= mesh                                                           
      do 15 i=4,mesh                                                    
      if(qq(m)) 14,15,15                                                
 14   ik=m+1                                                            
      go to  17                                                         
 15   m=m-1                                                             
 16   nstop =521                                                        
c     q is everywhere positive                                          
      go to 8                                                           
 17   if(mesh-ik) 18,18,21                                              
 18   eps = qq(mesh-40)                                                 
      e = e+eps                                                         
 19   continue                                                          
      do 20 i=4,mesh                                                    
 20   qq(i) = qq(i)-eps                                                 
      go to 13                                                          
 21   continue                                                          
 22   ncross=0                                                          
      sign=1.                                                           
      h=oc                                                              
      y=h+h                                                             
c     b= lam*(lam+1)                                                    
c     b1= -2.d*z                                                        
      b2=3. *z/h-e+2. *v(2)-v(3)                                        
c     b3=(v(3)-v(2))/h -z/hsq                                           
c     a1= -z/(lam+1)                                                    
      a2=(ab1+b2)/flps                                                  
c     a2=(a1*b1+b2)/(4*lam+6)                                           
      a3=(a2*b1+a1*b2+b3)/slpt                                          
c     a3=(a2*b1+a1*b2+b3)/(6*lam+12)                                    
      a4=(a3*b1+a2*b2+ab3)/elpt                                         
c     a4=(a3*b1+a2*b2+a1*b3)/(8*lam+20)                                 
      p(3)=(1. +h*(a1+h*(a2+h*(a3+h*a4))))*htl                          
c     p(3)=(1.d+a1*h+a2*h**2+a3*h**3+a4*h**4)*h**(xl+1.d)               
      p(4)=(1. +y*(a1+y*(a2+y*(a3+y*a4))))*ytl                          
c     p(4)=(1.d+a1*y+a2*y**2+a3*y**3+a4*y**4)*y**(xl+1.d)               
      q(3)=bq3+b2                                                       
c     q(3)=(b+b1*h+b2*h**2+b3*h**3)/h**2                                
      q(4)=bq4+b2                                                       
c     q(4)=(b+b1*y+b2*y**2+b3*y**3)/y**2                                
      snlo(2)=p(3)                                                      
      snlo(3)=p(4)                                                      
      i=3                                                               
      dx=oc                                                             
      h1=h**2                                                           
      h2=h1/12.                                                         
      t(3)=p(3)*(1. -h2*q(3))                                           
      t(4)=p(4)*(1. -h2*q(4))                                           
      d(4)=t(4)-t(3)                                                    
      ncount=3                                                          
      nint=2                                                            
 23   i=i+1                                                             
c     if end of mesh is reached, modify trial eigenvalue                
      if(i-mesh) 25,24,24                                               
 24   if(ndcr-ncross) 40,47,47                                          
c     return to beginning of outward integration if necessary           
 25   q(5) =qq(i)                                                       
      if(ik-i)  37,37,26                                                
 26   d(5)=d(4)+h1*q(4)*p(4)                                            
      t(5)=d(5)+t(4)                                                    
      if(1. - abs(h2*q(5)))24,24,27                                     
 27   p(5)=t(5)/(1. -h2*q(5))                                           
      snlo(i) = p(5)                                                    
      if(sign) 28,7,29                                                  
 28   if(p(5))  31,31,30                                                
 29   if(p(5)) 30,31,31                                                 
 30   ncross=ncross+1                                                   
c     count changes in sign                                             
      sign=-sign                                                        
 31   ncount=ncount+1                                                   
      if(7-ncount)7,32,33                                               
 32   ncount=2                                                          
 33   nint=nint+1                                                       
      if(40-nint)7,34,35                                                
 34   dx=dx+dx                                                          
      h=dx                                                              
      h1=h**2                                                           
      h2=h1/12.                                                         
      nint=0                                                            
      t(5)=p(5)*(1. -h2*q(5))                                           
      t(3)=p(3)*(1. -h2*q(3))                                           
      d(5)=t(5)-t(3)                                                    
 35   do 36 k=1,4                                                       
      p(k)=p(k+1)                                                       
      t(k)=t(k+1)                                                       
      d(k)=d(k+1)                                                       
 36   q(k)=q(k+1)                                                       
      go to 23                                                          
 37   if(ncount-2)7,38,26                                               
 38   if(nint-4)26,26,39                                                
c     matching radius has been reached going out                        
c     if ndcr not equal to ncross, modify trial eigenvalue              
 39   eigen=e                                                           
      if(ndcr-ncross) 40,55,47                                          
 40   more=1                                                            
c     too many crossings, increase absf(e)                              
      morev=morev+1                                                     
      if(morev-1) 41,43,42                                              
 41   nstop=50                                                          
      go to 8                                                           
 42   if(e -emore) 43,44,44                                             
 43   emore=e                                                           
 44   if (less) 45,46,54                                                
 45   nstop=55                                                          
      go to 8                                                           
 46   e=1.25 *eg                                                        
      go to 5                                                           
 47   less=1                                                            
c     too few  crossings, decrease absf(e)                              
      lessv=lessv+1                                                     
      if(lessv-1) 48,50,49                                              
 48   nstop=57                                                          
      go to 8                                                           
 49   if(eless- e) 50,51,51                                             
 50   eless=e                                                           
 51   if(more) 52,53,54                                                 
 52   nstop=62                                                          
      go to 8                                                           
 53   e=0.75 *eg                                                        
      go to 5                                                           
 54   e=.5 *(emore+eless)                                               
      go to 5                                                           
 55   if( abs(snlo(i-1))- abs(snlo(i-2))) 56,59,59                      
c     check to see that wave is in the damped region (absolute value    
c       decreasing and signs alike)                                     
 56   if (p(5)) 57,26 ,58                                               
 57   if(snlo(i-2)) 60,26,26                                            
 58   if(snlo(i-2)) 26,26,60                                            
 59   if(1.0e+25- abs(p(5)))47,47,26                                    
c     large absolute value of p in what should be the damped region     
c       indicates too few peaks, decrease absf(e)                       
c     now ndcr = ncross and matching radius lies in damped region       
 60   imatch=i-2                                                        
      xmatch=r(i-2)                                                     
c line 704
      ppout=(t(4)-t(2)-.5 *(p(4)-p(2)))/h                               
      s2=ppout/p(3)                                                     
c     integration is by 8 applications of newton-cotes closed           
c       quadrature for five intervals on each block                     
c       xifc =(5*h(block 1)/288)/2 ,h(1) =0.0025d*scale factor          
      sum1=0.0                                                          
      xif=xifc                                                          
      i=1                                                               
      value=0.0                                                         
 61   mm=8                                                              
      sum2=0.0                                                          
      xif=xif+xif                                                       
 62   y=value                                                           
      sum2=sum2+19. *(value+y)+75. *(snlo(i+4)**2+snlo(i+1)**2)         
     1+50. *(snlo(i+2)**2+snlo(i+3)**2)                                 
      i=i+5                                                             
      if (imatch-i)  7,65,63                                            
 63   mm=mm-1                                                           
      if(mm)7,64,62                                                     
 64   sum1=sum2*xif+sum1                                                
      go to 61                                                          
 65   sum1= sum1+sum2*xif                                               
 66   s1=sum1/p(3)**2                                                   
      pmatch=p(3)                                                       
      if(nn-1)7,67,68                                                   
 67   xinw=epl*xmatch                                                   
c     for n =1, start inward integration at(8+lam)*xmatch or x max      
      go to 69                                                          
 68   xinw=fpl*xmatch                                                   
c     for n not=1, start at (5+lam)*xmatch or x max (end of mesh)       
 69   do 71  i= 41,mesh,40                                              
      if(xinw-r(i))  70,70,71                                           
 70   kkk =i                                                            
      go to 72                                                          
 71   continue                                                          
      kkk =mesh                                                         
 72   i =kkk                                                            
      dx =r(i-1)-r(i)                                                   
      h =dx                                                             
      xif=.17361111e-01*dx                                              
      hsq=h*h                                                           
      hsq12=hsq/12.                                                     
      q(3)= qq(i)                                                       
      p(3)= exp(-r(i)* sqrt(q(3)))                                      
 73   sum3=p(3)/q(3)                                                    
      i=i-1                                                             
      q(4) =qq(i)                                                       
 74   p(4)= exp(-r(i)* sqrt(q(4)))                                      
      if( abs(p(4))-1.0e-35) 75,75,77                                   
 75   kkk=kkk -40                                                       
      if(kkk -imatch) 76,76,72                                          
 76   write (*,1003)z ,nn,lam,kkk                                       
 1003 format  (/'at z=',f6.0,'  nl =',i3,i1,' kkk =',i5,' is less tha',
     1'n imatch =',i5,/' inward integration will be tried at kkk+40')     
      kkk =kkk +40                                                      
      p(4 ) = 1.5e-35                                                   
      p(3) = 1.0e-35                                                    
 77   if(pmatch)78,7,79                                                 
 78   p(3)=-p(3)                                                        
      p(4)=-p(4)                                                        
 79   snlo(i+1)=p(3)                                                    
      snlo(i)=p(4)                                                      
      t(3)=p(3)*(1. -hsq12 *q(3))                                       
      t(4)=p(4)*(1. -hsq12*q(4))                                        
      d(4)=t(4)-t(3)                                                    
 80   do 82 m=2,40                                                      
      i=i-1                                                             
      q(5) =qq(i)                                          
      d(5)=hsq*q(4)*p(4)+d(4)                                
      t(5)=d(5)+t(4)                                                    
      p(5)=t(5)/(1. -hsq12*q(5))                           
      if(i-imatch+1)7,84,81                                
 81   snlo(i)=p(5)                                           
      do 82 k=1,4                                                       
      p(k)=p(k+1)                                                       
      t(k)=t(k+1)                                                       
      d(k)=d(k+1)                                                       
 82   q(k)=q(k+1)                                                       
      q(5) =qq(i-2)                                                     
      d(5)=hsq*q(4)*p(4)+d(4)                                           
      t(5)=d(5)+t(4)                                                    
      p(5)=t(5)/(1. -hsq12*q(5))                                        
      p(5)=1.09375 *p(4)+0.2734375 *p(5)-0.546875 *p(3)+0.21875 *p(2)-  
     1  0.0390625 *p(1)                                                   
      i=i-1                                                             
      dx=dx*.5                                                          
      q(5) =qq(i)                                                       
      h=dx                                                              
      hsq=h*h                                                           
      hsq12=hsq/12.                                                     
      t(5)=p(5)*(1. -hsq12*q(5))                                        
      t(4)=p(4)*(1. -hsq12*q(4))                                        
      d(5)=t(5)-t(4)                                                    
      snlo(i)=p(5)                                                      
      do 83 l=1,4                                                       
      p(l)=p(l+1)                                                       
      t(l)=t(l+1)                                                       
      d(l)=d(l+1)                                                       
 83   q(l)=q(l+1)                                                       
      go to 80                                                          
c     matching radius has been reached coming in                        
 84   k=kkk                                                             
      value=snlo(k)**2                                                  
      go to 87                                                          
 85   continue                                                          
 86   sum3=sum3+xif*sum4                                                
      xif =xif*0.5                                                      
 87   mm=8                                                              
      sum4 =0.0                                                         
 88   y=value                                                           
      value=snlo(k-5)**2                                                
      sum4=sum4+19. *(value+y)+75. *(snlo(k-1)**2+snlo(k-4)**2)         
     1+50. *(snlo(k-2)**2+snlo(k-3)**2)                                 
      k=k-5                                                             
      if(k-imatch) 7,90,89                                              
 89   mm=mm-1                                                           
      if(mm)  7,85,88                                                   
 90   sum3=sum3+xif*sum4                                                
 91   s3=sum3/p(4)**2                                                   
      ppin=(t(5)-t(3)-.5 *(p(5)-p(3)))/h                                
      s4=ppin/p(4)                                                      
      de=(s2-s4)/(s1-s3)                                                
      if(abs(de/e)-thresh) 94,92,92                                     
   92 e=e+de                                                            
      if(e) 5,93,93                                                     
   93 e=e-de                                                            
      de=de*.5                                                          
      go to 92                                                          
c     improve trial eigenvalue by perturbation theory if necessary      
c     calculate the normalized wave functions                           
 94   pop=pmatch/p(4)                                                   
      do 95 j=imatch,kkk                                                
 95   snlo(j)=snlo(j)*pop                                               
      sum1=0.0                                                          
      j=1                                                               
      xif=xifc                                                          
      value=0.0                                                         
 96   mm=8                                                              
      xif=xif+xif                                                       
      sum2=0.0                                                          
 97   y=value                                                           
      value=snlo(j+5)**2                                                
      sum2=sum2+19. *(value+y)+75. *(snlo(j+4)**2+snlo(j+1)**2)         
     1 +50. *(snlo(j+2)**2+snlo(j+3)**2)                                
      j=j+5                                                             
      mm=mm-1                                                           
      if(mm)7,98,97                                                     
 98   sum1=sum1+xif*sum2                                                
      if(kkk-j)7,99,96                                                  
 99   c1= sqrt(sum1)                                                    
      if(snlo(3))100,7,101                                              
 100  c1=-c1                                                            
 101  do 102 i=1,kkk                                                    
 102  snlo(i)=snlo(i)/c1                                                
      en=e                                                              
      return                                                            
      end                                                               
      subroutine crosym(m)                                              
c     simultaneous equation solver                                      
c     written by i.c. hanson, scientific computation department,        
c     lockheed missles and space company, sunnyvale, california         
c     solve m simultaneous equations by the method of crout             
      common LimMesh,xnum(521),rsvale(521),rsatom(521),v(521),vout(521)                 
      common snlo(521)  ,r(521)     ,ru(521)    ,ruexch(521),xi(521)    
      common xj(521)    ,snl(22,521),ruinl1(521),rufnl1(521),ruinl2(521)
      common rufnl2(521),ru2(521)   ,ru3(521)   ,nnlz(24)   ,wwnl(24)   
      common nkkk(24)   ,ee(24)     ,a(4,5)                             
      dimension x(521),rscore(521),denm(521),qq(521)                    
      equivalence (xnum,rscore),(rsvale,denm),(rsatom,qq)               
      n=m+1                                                             
      i1=1                                                              
 1    i3=i1                                                             
      sum= abs(a(i1,i1))                                                
      do3i=i1,m                                                         
      if(sum- abs(a(i,i1)))2,3,3                                        
 2    i3=i                                                              
      sum= abs(a(i,i1))                                                 
 3    continue                                                          
      if(i3-i1)4,6,4                                                    
 4    do 5j=1,n                                                         
      sum=-a(i1,j)                                                      
      a(i1,j)=a(i3,j)                                                   
 5    a(i3,j)=sum                                                       
 6    i3=i1+1                                                           
      do7i=i3,m                                                         
 7    a(i,i1)=a(i,i1)/a(i1,i1)                                          
 8    j2=i1-1                                                           
      i3=i1+1                                                           
      if(j2)9,11,9                                                      
 9    do10j=i3,n                                                        
      do10i=1,j2                                                        
 10   a(i1,j)=a(i1,j)-a(i1,i)*a(i,j)                                    
      if(i1-m)11,13,11                                                  
 11   j2=i1                                                             
      i1=i1+1                                                           
      do12i=i1,m                                                        
      do12j=1,j2                                                        
 12   a(i,i1)=a(i,i1)-a(i,j)*a(j,i1)                                    
      if(i1-m)1,8,1                                                     
 13   do15i=1,m                                                         
      j2=m-i                                                            
      i3=j2+1                                                           
      a(i3,n)=a(i3,n)/a(i3,i3)                                          
      if(j2)14,16,14                                                    
 14   do15j=1,j2                                                        
 15   a(j,n)=a(j,n)-a(i3,n)*a(j,i3)                                     
 16    return                                                           
       end

