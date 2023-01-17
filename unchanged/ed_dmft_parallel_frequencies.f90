program lisalanc
    use mpi
    use types
    use AuxRoutines
    use LatticeRoutines
    implicit none
    !TODO: search and remove unused variables
    !TODO: move global variables from here to subroutines 
    !TODO: blocks of definitions for: model,simulation,auxilliary,... variables

    integer blocksize1,i,j,k
    integer imaxmu,iteramax,imu,itera,Nitermax,ikx,iky,ikz
    integer omnumber,omrest,omoffset
    real(dp) deltamu,testdmft,piob,xmu0,chi2,chi22,kx,ky,kz
    include'init.h'  
    integer,  parameter :: nsp1 = ns+1 
    integer,  parameter :: nblock = nsp1*nsp1
    integer,  parameter :: ntot = 4**ns
    integer,  parameter :: hubb_dat_fh = 30
    integer,  parameter :: hubb_andpar_fh = 15
    real(dp), parameter :: Pi = 4 * atan(1.d0)
    real(dp) zpart,t,t1,t2
    real(dp) tpar(ns),epsk(ns),uhub,hmag,xmu,xmuold,beta,densimp
    real(dp) wlow,wup,deltino,range,sumdos,help1,help2,dens,sum
    real(dp) om(0:Iwmax), dos(0:Iwmaxreal),dospart(0:Iwmaxreal)
    real(dp) omm
    complex(dp) ni(-10*Iwmax:10*Iwmax)
    integer offsetnum(0:prozessoren-1),recvnumber(0:prozessoren-1)
    integer offsetnuml(0:prozessoren-1),recvnumberl(0:prozessoren-1)
    integer offsetnume(0:prozessoren-1),recvnumbere(0:prozessoren-1)
    integer lehmnumber,eigennumber
      
    complex(dp) omr(0:Iwmaxreal)
    complex(dp) Gw(0:Iwmax),self(0:Iwmax),Gwpart(0:Iwmax)
    complex(dp) Gww(0:2*Iwmax-1)
    complex(dp) Gwreal(0:Iwmaxreal),Gwrealpart(0:Iwmaxreal)
    complex(dp) sigre(0:Iwmaxreal)
    complex(dp) G0w(0:Iwmax),G0wand(0:Iwmax),G0wwand(0:2*Iwmax-1)
    complex(dp) G0wpart(0:Iwmax)
!   matrix & friends
    integer ioffset(nblock+1),b(ns),ntable(ntot),abas(ns)
    integer nleng(nblock),imat(nmax,nmax,nblock),idmat(ntot)
    integer iddomat(ntot),idouble(ntot)
!     for minimization
    real(dp) hess(nmpara**2),g(nmpara),xtemp(nmpara),w(nmpara**2),xprmt(nmpara)
!     for dspev
    real(dp) work(3*nmax)
    real(dp) eig(nmax),zeig(nmax,nmax)
    real(dp) eigtotal(nmax,nblock+1)
    real(dp) eigold(nmax),zold(nmax,nmax)
    real(dp) realmat(nmax*(nmax+1)/2)
    real(dp) xmat(nmax,nmax)
    real(dp) rlehm(nmax,nmax)
    real(dp) rlehmtotal(nmax,nmax,nblock+1)

    complex(dp) omega, cdummy1
    real(dp) Epskold(ns),tparold(ns)
    complex(dp) sig0
    real(dp) zed,dens_wanted,diffdens,chi,densold,emin,dostest,emin1
    real(dp) xmularge,xmusmall,denslarge,denssmall,densmix,emintot
    real(dp) doubletot,zparttot
    integer ifix,inew,iauto
    integer iattempt,idritto,ilarge,ismall
    integer imix,iteraok,iexp
    real(dp) threshold,th0,b1,b2,E
    real(dp) doublep(ntot)
    real(dp) double
      
    integer ierror,myid,nprocs,mpi_type,myblock
    integer sendstatus(MPI_STATUS_SIZE)
    integer recvstatus(MPI_STATUS_SIZE)
    integer request,next,previous
    real(dp) Energya(0:ksteps,0:ksteps,0:ksteps)
    real(dp) Energya2D(0:ksteps,0:ksteps)
     
    !TODO: t,t1,t2 should either be compile time parameters (like in init.h) or
    !....: run time parameters, like in hubb.andpar and read with a function
    include 'tpri.dat'
       
      
    call MPI_INIT(ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierror)
  
    myblock=myid+1
    next=myid+1
    previous=myid-1
      

    omnumber=(Iwmax+1)/nprocs
    omrest=mod(Iwmax+1,nprocs)
    do i=0,nprocs-1
         if (i.lt.omrest) then
            recvnumber(i)=omnumber+1
            offsetnum(i)=i*(omnumber+1)
         else
            recvnumber(i)=omnumber
            offsetnum(i)=omrest*(omnumber+1)+(i-omrest)*omnumber
         endif
      enddo
      if (myid.lt.omrest) then
         omnumber=omnumber+1
      endif
      if (myid.lt.omrest) then
         omoffset=myid*omnumber-1
      else
         omoffset=omrest*(omnumber+1)+(myid-omrest)*omnumber-1
      endif

      if (myid.eq.0) then
         if (lattice_type.eq.4) then
            write(6,*) 'BETHE lattice case'
         else
            if (lattice_type.eq.3) then
               write(6,*) '2D-simple cubic lattice'
            else
               if (lattice_type.eq.2) then   
                 write(6,*) '3D simple  cubic lattice'
               else
                 write(6,*) '3D face centered cubic lattice'
               endif
            endif
         endif
      endif 

      if (lattice_type.eq.1) then
         do ikx=0,ksteps
            kx=Pi*dfloat(ikx)/dfloat(ksteps)
            do iky=0,ksteps
               ky=Pi*dfloat(iky)/dfloat(ksteps)
               do ikz=0,ksteps
                  kz=Pi*dfloat(ikz)/dfloat(ksteps)
                  Energya(ikx,iky,ikz)=dispersion_sc3D(t,t1,t2,kx,ky,kz)
               enddo
            enddo
         enddo
      endif

      if (lattice_type.eq.2) then
         do ikx=0,ksteps
            kx=2*Pi*dfloat(ikx)/dfloat(ksteps)
            do iky=0,ksteps
               ky=2*Pi*dfloat(iky)/dfloat(ksteps)
               do ikz=0,ksteps
                  kz=2*Pi*dfloat(ikz)/dfloat(ksteps)
                  Energya(ikx,iky,ikz)=dispersion_fcc(t,t1,t2,kx,ky,kz)
               enddo
            enddo
         enddo
      endif

    if (lattice_type.eq.3) then
        do ikx=0,ksteps
            kx=2*Pi*dfloat(ikx)/dfloat(ksteps)
            do iky=0,ksteps
                ky=2*Pi*dfloat(iky)/dfloat(ksteps)
                Energya2D(ikx,iky)=dispersion_sc2D(t,t1,t2,kx,ky)
            enddo
        enddo
    endif

    if (myid.eq.0) then
        open(hubb_dat_fh,file='hubb.dat',form='formatted',status='old')
        open(hubb_andpar_fh,file='hubb.andpar',form='formatted',status='old')
         
        call datain(hubb_dat_fh,hubb_andpar_fh,ns,uhub,xmu,hmag,beta,&
                    wlow,wup,deltino,&
                imaxmu,deltamu,iteramax,testdmft,dens_wanted,ifix,   &
                inew,iexp,th0,iauto)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_BCAST(uhub,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(hmag,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(beta,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(wlow,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(wup,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(deltino,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(imaxmu,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(deltamu,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(iteramax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(testdmft,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(dens_wanted,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(ifix,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(inew,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(iexp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(th0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
    call MPI_BCAST(iauto,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)

    if (myid.eq.((ns+1)**2-1)) next=0
    if (myid.eq.0)             previous=(ns+1)**2-1
    if (myblock.gt.(ns+1)**2)  myblock=1
    !TODO: most of these can probably be eliminated
    piob=Pi/beta
    diffdens=10000.0d0
    double=0.0d0
    zpart=0.0d0
    dens=0.0d0
    emin1=10000.0d0

    if (myid.le.nblock-1) then
        lehmnumber=nmax**2
        eigennumber=nmax
    else 
        lehmnumber=1
        eigennumber=1
    endif
    do i=0,nprocs-1
        if (i.le.nblock-1) then
            recvnumberl(i)=nmax**2
            offsetnuml(i)=i*nmax**2
            recvnumbere(i)=nmax
            offsetnume(i)=i*nmax
        else
            recvnumberl(i)=1
            offsetnuml(i)=nblock*nmax**2+i-nblock
            recvnumbere(i)=nmax
            offsetnume(i)=nblock*nmax+i-nblock
        endif
    enddo
    rlehm = 0.d0
    zeig  = 0.d0
    tpar  = 0.d0
    epsk  = 0.d0
    b(1)=1
    do i=2,ns
        b(i)=b(i-1)*4
    enddo
    range=(wup-wlow)/dfloat(Iwmaxreal)
    do i=0,Iwmaxreal
        omr(i)=wlow+range*dfloat(i)+(0.d0,1.d0)*deltino
    enddo

    !TODO: IO should probably be moved to a function, since this msses up the code in several places
      if (myid.eq.0) then
         write(*,*) 'nmax = ', nmax
         write(*,*)'U  = ',uhub
         write(*,*)'mu = ',xmu
         write(*,*)'h = ',hmag
         write(*,*)'beta =',beta
         write(*,*) 'Iwmax =',Iwmax
         write(*,*)'for real frequencies: ',wlow,wup
         call flush(6) 
         if (ifix.eq.0) then
            write(6,*)'Working at fixed mu =',xmu
         elseif(ifix.eq.1) then
            write(6,*)'Working at fixed density =',dens_wanted  
            if(beta.ge.200.d0) then
               chi= 1.2d0/(1+0.6d0*dabs(uhub))
               xmu=(dens_wanted-1.d0)/chi
               xmu=0.5d0*uhub+xmu
            endif
         else
            write(6,*)'ifix should be 0, 1 !!'
            stop
         endif
         chi= 1.2d0/(1+0.6d0*dabs(uhub))
         call flush(6)
         write(6,*)'Initial mu =',xmu
         open(34,file='self-en_wim',form='formatted',status='unknown')
         open(90,file='gm_wim',form='formatted',status='unknown')
         open(91,file='gm_wre',form='formatted',status='unknown')
         open(92,file='g0m',form='formatted',status='unknown')
         open(93,file='g0mand',form='formatted',status='unknown')
         open(36,file='zpart.dat',form='formatted',status='unknown')
         open(37,file='doubleOcc.dat',form='formatted',status='unknown')
         open(38,file='densimp.dat',form='formatted',status='unknown')
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      
      call MPI_BCAST(chi,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
      call computematrix(ioffset,nblock,ntot,b,ns,abas,nsp1,nmax, &
          nleng,imat,idmat,iddomat,ntable,idouble,myid)
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      
      if (myblock.eq.1) then
         blocksize1=1
      else
         blocksize1=nleng(myblock-1)
      endif

      if (myid.eq.0) then
         write(6,*)'# of independent blocks:',nblock
      endif

    do i=0,Iwmax
        om(i)=(2.d0*dfloat(i)+1.d0)*piob
    enddo
    do i=-10*Iwmax,10*Iwmax
        ni(i)=dfloat(i)*piob*(0.d0,1.d0)
    enddo
    if (myid.eq.0) then
       call read_andpar(hubb_andpar_fh,ns,epsk,tpar,xmuold)
    !TODO: IO should probably be moved to a function, since this msses up the code in several places
       write(*,*)'------------------------------------------------'
       write(*,*) 'starting Anderson parameters '
       write(*,*)'------------------------------------------------'
       write(*,*) 'Eps(k) '
       do i=2,ns
           write(*,'(2f27.16)')epsk(i)
       enddo
       write(*,*)'V(k) '
       do i=1,ns-1
           write(*,'(2f27.16)')tpar(i)
       enddo
       if (ifix.eq.1.and.inew.eq.0) xmu=xmuold
       if (inew.eq.0.and.iauto.eq.0) xmu=xmuold
       write(*,'(f27.16,"   #chemical potential")')xmu
    endif

      iteraok=1
      
      do itera=1,iteramax
         if (myid.eq.0) then
            imix=0
            threshold=th0/dfloat(iteraok)**iexp
            if (threshold.lt.1.d-6) threshold=1.d-6
            epsk(1)=-xmu
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_BCAST(tpar,ns,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         call MPI_BCAST(epsk,ns,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         call MPI_BCAST(xmu,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

         do i=0,Iwmax
            call calcg0(om(i),cdummy1,tpar,epsk,ns,Iwmax)
            g0wand(i)=cdummy1
         end do
         if (myid.eq.0) then   
    !TODO: IO should probably be moved to a function, since this msses up the code in several places
            write(6,*)'------------------------------------------------'
            write(6,*)'   Iteration : ',itera
            write(6,*)'------------------------------------------------'
            write(6,*)'Threshold for n=',threshold
            if (ifix.ne.0) write(6,*)'initial mu = ',xmu
            call flush(6)
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_BCAST(threshold,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)

         iattempt=0 
         idritto=0
         ilarge=0
         ismall=0
         double=0.d0

         if (myid.lt.(ns+1)**2) then

         call diag(ns,nsp1,nblock,ntot,nmax,myblock,nleng(myblock), &
                    b,ntable,ioffset,imat(1,1,myblock),xmu,uhub, &
                    hmag,tpar,epsk,emin1,eig,zeig)                     
  
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_REDUCE(emin1,emintot,1,MPI_REAL8,MPI_MIN, &
                           0,MPI_COMM_WORLD,ierror)
         if (myid.eq.0) then
            emin=emintot
            write(6,*)emin
            call flush(6)
         endif

         call MPI_BCAST(emin,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         do i=1,nleng(myblock)
            eig(i)=eig(i)-emin
         enddo

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)

         if (myid.lt.(ns+1)**2) then

         call MPI_ISSEND(eig,nmax,MPI_REAL8,next, &
                         98,MPI_COMM_WORLD,request,ierror)
         call MPI_RECV(eigold,nmax,MPI_REAL8,previous, &
                         98,MPI_COMM_WORLD,recvstatus,ierror)
         call MPI_WAIT(request,sendstatus,ierror)
         call MPI_ISSEND(zeig,nmax**2,MPI_REAL8,next, &
                         99,MPI_COMM_WORLD,request,ierror)
         call MPI_RECV(zold,nmax**2,MPI_REAL8,previous, &
                         99,MPI_COMM_WORLD,recvstatus,ierror)
         call MPI_WAIT(request,sendstatus,ierror)

         call lehmann(ns,nsp1,nblock,ntot,nmax,myblock, &         
                     nleng(myblock),blocksize1, &
                     ioffset,idouble,idmat,     &               
                     xmu,uhub,hmag,beta,        &              
                     eig,zeig,eigold,zold,      &             
                     double,zpart,dens,         &            
                     rlehm)                                
         
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         
         call MPI_REDUCE(double,doubletot,1,MPI_REAL8,MPI_SUM, &
                        0,MPI_COMM_WORLD,ierror)
         call MPI_REDUCE(zpart,zparttot,1,MPI_REAL8,MPI_SUM, &
                        0,MPI_COMM_WORLD,ierror)
         call MPI_REDUCE(dens,densimp,1,MPI_REAL8,MPI_SUM, &
                        0,MPI_COMM_WORLD,ierror)
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         if (myid.eq.0) then
            doubletot=doubletot/zparttot
            densimp=densimp/zparttot
            if (ifix.eq.1) then
               diffdens=densimp-dens_wanted
            endif
         endif

         if (ifix.eq.1) then
            call MPI_BCAST(diffdens,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         endif
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         
         if (ifix.eq.0) then

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   Computation of the Green-Function on the imaginary axis           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            if (myid.eq.0) then
               write(6,*)'Green-Function started'
               call flush(6)
            endif

            call MPI_BARRIER(MPI_COMM_WORLD,ierror)
            call MPI_ALLGATHERV(rlehm,lehmnumber,MPI_REAL8, &
                               rlehmtotal,recvnumberl,offsetnuml, &
                               MPI_REAL8,MPI_COMM_WORLD,ierror)   
            call MPI_ALLGATHERV(eig,eigennumber,MPI_REAL8, &
                               eigtotal,recvnumbere,offsetnume, &
                               MPI_REAL8,MPI_COMM_WORLD,ierror)   
            call computegimag(omoffset,omnumber,  &    
                            nsp1,nmax,Iwmax,prozessoren,nleng, &                 
                            beta,eigtotal,rlehmtotal,om,Gwpart)             
            call MPI_BARRIER(MPI_COMM_WORLD,ierror)
            call MPI_GATHERV(Gwpart(omoffset+1:omoffset+omnumber), &
                            omnumber,MPI_DOUBLE_COMPLEX,Gw,recvnumber, &
                            offsetnum, MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierror)
            call MPI_BARRIER(MPI_COMM_WORLD,ierror)

            if (myid.eq.0) then
               do i=0,Iwmax
                  Gw(i)=Gw(i)/zparttot
               enddo
               write(6,*)'Green-Function finished'
               call flush(6)
            endif
            call MPI_BARRIER(MPI_COMM_WORLD,ierror)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   End Green-Function Calculation                                    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         endif
         if ((ifix.eq.0).and.(myid.eq.0)) then
            write(6,*)'densimp = ',densimp
         endif
 
         if (ifix.eq.1) then 
 67         continue
            if (myid.eq.0) then
               write(6,*)'<n>=',densimp,'  <n>-n0=',diffdens
            endif
            if (dabs(diffdens).lt.threshold) then               
               if (myid.eq.0) then
                  write(6,*)'------------------------'
                  write(6,*)'converged at step ',itera
                  write(6,*)'------------------------'
                  call flush(6)
               endif
!           converged now
               if (myid.eq.0) then
                  write(6,*) 'Computing G_imag'
                  call flush(6)
               endif

               

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   Computation of the Green-Function on the imaginary axis           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               if (myid.eq.0) then
                  write(6,*)'Green-Function started'
                  call flush(6)
               endif

               call MPI_BARRIER(MPI_COMM_WORLD,ierror)
               call MPI_ALLGATHERV(rlehm,lehmnumber,MPI_REAL8,                      &
                                  rlehmtotal,recvnumberl,offsetnuml,                &
                                  MPI_REAL8,MPI_COMM_WORLD,ierror)   
               call MPI_ALLGATHERV(eig,eigennumber,MPI_REAL8,                       &
                                  eigtotal,recvnumbere,offsetnume,                  &
                                  MPI_REAL8,MPI_COMM_WORLD,ierror)
               call computegimag(omoffset,omnumber,                                 & 
                                  nsp1,nmax,Iwmax,prozessoren,                      &        
                                  nleng,beta,eigtotal,rlehmtotal,om,                 &
                                  Gwpart)                   
               call MPI_BARRIER(MPI_COMM_WORLD,ierror) 
               call MPI_GATHERV(Gwpart(omoffset+1:omoffset+omnumber),               &
                                  omnumber,MPI_DOUBLE_COMPLEX,Gw,recvnumber,        &
                                  offsetnum,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,    &
                                  ierror)
               call MPI_BARRIER(MPI_COMM_WORLD,ierror)

               if (myid.eq.0) then
                  do i=0,Iwmax
                     Gw(i)=Gw(i)/zparttot
                  enddo
                  write(6,*)'Green-Function finished'
                  call flush(6)
               endif
               call MPI_BARRIER(MPI_COMM_WORLD,ierror)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   End Green-Function Calculation                                    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               iteraok=iteraok+1
               if (myid.eq.0) then
                  open(101)
                  do i=0,Iwmax
                     write(101,*)i,om(i),Gw(i)
                  enddo
               endif

            else
!           not converged now
    !TODO: IO should probably be moved to a function, since this msses up the code in several places
               if (myid.eq.0) then
                  write(6,*)'-Not yet converged-'
                  call flush(6)
               endif
               if (diffdens.gt.0.d0) then
                  ilarge=1
                  xmularge=xmu
                  denslarge=densimp
               endif
               if (diffdens.lt.0.d0) then
                  ismall=1
                  xmusmall=xmu
                  denssmall=densimp
               endif
               if (ilarge*ismall.eq.0) then
                  if (myid.eq.0) then
                     write(6,*)'Still lacking a bracket'
                  endif
                  xmu=xmu-chi*diffdens
                  if (myid.eq.0) then
                     write(6,*)'Try with xmu=',xmu
                  endif
                  idritto=idritto+1
                  if (idritto.eq.6) then
                     chi=chi*2.d0
                  endif

                  if (myid.lt.(ns+1)**2) then
                  
                  call diag(ns,nsp1,nblock,ntot,nmax,myblock,nleng(myblock), &
                            b,ntable,ioffset,imat(1,1,myblock),xmu,uhub, &
                            hmag,tpar,epsk,emin1,eig,zeig)
  
                  endif

    !TODO: all these mpi calls should be one line. some of them are proably not used anyway
                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)


                  call MPI_REDUCE(emin1,emintot,1,MPI_REAL8,MPI_MIN, &
                                    0,MPI_COMM_WORLD,ierror)
                  if (myid.eq.0) then
                     emin=emintot
                     write(6,*)emin
                     call flush(6)
                  endif

                  call MPI_BCAST(emin,1,MPI_REAL8,0, &
                                MPI_COMM_WORLD,ierror)
                  
                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                  do i=1,nleng(myblock)
                     eig(i)=eig(i)-emin
                  enddo

                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                  
                  if (myid.lt.(ns+1)**2) then

                  call MPI_ISSEND(eig,nmax,MPI_REAL8,next, &
                               98,MPI_COMM_WORLD,request,ierror)     
             call MPI_RECV(eigold,nmax,MPI_REAL8,previous, &
                               98,MPI_COMM_WORLD,recvstatus,ierror)
                  call MPI_WAIT(request,sendstatus,ierror)
                  call MPI_ISSEND(zeig,nmax**2,MPI_REAL8,next, &
                               99,MPI_COMM_WORLD,request,ierror)
                  call MPI_RECV(zold,nmax**2,MPI_REAL8,previous, &
                               99,MPI_COMM_WORLD,recvstatus,ierror)
                  call MPI_WAIT(request,sendstatus,ierror)
         
                  call lehmann(ns,nsp1,nblock,ntot,nmax,myblock, &
                              nleng(myblock),blocksize1, &
                              ioffset,idouble,idmat, &    
                              xmu,uhub,hmag,beta, &       
                              eig,zeig,eigold,zold, &         
                              double,zpart,dens, &      
                              rlehm)                      
         
                  endif

                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                  
                  call MPI_REDUCE(double,doubletot,1,MPI_REAL8,MPI_SUM, &
                      0,MPI_COMM_WORLD,ierror)
                  call MPI_REDUCE(zpart,zparttot,1,MPI_REAL8,MPI_SUM, &
                      0,MPI_COMM_WORLD,ierror)
                  call MPI_REDUCE(dens,densimp,1,MPI_REAL8,MPI_SUM, &
                      0,MPI_COMM_WORLD,ierror)

                  if (myid.eq.0) then
                     doubletot=doubletot/zparttot
                     densimp=densimp/zparttot
                     diffdens=densimp-dens_wanted
                  endif
                  
                  call MPI_BCAST(diffdens,1,MPI_REAL8,0, &
                                MPI_COMM_WORLD,ierror)
                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                  go to 67
               else
                  call MPI_BCAST(denssmall,1,MPI_REAL8,0, &
                                MPI_COMM_WORLD,ierror)
                  call MPI_BCAST(denslarge,1,MPI_REAL8,0, &
                                MPI_COMM_WORLD,ierror)
                  if (myid.eq.0) then
                     write(6,*)'Interpolation between:'
                     write(6,*)'mu, n (1):',xmusmall,denssmall
                     write(6,*)'mu, n (2):',xmularge,denslarge
                  endif
                  xmu=xmusmall+(xmularge-xmusmall)* &
                      (denssmall-dens_wanted)/ &
                      (denssmall-denslarge)
                  if (myid.eq.0) then
                     write(6,*)'interpolated xmu',xmu
                  endif

                  if (myid.lt.(ns+1)**2) then

                  call diag(ns,nsp1,nblock,ntot,nmax,myblock, nleng(myblock), &
                            b,ntable,ioffset,imat(1,1,myblock), &
                            xmu,uhub,hmag,tpar,epsk,emin1,eig,zeig)
  
                  endif

                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                  call MPI_REDUCE(emin1,emintot,1,MPI_REAL8,MPI_MIN, &
                                 0,MPI_COMM_WORLD,ierror)
                  if (myid.eq.0) then
                     emin=emintot
                     write(6,*)emin
                     call flush(6)
                  endif

                  call MPI_BCAST(emin,1,MPI_REAL8,0, &
                                MPI_COMM_WORLD,ierror)
                  
                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                  do i=1,nleng(myblock)
                     eig(i)=eig(i)-emin
                  enddo

                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

                  if (myid.lt.(ns+1)**2) then

                  call MPI_ISSEND(eig,nmax,MPI_REAL8,next, &
                               98,MPI_COMM_WORLD,request,ierror)
                  call MPI_RECV(eigold,nmax,MPI_REAL8,previous, &
                               98,MPI_COMM_WORLD,recvstatus,ierror)
                  call MPI_WAIT(request,sendstatus,ierror)
                  call MPI_ISSEND(zeig,nmax**2,MPI_REAL8,next, &
                               99,MPI_COMM_WORLD,request,ierror)
                  call MPI_RECV(zold,nmax**2,MPI_REAL8,previous, &
                               99,MPI_COMM_WORLD,recvstatus,ierror)
                  call MPI_WAIT(request,sendstatus,ierror)

                  call lehmann(ns,nsp1,nblock,ntot,nmax,myblock, &
                              nleng(myblock),blocksize1, &
                              ioffset,idouble,idmat, &
                              xmu,uhub,hmag,beta, &
                              eig,zeig,eigold,zold, &
                              double,zpart,dens,rlehm)

                  endif
         
                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                  
                  call MPI_REDUCE(double,doubletot,1,MPI_REAL8,MPI_SUM, &
                      0,MPI_COMM_WORLD,ierror)
                  call MPI_REDUCE(zpart,zparttot,1,MPI_REAL8,MPI_SUM, &
                      0,MPI_COMM_WORLD,ierror)
                  call MPI_REDUCE(dens,densimp,1,MPI_REAL8,MPI_SUM, &
                      0,MPI_COMM_WORLD,ierror)

                  if (myid.eq.0) then
                     doubletot=doubletot/zparttot
                     densimp=densimp/zparttot
                     diffdens=densimp-dens_wanted
                  endif
                 
                  call MPI_BCAST(diffdens,1,MPI_REAL8,0, &
                                MPI_COMM_WORLD,ierror)
                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                          
                  goto 67
               endif
            endif
         endif


         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_BCAST(Gw,Iwmax+1,MPI_DOUBLE_COMPLEX,0, &
                       MPI_COMM_WORLD,ierror)
         call MPI_BCAST(xmu,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
 68      if (myid.eq.0) then
            write(6,*) 'SELF_CONSISTENCY LOOP'
         endif

         if (lattice_type.eq.1) then
            write(6,*) 'for the 3D-cubic' 
             call selfconst3(omnumber,omoffset,Iwmax,ksteps,beta,xmu, &
                        om,Energya,g0wand,Gw,G0wpart)
         elseif (lattice_type.eq.2) then              
            write(6,*) 'for the FCC' 
             call selfconst_fcc(omnumber,omoffset,Iwmax,ksteps,beta,xmu, &
                        om,Energya,g0wand,Gw,G0wpart)
         elseif (lattice_type.eq.3) then              
            write(6,*) 'for the 2D-square' 
             call selfconst2(omnumber,omoffset,Iwmax,ksteps, &
                      beta,xmu,om,Energya2D,g0wand,Gw,G0wpart)
         else
             error stop          
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_GATHERV(G0wpart(omoffset+1:omoffset+omnumber),omnumber, &
                         MPI_DOUBLE_COMPLEX,G0w,recvnumber,offsetnum, &
                         MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierror)
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         if (myid.eq.0) then
            if (iauto.eq.0) then 
               write(6,*)'Average Double Occupancy :',doubletot
               goto 777
            endif
            
            Nitermax=4000
            write(6,*)'minimization'
            call flush(6)

            call search(difference,chi2,Nitermax,hess,g,xtemp,w,xprmt,nmpara, &
                tpar,epsk,ns,piob,Iwmax,G0w,G0wand,om,symm)
            call flush(6)
            call read_andpar(hubb_andpar_fh,ns,Epskold,tparold,xmuold)

            chi22=0.d0
            do i=1,ns-1
               chi22=chi22+(tpar(i)-tparold(i))**2
               chi22=chi22+(Epsk(i+1)-Epskold(i+1))**2
            enddo
            chi22=chi22/(2*ns-2)

            call write_andpar(hubb_andpar_fh,ns,Epsk,tpar,xmu,Iwmax)
            write(*,*)'------------------------------------------------'
            write(*,'(" convergence parameter :",e18.10)')chi22
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_BCAST(chi22,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)

         if (chi22.lt.testdmft) goto 777
      enddo
!========+=========+=========+=========+=========+=========+=========+=$
!                               OUTPUT
!========+=========+=========+=========+=========+=========+=========+=$
         
777      continue

      if (iauto.ne.0) then
            
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_BCAST(tpar,ns,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         call MPI_BCAST(epsk,ns,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         call MPI_BCAST(xmu,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         
         if (myid.lt.(ns+1)**2) then

         call diag(ns,nsp1,nblock,ntot,nmax,myblock,nleng(myblock), &
                    b,ntable,ioffset,imat(1,1,myblock),xmu,uhub,hmag, &
                    tpar,epsk,emin1,eig,zeig)
  
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_REDUCE(emin1,emintot,1,MPI_REAL8,MPI_MIN, &
                        0,MPI_COMM_WORLD,ierror)
         if (myid.eq.0) then
            emin=emintot
            write(6,*)emin
            call flush(6)
         endif

         call MPI_BCAST(emin,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         do i=1,nleng(myblock)
            eig(i)=eig(i)-emin
         enddo

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)

         if (myid.lt.(ns+1)**2) then

         call MPI_ISSEND(eig,nmax,MPI_REAL8,next, &
                      98,MPI_COMM_WORLD,request,ierror)
         call MPI_RECV(eigold,nmax,MPI_REAL8,previous, &
                      98,MPI_COMM_WORLD,recvstatus,ierror)
         call MPI_WAIT(request,sendstatus,ierror)
         call MPI_ISSEND(zeig,nmax**2,MPI_REAL8,next, &
                      99,MPI_COMM_WORLD,request,ierror)
         call MPI_RECV(zold,nmax**2,MPI_REAL8,previous, &
                      99,MPI_COMM_WORLD,recvstatus,ierror)
         call MPI_WAIT(request,sendstatus,ierror)

         
         call lehmann(ns,nsp1,nblock,ntot,nmax,myblock, &
                     nleng(myblock),blocksize1, &
                     ioffset,idouble,idmat, &
                     xmu,uhub,hmag,beta, &
                     eig,zeig,eigold,zold, &
                     double,zpart,dens, &
                     rlehm)

         endif
         
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         
         call MPI_REDUCE(double,doubletot,1,MPI_REAL8,MPI_SUM, &
                        0,MPI_COMM_WORLD,ierror)
         call MPI_REDUCE(zpart,zparttot,1,MPI_REAL8,MPI_SUM, &
                        0,MPI_COMM_WORLD,ierror)
         call MPI_REDUCE(dens,densimp,1,MPI_REAL8,MPI_SUM, &
                        0,MPI_COMM_WORLD,ierror)

         if (myid.eq.0) then
            doubletot=doubletot/zparttot
            densimp=densimp/zparttot
         endif
         
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   Computation of the Green-Function on the imaginary axis           c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         if (myid.eq.0) then
            write(6,*)'Green-Function started'
            call flush(6)
         endif
         
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_ALLGATHERV(rlehm,lehmnumber,MPI_REAL8, &
                            rlehmtotal,recvnumberl,offsetnuml, &
                            MPI_REAL8,MPI_COMM_WORLD,ierror)   
         call MPI_ALLGATHERV(eig,eigennumber,MPI_REAL8, &
                            eigtotal,recvnumbere,offsetnume, &
                            MPI_REAL8,MPI_COMM_WORLD,ierror)   
         call computegimag(omoffset,omnumber,nsp1,nmax,Iwmax,prozessoren, &       
                            nleng,beta,eigtotal,rlehmtotal,om,Gwpart) 
         call MPI_BARRIER(MPI_COMM_WORLD,ierror) 
         call MPI_GATHERV(Gwpart(omoffset+1:omoffset+omnumber), &
                         omnumber,MPI_DOUBLE_COMPLEX,Gw, &
                         recvnumber,offsetnum, &
                         MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD, &
                         ierror)
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         
         if (myid.eq.0) then
            do i=0,Iwmax
               Gw(i)=Gw(i)/zparttot
            enddo
            write(6,*)'Green-Function finished'
            call flush(6)
         endif
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   End Green-Function Calculation                                    c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         if (myid.eq.0) then
            do i=0,Iwmax-1
               Gww(i)=Gw(i)
               Gww(i+Iwmax)=1.d0/ni(2*(Iwmax+i)+1)
            enddo
           
            do i=0,2*Iwmax-1
               omm=dimag(ni(2*i+1))
               call calcg0(omm,cdummy1,tpar,epsk,ns,Iwmax)
               cdummy1=1.d0/cdummy1
               g0wwand(i)=cdummy1
            end do
         endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   Computation of the Green-Function on the real axis                c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         if (myid.eq.0) then
            write(6,*)'Real Green-Function started'
            call flush(6)
         endif

         omnumber=(Iwmaxreal+1)/nprocs
         omrest=mod(Iwmaxreal+1,nprocs)
         do i=0,nprocs-1
            if (i.lt.omrest) then
               recvnumber(i)=omnumber+1
               offsetnum(i)=i*(omnumber+1)
            else
               recvnumber(i)=omnumber
               offsetnum(i)=omrest*(omnumber+1)+(i-omrest)*omnumber
            endif
         enddo
         if (myid.lt.omrest) then
            omnumber=omnumber+1
         endif
         if (myid.lt.omrest) then
            omoffset=myid*omnumber-1
         else
            omoffset=omrest*(omnumber+1)+(myid-omrest)*omnumber-1
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call computegreal(omoffset,omnumber,nsp1,nmax,Iwmaxreal,prozessoren, &        
                            nleng,beta,eigtotal,rlehmtotal,omr,Gwrealpart)

         call MPI_BARRIER(MPI_COMM_WORLD,ierror) 
         call MPI_GATHERV(Gwrealpart(omoffset+1:omoffset+omnumber), &
                         omnumber,MPI_DOUBLE_COMPLEX,Gwreal, &
                         recvnumber,offsetnum, &
                         MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD, &
                         ierror)
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         
         if (myid.eq.0) then
            do i=0,Iwmaxreal
               Gwreal(i)=Gwreal(i)/zparttot
            enddo
            write(6,*)'Real Green-Function finished'
            call flush(6)
         endif
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   End Real-Green-Function Calculation                               c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         if (myid.eq.0) then
            write (6,*) 'Now computing Matsubara self-energy'       
            write(6,*) 'Pi is equal to ', Pi
         
            do i=0, Iwmax
               self(i) =g0wwand(i)-1.d0/Gw(i)
               write(34,'(3f17.10)') om(i),dreal(self(i)),dimag(self(i))
            enddo
         
            write (6,*) 'Now computing self-energy on the real axis'       
            
            do i=0,Iwmaxreal
               sigre(i)=(0.d0,0.d0)
               do j=1,ns-1
                  sigre(i)=sigre(i)-tpar(j)**2/(omr(i)-Epsk(1+j))
               end do
               sigre(i)=sigre(i)+omr(i)-Epsk(1)- &
                       dconjg(Gwreal(i))/(Gwreal(i)*dconjg(Gwreal(i)))
            enddo
         endif
      endif

      if (myid.eq.0) then
         write(6,*)'Average Double Occupancy :',doubletot
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_BCAST(sigre,Iwmaxreal+1,MPI_DOUBLE_COMPLEX,0, &
                    MPI_COMM_WORLD,ierror)
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)

      if(lattice_type.eq.3) then
         open(57,file='dos.dat',form='formatted',status='unknown')
         call dos2(omnumber,omoffset,Iwmaxreal,ksteps,beta,xmu, &
                    Energya,omr,sigre,dospart)
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_GATHERV(dospart(omoffset+1:omoffset+omnumber), &
                         omnumber,MPI_REAL8,dos, &
                         recvnumber,offsetnum, &
                         MPI_REAL8,0,MPI_COMM_WORLD, &
                         ierror)
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         
         if (myid.eq.0) then
            sum=0.0d0
            do i=0,Iwmaxreal
               sum=sum+2.0d0*dos(i)*dreal(omr(1)-omr(0))/ &
                 (dexp(beta*dreal(omr(i)))+1.0d0)
               write(57,'(3f17.10)') dreal(omr(i)), dos(i),sum
            enddo
            write(6,*) 'Final Spectral weight check=',sum
         endif
      endif
        
      if(lattice_type.eq.4) then
         open(57,file='dos.dat',form='formatted',status='unknown')
         do i=0,Iwmaxreal
            dos(i)=0.d0
            do j=0,4999
               E=-1.d0+2.d0*dfloat(j)/dfloat(5000) +1.d0/dfloat(5000)
                  dos(i)=dos(i)-2.d0/(pi)**2*dsqrt(1.d0-E**2)* &
                      dimag(1.d0/(omr(i)-Epsk(1)-E-sigre(i))) &
                      /dfloat(5000)
            enddo
               if (myid.eq.0) then
                  write(57,'(2f17.10)') dreal(omr(i)),dos(i)
               endif
         enddo
      endif

      if (myid.eq.0) then
         write(*,*)'------------------------------------------------'
         write(*,*) 'final Anderson parameters '
         write(*,*)'------------------------------------------------'
         write(*,*)'cut from the following line...'
         write(*,*) 'Eps(k) '
         do i=2,ns
           write(*,'(2f27.16)')epsk(i)
         enddo
         write(*,*)'V(k) '
         do i=1,ns-1
            write(*,'(2f27.16)')tpar(i)
         enddo
         write(*,'(f27.16,"   #chemical potential")')xmu      
         write(*,*)'...to the line above and copy into hubb.andpar'
         write(*,*)'------------------------------------------------'
 
         do i=0,2*Iwmax-1
            write(90,'(3f27.16)') dimag(ni(2*i+1)) &
                ,dreal(Gww(i)),dimag(Gww(i))
         enddo
         
         do i=0,Iwmax
            G0w(i)=1.d0/G0w(i)
            write(92,'(3f27.16)') dimag(ni(2*i+1)) &
                ,dreal(G0w(i)),dimag(G0w(i))
         end do
               
         do i=0,2*Iwmax-1
            write(93,'(3f27.16)')dimag(ni(2*i+1)),dreal(G0wwand(i)), &
                dimag(G0wwand(i))
         end do
         
         do i=0,Iwmaxreal
            write(91,'(5f27.16)')dreal(omr(i)), &
                dreal(Gwreal(i)),dimag(Gwreal(i)),dreal(sigre(i)), &
                dimag(sigre(i))
         end do
         
         sig0=g0w(0)-1.d0/gw(0)
         write(6,*) 'sig0 = ',sig0
         zed=1.d0-dimag(sig0)/piob
         zed=1.d0/zed
         write(6,*) 'zed = ', zed
         call flush(6)
          
         dostest=0.d0
         do i=1,ns-1
            dostest=dostest+tpar(i)**2
         enddo
         write(*,*) 'DOSTEST=', dostest
         write(*,*) 'U | xmu | densimp | densold | emin | zpart'
         write(*,*) uhub,xmu,densimp,densold,emin,zpart
         write(36,*) zparttot
         write(37,*) doubletot
         write(38,*) densimp
!========+=========+=========+=========+=========+=========+=========+=$
!     END OF Iteration 
!========+=========+=========+=========+=========+=========+=========+=$
         write(6,'(a20)')'     End lisalanc '
         write(6,'(a60)')'========================================'
      endif
      
      call MPI_FINALIZE(ierror)

end







      !TODO: move these routines somewhere else


subroutine calcg0(omega,g0and,tpar,epsk,ns,Iwmax)
    use types
    implicit none
    integer,     intent(in)  :: ns,Iwmax
    real(dp),    intent(in)  :: omega,tpar(ns),epsk(ns)
    complex(dp), intent(out) :: g0and
    integer i

    g0and=(0.d0,1.d0)*omega-Epsk(1)
    do i=1,ns-1
        g0and=g0and-tpar(i)**2/((0.d0,1.d0)*omega-Epsk(1+i))
    end do
    g0and=dconjg(g0and)/(dconjg(g0and)*g0and)
end

!========+=========+=========+=========+=========+=========+=========+=$
      subroutine search(difference,fmin,Nitermax,hess,g,xtemp,w,xprmt,nmpara, &
          tpar,epsk,ns,piob,Iwmax,G0w,G0wand,om,symm)
      use types
      implicit none
      integer nmpara,ns,Iwmax,nbparm,icount,i,mode,Nitermax,iexit,nsym
      real(dp) dfn,deps,fmin,hh
      real(dp) hess(nmpara**2),g(nmpara),xtemp(nmpara), &
          w(nmpara**2),xprmt(nmpara)
      real(dp) tpar(ns),epsk(ns),piob  
      real(dp) om(0:Iwmax)
      logical symm
      complex*16 G0w(0:Iwmax),G0wand(0:Iwmax)
      data iexit/0/
      data hh/1.e-5/
      external difference
 
!     number of parameters to be optimized (per spin):
!     eps(2)...... eps(ns) --->  ns-1
!     tpar(1)........ tpar(ns-1) --->  ns-1

      if(symm) then
         nsym=ns/2
         nbparm=2*(nsym-1)
         if (nbparm.gt.nmpara) stop ' nmpara too small'
         icount=0
!         nsym=(ns-1)/2
         nbparm=2*(nsym)
         do i=2,nsym+1
            icount=icount+1	
            xtemp(icount)=Epsk(i)
 11         continue
         end do
         
         do i=1,nsym
            icount=icount+1
            xtemp(icount)=tpar(i)
 12         continue
         end do
         
         do i=nbparm+1,nmpara
            xtemp(i)=0.d0
         end do
      
         do i=1,nbparm
            xprmt(i)=dabs(xtemp(i))+1.d-15
         end do
      else  
         nbparm=2*(ns-1)
         if (nbparm.gt.nmpara) stop ' nmpara too small'
         icount=0
         do i=2,ns
            icount=icount+1	
            xtemp(icount)=Epsk(i)
 111        continue
         end do
         
         do i=1,ns-1
            icount=icount+1
            xtemp(icount)=tpar(i)
 112         continue
         end do
         
         do i=nbparm+1,nmpara
            xtemp(i)=0.d0
         end do
         
         do i=1,nbparm
            xprmt(i)=dabs(xtemp(i))+1.d-15
         end do
      endif

      

      mode=1
      dfn=-.5d0
      deps=.00001d0

      call minimize(difference,nbparm,xtemp,fmin,g,hess,w, &
          dfn,xprmt,hh,deps,mode,Nitermax,iexit, &
          tpar,epsk,piob,Iwmax,ns,G0wand,G0w,om,symm)
      write (6,30) iexit,fmin
 30   format(' iexit fmin ',i5,e14.6)
      end
!========+=========+=========+=========+=========+=========+=========+=$
subroutine minimize(funct, n, x, f, g, h, w, dfn, xm, &
       hh, eps, mode, maxfn, iexit,tpar,epsk,piob, &
          Iwmax,ns,G0wand,G0w,om,symm)
      use types
      implicit none
      integer ns,Iwmax
      integer np,n,n1,nn,is,iu,iv,ib,idiff,iexit,mode,ij,maxfn,i,j, &
          i1,jk,ik,k,itn,ifn,link,int
      real(dp) z,zz,dmin,f,df,dfn,aeps,eps,alpha,ff,tot,f1,f2,half, &
          gys,dgs,sig,hh,gs0
      real(dp) tpar(ns),epsk(ns),piob
      complex*16 G0w(0:Iwmax),G0wand(0:Iwmax)
      real(dp)  x(*), g(*), h(*), w(*), xm(*)
      real(dp) om(0:Iwmax)
      logical symm
      external funct
      data half /0.5d0/


      np = n + 1
      n1 = n - 1
      nn=(n*np)/2
      is = n
      iu = n
      iv = n + n
      ib = iv + n
      idiff = 1
      iexit = 0
      if (mode .eq. 3) go to 15
      if (mode .eq. 2) go to 10
      ij = nn + 1
      do 5 i = 1, n
      do 6 j = 1, i
      ij = ij - 1
   6  h(ij) = 0.d0
   5  h(ij) = 1.d0
      go to 15
  10  continue
      ij = 1
      do 11 i = 2, n
      z = h(ij)
      if (z .le. 0.d0) return
      ij = ij + 1
      i1 = ij
      do 11 j = i, n
      zz = h(ij)
      h(ij) = h(ij) / z
      jk = ij
      ik = i1
      do 12 k = i, j
      jk = jk + np - k
      h(jk) = h(jk) - h(ik) * zz
      ik = ik + 1
  12  continue
      ij = ij + 1
  11  continue
      if (h(ij) .le. 0.d0) return
  15  continue
      ij = np
      dmin = h(1)
      do 16 i = 2, n
      if (h(ij) .ge. dmin) go to 16
      dmin = h(ij)
  16  ij = ij + np - i
      if (dmin .le. 0.d0) return
      z = f
      itn = 0
      call funct (n, x, f,epsk,tpar,ns,piob,Iwmax,G0w,G0wand,om,symm)
      ifn = 1
      df = dfn
      if (dfn .eq. 0.d0) df = f - z
      if (dfn .lt. 0.d0) df = abs (df * f)
      if (df .le. 0.d0) df = 1.d0
  17  continue
      do 19 i = 1, n
      w(i) = x(i)
  19  continue
      link = 1
      if (idiff - 1) 100, 100, 110
  18  continue
      if (ifn .ge. maxfn) go to 90
  20  continue
  21  continue
      itn = itn + 1
      w(1) = -g(1)
      do 22 i = 2, n
      ij = i
      i1 = i - 1
      z = -g(i)
      do 23 j = 1, i1
      z = z - h(ij) * w(j)
      ij = ij + n - j
  23  continue
  22  w(i) = z
      w(is+n) = w(n) / h(nn)
      ij = nn
      do 25 i = 1, n1
      ij = ij - 1
      z = 0.d0
      do 26 j = 1, i
      z = z + h(ij) * w(is+np-j)
      ij = ij - 1
  26  continue
  25  w(is+n-i) = w(n-i) / h(ij) - z
      z = 0.d0
      gs0 = 0.d0
      do 29 i = 1, n
      if (z * xm(i) .ge. abs (w(is+i))) go to 28
      z = abs (w(is+i)) / xm(i)
  28  gs0 = gs0 + g(i) * w(is+i)
  29  continue
      aeps = eps / z
      iexit = 2
      if (gs0 .ge. 0.d0) go to 92
      alpha = -2.d0 * df / gs0
      if (alpha .gt. 1.d0) alpha = 1.d0
      ff = f
      tot = 0.d0
      int = 0
      iexit = 1
  30  continue
      if (ifn .ge. maxfn) go to 90
      do 31 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  31  continue
      call funct (n, w, f1,epsk,tpar,ns,piob,Iwmax,G0w,G0wand,om)
      ifn = ifn + 1
      if (f1 .ge. f) go to 40
      f2 = f
      tot = tot + alpha
  32  continue
      do 33 i = 1, n
      x(i) = w(i)
  33  continue
      f = f1
      if (int - 1) 35, 49, 50
  35  continue
      if (ifn .ge. maxfn) go to 90
      do 34 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  34  continue
      call funct (n, w, f1,epsk,tpar,ns,piob,Iwmax,G0w,G0wand,om)
      ifn = ifn + 1
      if (f1 .ge. f) go to 50
      if ((f1 + f2 .ge. f + f) .and. (7.0d0 * f1 + 5.0d0 * f2 .gt. 12.0d0 * f)) int = 2
      tot = tot + alpha
      alpha = 2.d0 * alpha
      go to 32
  40  continue
      if (alpha .lt. aeps) go to 92
      if (ifn .ge. maxfn) go to 90
      alpha = half * alpha
      do 41 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  41  continue
      call funct (n, w, f2,epsk,tpar,ns,piob,Iwmax,G0w,G0wand,om)
      ifn = ifn + 1
      if (f2 .ge. f) go to 45
      tot = tot + alpha
      f = f2
      do 42 i = 1, n
      x(i) = w(i)
  42  continue
      go to 49
  45  continue
      z = 0.1d0
      if (f1 + f .gt. f2 + f2) z = 1.d0 + half * (f - f1) / (f + f1 - f2 - f2)
      if (z .lt. 0.1d0) z = 0.1d0
      alpha = z * alpha
      int = 1
      go to 30
  49  continue
      if (tot .lt. aeps) go to 92
  50  continue
      alpha = tot
      do 56 i = 1, n
      w(i) = x(i)
      w(ib+i) = g(i)
  56  continue
      link = 2
      if (idiff - 1) 100, 100, 110
  54  continue
      if (ifn .ge. maxfn) go to 90
      gys = 0.d0
      do 55 i = 1, n
      w(i) = w(ib+i)
      gys = gys + g(i) * w(is+i)
  55  continue
      df = ff - f
      dgs = gys - gs0
      if (dgs .le. 0.d0) go to 20
      link = 1
      if (dgs + alpha * gs0 .gt. 0.d0) go to 52
      do 51 i = 1, n
      w(iu + i) = g(i) - w(i)
  51  continue
      sig = 1.d0 / (alpha * dgs)
      go to 70
  52  continue
      zz = alpha / (dgs - alpha * gs0)
      z = dgs * zz - 1.d0
      do 53 i = 1, n
      w(iu+i) = z * w(i) + g(i)
  53  continue
      sig = 1.d0 / (zz * dgs * dgs)
      go to 70
  60  continue
      link = 2
      do 61 i = 1, n
      w(iu+i) = w(i)
  61  continue
      if (dgs + alpha * gs0 .gt. 0.d0) go to 62
      sig = 1.d0 / gs0
      go to 70
  62  continue
      sig = -zz
  70  continue
      w(iv+1) = w(iu+1)
      do 71 i = 2, n
      ij = i
      i1 = i - 1
      z = w(iu+i)
      do 72 j = 1, i1
      z = z - h(ij) * w(iv+j)
      ij = ij + n - j
  72  continue
      w(iv+i) = z
  71  continue
      ij = 1
      do 75 i = 1, n
      z = h(ij) + sig * w(iv+i) * w(iv+i)
      if (z .le. 0.d0) z = dmin
      if (z .lt. dmin) dmin = z
      h(ij) = z
      w(ib+i) = w(iv+i) * sig / z
      sig = sig - w(ib+i) * w(ib+i) * z
      ij = ij + np - i
  75  continue
      ij = 1
      do 80 i = 1, n1
      ij = ij + 1
      i1 = i + 1
      do 80 j = i1, n
      w(iu+j) = w(iu+j) - h(ij) * w(iv+i)
      h(ij) = h(ij) + w(ib+i) * w(iu+j)
      ij = ij + 1
  80  continue
      go to (60, 20), link
  90  continue
      iexit = 3
      go to 94
  92  continue
      if (idiff .eq. 2) go to 94
      idiff = 2
      go to 17
  94  continue
      return
 100  continue
      do 101 i = 1, n
         z = hh * xm(i)
         w(i) = w(i) + z
         call funct (n, w, f1,epsk,tpar,ns,piob,Iwmax,G0w,G0wand,om)
      g(i) = (f1 - f) / z
      w(i) = w(i) - z
 101  continue
      ifn = ifn + n
      go to (18, 54), link
 110  continue
      do 111 i = 1, n
      z = hh * xm(i)
      w(i) = w(i) + z
      call funct (n, w, f1,epsk,tpar,ns,piob,Iwmax,G0w,G0wand,om)
      w(i) = w(i) - z - z
      call funct (n, w, f2, epsk,tpar,ns,piob,Iwmax,G0w,G0wand,om)
      g(i) = (f1 - f2) / (2.d0 * z)
      w(i) = w(i) + z
 111  continue
      ifn = ifn + n + n
      go to (18, 54), link
      end 
!========+=========+=========+=========+=========+=========+=========+=$
      subroutine diag(ns,nsp1,nblock,ntot,nmax,myblock,blocksize, &
                    b,ntable,ioffset,imat,xmu,uhub,hmag,tpar,epsk, &
                    emin,eig,zeig)
      use AuxRoutines
      implicit none
!Input-Variables
      integer ns,nsp1,nblock,ntot,nmax,myblock,blocksize
      integer b(ns),ntable(ntot),ioffset(nblock+1)
      integer imat(nmax,nmax)
      real*8 xmu,uhub,hmag
      real*8 tpar(ns),epsk(ns)
!Output-Variables
      real*8 emin
      real*8 eig(nmax),zeig(nmax,nmax)
!Subroutine-Internal-Variables
      integer j,k,is,neig,nup,ndo,icount,info
      integer abas(ns)
      real*8 work(3*nmax)
      real*8 realmat(nmax*(nmax+1)/2),xmat(nmax,nmax)
      
      epsk(1)=-xmu
      emin=1.d16

      do j=1,nmax
         do k=j+1,nmax
            xmat(j,k)=0.d0
            xmat(k,j)=0.d0
         enddo
         xmat(j,j)=100000.d0
      enddo

      do j=1,blocksize
         call findnupndo(ntable(ioffset(myblock)+j)-1,nup,ndo, &
             b,abas,ns)
         xmat(j,j)=0.d0
         do is=1,ns
            xmat(j,j)=xmat(j,j)+epsk(is)*dabs(dfloat(abas(is)))
         enddo
         if (abas(1).eq.2) then
            xmat(j,j)=xmat(j,j)+uhub
         else
            xmat(j,j)=xmat(j,j)-hmag*dfloat(abas(1))
         endif
         do k=j+1,blocksize
            if (imat(j,k).ne.0) then
               xmat(j,k)=xmat(j,k)+tpar(abs(imat(j,k))-1)* &
                                  dfloat(imat(j,k))/ &
                                  abs(dfloat(imat(j,k)))
               xmat(k,j)=xmat(k,j)+tpar(abs(imat(j,k))-1)* &
                                  dfloat(imat(j,k))/ &
                                  abs(dfloat(imat(j,k)))
            endif
         enddo
      enddo
         
      icount=0
         
      do j=1,nmax
         do k=j,nmax
            icount=icount+1
            realmat(icount)=xmat(k,j)
         enddo
      enddo
            
       
      call dspev('V','L',nmax,realmat,eig,zeig, &
                nmax,work,info)

      do neig=1,blocksize
         if (eig(neig).lt.emin) emin=eig(neig)
      enddo

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lehmann(ns,nsp1,nblock,ntot,nmax,myblock, &
                        blocksize,blocksize1,ioffset,idouble,idmat, &
                        xmu,uhub,hmag,beta,eig,zeig,eigold,zold, &
                        double,zpart,dens,rlehm)

      implicit none
!Input-Variables
      integer ns,nsp1,nblock,ntot,nmax,myblock,blocksize
      integer blocksize1
      integer ioffset(nblock+1),idmat(ntot)
      integer idouble(ntot)
      real*8 xmu,uhub,hmag,beta
      real*8 eig(nmax),zeig(nmax,nmax)
      real*8 eigold(nmax),zold(nmax,nmax)
!Output-Variables
      real*8 double,zpart,dens
      real*8 rlehm(nmax,nmax)
!Subroutine-Internal-Variables
      integer i,j,k,meig,neig
      real*8 rmel

      double=0.0d0
      zpart=0.0d0
      dens=0.0d0

      do i=1,nmax
         do j=1,nmax
            rlehm(i,j)=0.d0
         enddo
      enddo

      if (mod(myblock-1,nsp1).ne.0) then
!     loop over |ndo,nup-1> basis states
         do neig=1,blocksize
            do meig=1,blocksize1
               rmel=0.d0
               do j=1,blocksize1
                  if (idmat(ioffset(myblock-1)+j).ne.0) then
                     rmel=rmel+zold(j,meig)* &
                              zeig(idmat(ioffset(myblock-1)+j),neig)
                  endif
               enddo
               rlehm(neig,meig)=rmel**2
               dens=dens+rlehm(neig,meig)*dexp(-beta*eig(neig))
            enddo
         enddo
      endif
      
      dens=2.0d0*dens

      do neig=1,blocksize
         do j=1,blocksize
            double=double+ &
                zeig(j,neig)*zeig(j,neig)* &
                dfloat(idouble(ioffset(myblock)+j))* &
                dexp(-beta*eig(neig))
         enddo
         zpart=zpart+dexp(-beta*eig(neig))
      enddo

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine computegimag(omoffset,omnumber,           &
                                nsp1,nmax,Iwmax,prozessoren,nleng, &
                                beta,eig,rlehm,om,Gw) 

      implicit none

!Input-Variables
      integer omoffset,omnumber,nsp1,nmax,Iwmax,prozessoren
      integer nleng(nsp1**2)
      real*8 beta
      real*8 eig(nmax,prozessoren)
      real*8 rlehm(nmax,nmax,prozessoren),om(0:Iwmax)
!Output-Variables
      complex*16 Gw(0:Iwmax)
!Subroutine-Internal-Variables
      integer i,iomega,meig,neig
      
      do i=0,Iwmax
         Gw(i)=dcmplx(0.d0,0.d0)
      enddo

      do i=1,nsp1**2
      if (mod(i-1,nsp1).ne.0) then
         do iomega=omoffset+1,omoffset+omnumber
            do neig=1,nleng(i)
               do meig=1,nleng(i-1)
                  Gw(iomega)=Gw(iomega)+rlehm(neig,meig,i) &
                      /((0.d0,1.d0)*om(iomega)- &
                      (eig(neig,i)-eig(meig,i-1)))* &
                      (dexp(-beta*eig(neig,i))+ &
                       dexp(-beta*eig(meig,i-1))) 
               enddo                  
            enddo
         enddo
      endif
      enddo

      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine computegreal(omoffset,omnumber, &
                                nsp1,nmax,Iwmaxreal,prozessoren, &        
                                nleng,beta,eig,rlehm,omr,Gw)

      implicit none

!Input-Variables
      integer omoffset,omnumber,nsp1,nmax,Iwmaxreal,prozessoren
      integer nleng(nsp1**2)
      real*8 beta
      real*8 eig(nmax,prozessoren)
      real*8 rlehm(nmax,nmax,prozessoren)
      complex*16 omr(0:Iwmaxreal)
!Output-Variables
      complex*16 Gw(0:Iwmaxreal)
!Subroutine-Internal-Variables
      integer i,iomega,meig,neig
      
      do i=0,Iwmaxreal
         Gw(i)=dcmplx(0.d0,0.d0)
      enddo
      
      do i=1,nsp1**2
      if (mod(i-1,nsp1).ne.0) then
         do iomega=omoffset+1,omoffset+omnumber
            do neig=1,nleng(i)
               do meig=1,nleng(i-1)
                  Gw(iomega)=Gw(iomega)+rlehm(neig,meig,i) &
                      /(omr(iomega)- &
                      (eig(neig,i)-eig(meig,i-1)))* &
                      (dexp(-beta*eig(neig,i))+ &
                       dexp(-beta*eig(meig,i-1)))
               enddo                  
            enddo
         enddo
      endif
      enddo

      return
      end
