module LatticeRoutines
    use types
    implicit none
    contains

pure real(dp) function dispersion_sc3D(t,t1,t2,x,y,z)
    implicit none
    real(dp),   intent(in)  :: t,t1,t2,x,y,z
    dispersion_sc3D =-2.0d0*t*(dcos(x)+dcos(y)+dcos(z))
end

pure real(dp) function dispersion_sc2D(t,t1,t2,x,y)
    implicit none
    real(dp),   intent(in)  :: t,t1,t2,x,y
    dispersion_sc2D =-2.0d0*t*(dcos(x)+dcos(y)) + 4.0*t1*(dcos(x)*dcos(y)) - 2*t2*(dcos(2*x)+dcos(2*y))
end

pure real(dp) function dispersion_fcc(t,t1,t2,x,y,z)
    implicit none
    real(dp),   intent(in)  :: t,t1,t2,x,y,z
    dispersion_fcc= -4.d0*t*(cos((x+y-z)/2)*cos((x-y+z)/2)+cos((x+y-z)/2)*cos((-x+y+z)/2)+cos((x-y+z)/2)*cos((-x+y+z)/2))
end

! TODO: documentation here
subroutine dos2(omnumber,omoffset,Iwmaxreal,ksteps, &
                beta,xmu,Energya,omr,sigre,dos)
    implicit none
    integer,    intent(in)  :: omnumber,omoffset,Iwmaxreal,ksteps
    real(dp),   intent(in)  :: beta,xmu
    real(dp),   intent(in)  :: Energya(0:ksteps,0:ksteps)
    complex(dp),intent(in)  :: omr(0:Iwmaxreal),sigre(0:Iwmaxreal)
    real(dp),   intent(out) :: dos(0:Iwmaxreal)
    real(dp), parameter :: Pi = 4 * atan(1.d0)
    integer i,ikx,iky
    real(dp) weightx,weighty,weightm
 
    if (omoffset.eq.-1)then
        write(6,*) 'Computing the interacting DOS'
        call flush(6)
    endif

    do i=omoffset+1,omoffset+omnumber
        dos(i)=0.d0
        do ikx=0,ksteps
            weightx=1.0d0
            weightm=1.0d0
            if ((ikx.eq.0).or.(ikx.eq.ksteps)) then
                weightx=0.5d0
            endif
            do iky=0,ikx
                weighty=1.0d0
                if ((iky.eq.0).or.(iky.eq.ikx)) then
                    weighty=0.5d0
                endif
                if (((iky.eq.0).and.(ikx.eq.0)).or. &
                    ((iky.eq.ksteps).and.(ikx.eq.ksteps))) then
                    weightm=0.5d0
                endif
                dos(i)=dos(i)-(2.0d0*weightx*weighty*weightm/Pi)* &
                         dimag(1.d0/(omr(i)-Energya(ikx,iky)+xmu &
                         -sigre(i)))
            enddo
        enddo
        dos(i)=dos(i)/(dfloat(ksteps)*dfloat(ksteps))
    enddo
    return
end

subroutine selfconst2_full(omnumber,omoffset,Iwmax,ksteps, &
                      beta,xmu,om,Energya,g0wand,Gw,Gloc)
      
      implicit none
!Input-Variables
    integer,     intent(in)  :: omnumber,omoffset,Iwmax,ksteps
    real(dp),    intent(in)  :: beta,xmu
    real(dp),    intent(in)  :: om(0:Iwmax),Energya(0:ksteps,0:ksteps)
    complex(dp), intent(inout) :: g0wand(0:Iwmax),Gw(0:Iwmax)
    complex(dp), intent(out) :: Gloc(0:Iwmax)
    integer i,ikx,iky
    real(dp) weightx,weighty,weightm
    complex(dp) W(omoffset+1:omoffset+omnumber)
    complex(dp) gand(omoffset+1:omoffset+omnumber)

    do i=omoffset+1,omoffset+omnumber
        gand(i)=1.d0/g0wand(i)
        Gw(i)=1.d0/Gw(i)
        W(i)= (0.d0,1.d0)*om(i)-gand(i)+Gw(i)+xmu!-tpri
    enddo 
      
    do i=omoffset+1,omoffset+omnumber
        Gloc(i)=dcmplx(0.d0,0.d0)
        do ikx=0,ksteps
            do iky=0,ksteps
                Gloc(i)=Gloc(i)+1.0d0/(W(i)-Energya(ikx,iky))
            enddo
        enddo
        Gloc(i)=Gloc(i)/(dfloat(ksteps)*dfloat(ksteps))
        Gloc(i)=1.d0/Gloc(i)+gand(i)-Gw(i)
        Gloc(i)=0.5d0/Gloc(i)+0.5d0*g0wand(i)
    enddo
    return 
end

subroutine selfconst2(omnumber,omoffset,Iwmax,ksteps, &
                      beta,xmu,om,Energya,g0wand,Gw,Gloc)
      
      implicit none
!Input-Variables
    integer,     intent(in)  :: omnumber,omoffset,Iwmax,ksteps
    real(dp),    intent(in)  :: beta,xmu
    real(dp),    intent(in)  :: om(0:Iwmax),Energya(0:ksteps,0:ksteps)
    complex(dp), intent(inout) :: g0wand(0:Iwmax),Gw(0:Iwmax)
    complex(dp), intent(out) :: Gloc(0:Iwmax)
    integer i,ikx,iky
    real(dp) weightx,weighty,weightm
    complex(dp) W(omoffset+1:omoffset+omnumber)
    complex(dp) gand(omoffset+1:omoffset+omnumber)

    do i=omoffset+1,omoffset+omnumber
        gand(i)=1.d0/g0wand(i)
        Gw(i)=1.d0/Gw(i)
        W(i)= (0.d0,1.d0)*om(i)-gand(i)+Gw(i)+xmu!-tpri
    enddo 
      
    do i=omoffset+1,omoffset+omnumber
        Gloc(i)=dcmplx(0.d0,0.d0)
        do ikx=0,ksteps
            weightx=1.0d0
            weightm=1.0d0
            if ((ikx.eq.0).or.(ikx.eq.ksteps)) then
                weightx=0.5d0
            endif
            do iky=0,ikx
                weighty=1.0d0
                if ((iky.eq.0).or.(iky.eq.ikx)) then
                    weighty=0.5d0
                endif
                if (((iky.eq.0).and.(ikx.eq.0)).or. &
                    ((iky.eq.ksteps).and.(ikx.eq.ksteps))) then
                    weightm=0.5d0
                endif
                Gloc(i)=Gloc(i)+2.0d0*weightx*weighty*weightm/(W(i)- &
                        Energya(ikx,iky))
            enddo
        enddo
        Gloc(i)=Gloc(i)/(dfloat(ksteps)*dfloat(ksteps))
        Gloc(i)=1.d0/Gloc(i)+gand(i)-Gw(i)
        Gloc(i)=0.5d0/Gloc(i)+0.5d0*g0wand(i)
    enddo
    return 
end


subroutine selfconst3(omnumber,omoffset,Iwmax,ksteps, &
                      beta,xmu,om,Energya,g0wand,Gw,Gloc)
      
    integer,     intent(in)  :: omnumber,omoffset,Iwmax,ksteps
    real(dp),    intent(in)  :: beta,xmu
    real(dp),    intent(in)  :: om(0:Iwmax),Energya(0:ksteps,0:ksteps,0:ksteps)
    complex(dp), intent(inout) :: g0wand(0:Iwmax),Gw(0:Iwmax)
    complex(dp), intent(out) :: Gloc(0:Iwmax)
    integer i,ikx,iky,ikz
    real(dp) weightx,weighty,weightz,weightm
    complex(dp) W(omoffset+1:omoffset+omnumber)
    complex(dp) gand(omoffset+1:omoffset+omnumber)

    do i=omoffset+1,omoffset+omnumber
        gand(i)=1.d0/g0wand(i)
        Gw(i)=1.d0/Gw(i)
        W(i)= (0.d0,1.d0)*om(i)-gand(i)+Gw(i)+xmu!-tpri
    enddo 
      
    do i=omoffset+1,omoffset+omnumber
        Gloc(i)=dcmplx(0.d0,0.d0)
        do ikx=0,ksteps
            weightx=1.0d0
            if ((ikx.eq.0).or.(ikx.eq.ksteps)) then
                weightx=weightx*0.5d0
            endif
            do iky=0,ikx
                weighty=weightx
                if ((iky.eq.0).or.(iky.eq.ksteps)) then
                    weighty=weighty*0.5d0
                endif
                do ikz=0,iky
                    weightz=weighty
                    if ((ikz.eq.0).or.(ikz.eq.ksteps)) then
                        weightz=weightz*0.5d0
                    endif
                    weightm=weightz*dfloat(6/ &
                       ((1+iky)/(1+ikx)+(1+ikz)/(1+iky)+ &
                       3*((1+ikz)/(1+ikx))+1))
                    Gloc(i)=Gloc(i)+weightm/(W(i)- &
                       Energya(ikx,iky,ikz))
                enddo
            enddo
        enddo
        Gloc(i)=Gloc(i)/(dfloat(ksteps)**3)
        Gloc(i)=1.d0/Gloc(i)+gand(i)-Gw(i)
        Gloc(i)=1.0d0/Gloc(i)   !without damping
    enddo
    return
 end


 subroutine selfconst2_magnetic(omnumber,omoffset,Iwmax,ksteps,L,   &  !Integer-Input     
                                   beta,xmu,B,                    &      !Real-Input  
                                   om,Energya,                    &      !Real-Array-Input
                                   Xi,                            &      !Complex-Input      
                                   g0wand,Gw,                     &      !Complex-Array-Input
                                   Gloc)                                !Complex-Array-Output
     
     implicit none

!Input-Variables
     integer omnumber,omoffset                                             !omnumber is the number of all omegas, omoffset normally -omnumber/2
     integer Iwmax                                                         !maximum value of omega frequency
     integer ksteps                                                        !twice many steps in the kGrid, e.g. ksteps = 2 gives kx = -pi/2, 0, pi/2, pi
     integer kstepsY                                                       !since Brillouin in y direction smaller but same finesse, have less ksteps
     integer myblock
     integer L                                                             !denominator in B, in Georgs paper it is q
     real*8 beta                                                           !Beta as 1/T with temperature T
     real*8 xmu                                                            !chemical potential
     real*8 B                                                              !magnetic field is given from outside and can have p as well
!B...magnetic field in units of [Phi_0 (Fluxquantum)/area of unit cell(=a^2=1)]
!L...Size of unit cell in y-direction -> normally choosen such that
!for testing: set L finite and B=0!
     real*8 om(0:Iwmax),Energya(0:ksteps,0:ksteps)                         !om are the frequencies omega, Energya are tight binding energies
     complex*16 Xi                                                         !just the imaginary number i
     complex*16 g0wand(0:Iwmax),Gw(0:Iwmax)                                !g0wand is G_0 of anderson Gw is full Greens function
!Output-Variables
     complex*16 gloc(0:Iwmax,1:L,1:L)                                      !the output function which is the local Greens function G_{L, L^'}(omega) of each band L in dependence of omega 
!Subroutine-Internal-Variables
     integer i,j,k,ikx,iky                                                 !running variables in do loops
     real*8 weightx,weighty,weightm                                        !since kx and ky are symmetric to 0 we do not have to calculate all, but can multiply them with weights, eg 4 at the diagonale
     real*8 kx,ky                                                          !real values of kx and ky in a Brillouin zone and not the running value ikx and iky
     real*8 Pi                                                             !will store the value of pi
     real*8 t,t1,t2,checkone                                               !t is hopping parameter, checkone to check the momentum sum

     complex*16 W(omoffset+1:omoffset+omnumber)                            !basically i\nu + \mu - (G_0^{-1} - G^{-1}) = i\nu + \mu - \Sigma
     complex*16 gand(omoffset+1:omoffset+omnumber)                         !anderson Greens function
     complex*16 Gw1(omoffset+1:omoffset+omnumber)                          !G^{-1}
     complex*16 epsmat(1:L,1:L)                                            !the epsilon matrix I give
     complex*16 ginv(1:L,1:L)                                              !only used to store G^{-1}
     complex*16 gk                                                         !never needed
!Varaibles for inversion subroutines
     INTEGER infoinv                                                       !written by Lapack, if exited succesfull or not
     INTEGER ipiv(1:L)                                                     !contaions some pivot indices for Lapack
     COMPLEX*16 workinv(1:10*L)                                            !only needed in Lapack and has those dimensions, can return optimal WORK
     include 'tpri.dat'
     Pi=dacos(-1.0d0)                                                      !get value of Pi by the arccos(-1)

!      write(6,*)B

     do i=omoffset+1,omoffset+omnumber                                    !We basically go from -omnumber + 1 to omnumber and therefore we look at each frequency by using the running variable of omega
        gand(i)=1.d0/g0wand(i)                                             !G_0^{-1}
        Gw1(i)=1.d0/Gw(i)                                                  !G^{-1}
        W(i)= Xi*om(i)-gand(i)+Gw1(i)+xmu!-tpri                            !basically i\nu + \mu - (G_0^{-1} - G^{-1}) = i\nu + \mu - \Sigma
     enddo

!Initialize dispersion matrix
     do j=1,L
        do k=1,L
           epsmat(j,k)=dcmplx(0.0d0,0.0d0)                             !Set each entry of the epsilon matrix 0
           do i=0,Iwmax
              gloc(i,j,k)=dcmplx(0.0d0,0.0d0)                          !Set each entry in G_{L, L^'}(omega) 0
           enddo
        enddo
     enddo

     checkone=0.0d0

     kstepsY=ksteps/L
     
     do ikx=-ksteps+1,ksteps
        kx=Pi*dfloat(ikx)/dfloat(ksteps)                               !get kx from the running variable ikx
        do iky=-kstepsY+1,kstepsY                                      !it is totally fine to have less points in y direction since the finesse should stay the same and ky is shorter in generell
           ky=Pi*dfloat(iky)/dfloat(ksteps)                            !therefore do not need changes here, since it is the same besides less k points
           checkone=checkone+1.0d0                                     !by this we sum over all (kx,ky) tuples which should be 2*ksteps * 2*ksteps/L = 2*ksteps * 2*kstepsY using integer division 
  
!Construct dispersion matrix for given values of kx and ky
           if (L.gt.2) then                                            !if epsilon matrix is bigger than a 2x2 matrix, L is bigger than 2
              epsmat(1,1)=dcmplx(-2.0d0*t*dcos(kx),0.0d0)
              epsmat(1,2)=dcmplx(-t,0.0d0)
              epsmat(1,L)=dcmplx(-t*dcos(dfloat(L)*ky),0.0d0)+    &     !top right entry
                  dcmplx(0.0d0,-t*dsin(dfloat(L)*ky))
              epsmat(L,L)=dcmplx(-2.0d0*t*                         &    !bottom right entry
                  dcos(kx+dfloat(L-1)*B*2.0d0*Pi),0.0d0)
              epsmat(L,L-1)=dcmplx(-t,0.0d0)                           !bottom right but one to the left entry
              epsmat(L,1)=dcmplx(-t*dcos(dfloat(L)*ky),0.0d0)-    &     !bottom left entry
                  dcmplx(0.0d0,-t*dsin(dfloat(L)*ky))                 !used cos and sin instead of Euler because maybe faster
              
              do j=2,L-1
                 epsmat(j,j)=dcmplx(-2.0d0*t*                    &      !diagonale entries
                     dcos(kx+dfloat(j-1)*B*2.0d0*Pi),0.0d0)
                 epsmat(j,j-1)=dcmplx(-t,0.0d0)                        !lower diagonale entries
                 epsmat(j,j+1)=dcmplx(-t,0.0d0)                        !upper diagonale entries
              enddo
           else                                                        !if epsilon matrix is a 2x2 matrix
              epsmat(1,1)=dcmplx(-2.0d0*t*dcos(kx),0.0d0)              !top left entry
              epsmat(1,2)=dcmplx(-t-t*dcos(dfloat(L)*ky),0.0d0)+   &    !top right entry
                  dcmplx(0.0d0,-t*dsin(dfloat(L)*ky))
              epsmat(2,2)=dcmplx(-2.0d0*t*dcos(kx+B*2.0d0*Pi),0.0d0)   !bottom right entry
              epsmat(2,1)=dcmplx(-t-t*dcos(dfloat(L)*ky),0.0d0)-   &   !bottom left entry
                  dcmplx(0.0d0,-t*dsin(dfloat(L)*ky))
           endif

           do i=omoffset+1,omoffset+omnumber         
              do j=1,L
                 do k=1,L
                    ginv(j,k)=-epsmat(j,k)                             !do we need i loop here? Still correct but needless but you need less memory
                 enddo
              enddo
              do j=1,L
                 ginv(j,j)=ginv(j,j)+W(i)                              ![-\epsilon(k)_{j,j} + i\nu + \mu - \Sigma] = G_{j,j}^{-1}(\nu, k), but have \epsilon_{l,lPrime} outside diagonale as well
              enddo
           
              CALL ZGETRF(L,L,ginv,L,ipiv,infoinv)                     !Lapack routine for "These subroutines factor general matrix A using Gaussian elimination with partial pivoting"
              CALL ZGETRI(L,ginv,L,ipiv,workinv,10*L,infoinv)          !Lapack routine to "ZGETRI computes the inverse of a matrix using the LU factorization computed by ZGETRF" Now ginv is G(\nu, k), inverts ginv
              
              do j=1,L
                 do k=1,L
                    gloc(i,j,k)=gloc(i,j,k)+ginv(j,k)                  !Now we get Glocal_{j,k}(omega) = sum over k (G_{l,lPrime} (\nu, k))
                 enddo
              enddo
           enddo
        enddo
     enddo 

     !checkone=checkone*dfloat(L)/(4.0d0*dfloat(ksteps**2))
     checkone=checkone/((2.0d0*dfloat(ksteps))*(2.0d0*dfloat(kstepsY)))                     !if this is not 1.0 or 10^(-15) or so away, then the normalization is wrong
     if (omoffset.le.2) then
        write(6,*)"Check momentum sum: ",checkone                                           !should be exactly 1
     endif

     do i=omoffset+1,omoffset+omnumber                                                      !basically from -omnumber/2 to omnumber/2
        do j=1,L
           do k=1,L
              !gloc(i,j,k)=gloc(i,j,k)*dfloat(L)/(4.0d0*dfloat(ksteps**2))                  this is the wrong normalization, since it does not account for the integer division when ksteps/L
              gloc(i,j,k)=gloc(i,j,k)/(4.0d0*dfloat(ksteps)*dfloat(kstepsY))                                   !better do this where we consider integer division by using kstepsY
           enddo
        enddo
     enddo
  
     if (omoffset.le.2) then
        do j=1,L
           do i=omoffset+1,omoffset+omnumber
              write(1000+j,*)om(i),                                 &
                  dreal(gloc(i,j,j)),dimag(gloc(i,j,j))
           enddo
        enddo
     endif


     return 
     end


subroutine selfconst_fcc(omnumber,omoffset,Iwmax,ksteps, &
                      beta,xmu,om,Energya,g0wand,Gw,Gloc)
      
    integer,     intent(in)  :: omnumber,omoffset,Iwmax,ksteps
    real(dp),    intent(in)  :: beta,xmu
    real(dp),    intent(in)  :: om(0:Iwmax),Energya(0:ksteps,0:ksteps,0:ksteps)
    complex(dp), intent(inout) :: g0wand(0:Iwmax),Gw(0:Iwmax)
    complex(dp), intent(out) :: Gloc(0:Iwmax)
    integer i,ikx,iky,ikz
    real(dp) weightx,weighty,weightz,weightm
    complex(dp) W(omoffset+1:omoffset+omnumber)
    complex(dp) gand(omoffset+1:omoffset+omnumber)

    do i=omoffset+1,omoffset+omnumber
        gand(i)=1.d0/g0wand(i)
        Gw(i)=1.d0/Gw(i)
        W(i)= (0.d0,1.d0)*om(i)-gand(i)+Gw(i)+xmu!-tpri
    enddo 
    
    do i=omoffset+1,omoffset+omnumber
        Gloc(i)=dcmplx(0.d0,0.d0)
        do ikx=0, ksteps-1
            do iky=0, ksteps-1
                do ikz=0, ksteps-1
                    Gloc(i)=Gloc(i)+1.d0/(W(i) - Energya(ikx,iky,ikz))
                enddo
            enddo
        enddo
        Gloc(i) = Gloc(i)/(dfloat(ksteps)**3)
        Gloc(i)=1.d0/Gloc(i)+gand(i)-Gw(i)
        Gloc(i)=1.0d0/Gloc(i)   !without damping
    enddo
  
    return
 end


end
