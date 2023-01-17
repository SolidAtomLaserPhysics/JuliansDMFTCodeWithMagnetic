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
