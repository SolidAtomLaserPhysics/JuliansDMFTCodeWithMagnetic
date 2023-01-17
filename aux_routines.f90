module types
    use ISO_Fortran_env
    integer, parameter :: dp = real64
    integer, parameter :: qp = REAL128
    integer, parameter :: id = int64
end module

module AuxRoutines
    use types
    implicit none
    contains

subroutine rheader(fh)
    implicit none
    integer i,fh
    do i=1,9
        read(fh,*)
    enddo
end 

subroutine read_andpar(hubb_andpar_fh,ns,epsk,tpar,xmu)
    implicit none
    integer,  intent(in) :: hubb_andpar_fh, ns
    real(dp), intent(out) :: xmu
    real(dp), intent(out) :: epsk(ns),tpar(ns)
    integer i

    rewind(hubb_andpar_fh)
    call rheader(hubb_andpar_fh)
    do i=2,ns
       read(hubb_andpar_fh,*) epsk(i)
    end do
    read(hubb_andpar_fh,*)
    do i=1,ns-1
       read(hubb_andpar_fh,*) tpar(i)
    end do
    read(hubb_andpar_fh,*) xmu
    rewind(hubb_andpar_fh)
end

subroutine initial(hubb_andpar_fh,epsk,tpar,ns)
    implicit none
    integer, intent(in)  :: hubb_andpar_fh,ns
    real(dp) , intent(out) :: epsk(ns),tpar(ns)
    integer i
    call rheader(hubb_andpar_fh)
    do i=2,ns
       read(hubb_andpar_fh,*) epsk(i)
    end do
    read(hubb_andpar_fh,*)
    do i=1,ns-1
       read(hubb_andpar_fh,*) tpar(i)
    end do
end 

subroutine write_andpar(hubb_andpar_fh,ns,epsk,tpar,xmu,Iwmax)
    implicit none
    integer,  intent(in) :: hubb_andpar_fh,ns,Iwmax
    real(dp), intent(in) :: xmu
    real(dp), intent(in) :: epsk(ns),tpar(ns)
    integer i

    write(*,*)'------------------------------------------------'
    write(*,*) 'new Anderson parameters '
    write(*,*)'------------------------------------------------'
    write(*,*) 'Eps(k) '
    do i=2,ns
        write(*,'(2f27.16)')epsk(i)
    end do
    write(*,*)'V(k) '
    do i=1,ns-1
        write(*,'(2f27.16)')tpar(i)
    end do
    write(*,'(f27.16,"   #chemical potential")')xmu
    call flush(6)            

    call rheader(hubb_andpar_fh)
    do i=2,ns
        write(hubb_andpar_fh,*) epsk(i)
    end do
    write(hubb_andpar_fh,'(9a)')' tpar(k)'
    do i=1,ns-1
        write(hubb_andpar_fh,*) tpar(i)
    end do
    write(hubb_andpar_fh,*) xmu,'    #chemical potential'
    call flush(hubb_andpar_fh)

end

subroutine datain(fh,fh2,ns,uhub,xmu,hmag,beta,wlow,wup,deltino,imaxmu, &
                  deltamu,iteramax,testdmft,dens_wanted,ifix, &
                  inew,iexp,th0,iauto)
    implicit none
    integer,  intent(in) :: fh,fh2,ns
    real(dp), intent(out) :: uhub,xmu,hmag,beta,deltamu,dens_wanted
    real(dp), intent(out) :: wlow,wup,deltino,testdmft,th0
    integer,  intent(out) :: imaxmu,iteramax,ifix,inew,iexp,iauto
    integer ns_in,i
    character(len=80) text
    !hubb.andpar
    read(fh,'(a)') text
    read(fh,*) uhub,hmag
    read(fh,'(a)') text
    read(fh,*) beta, wlow, wup, deltino
    read(fh,'(a)') text
    read(fh,*) ns_in,imaxmu,deltamu,iteramax,testdmft   
    read(fh,'(a)') text
    read(fh,*)ifix,dens_wanted,inew,iauto
    read(fh,'(a)') text
    read(fh,*)th0,iexp
    write(*,*)'number of conduction bath levels =',ns-1
    write(*,*)'WARNING: ns given in hubb.dat is NOT used!!'
    write(*,*)'number of iterations :',iteramax
    call rheader(fh2)
    do i=2,ns
        read(fh2,*)
    end do
    read(fh2,*)
    do i=1,ns-1
        read(fh2,*)
    end do
    read(fh2,*) xmu 
    rewind(fh)
    rewind(fh2)
    return
end

pure integer function nchoos(n1,n2)
    implicit none
    integer, intent(in)  :: n1,n2
    integer i
    real(dp)  xh

    xh = 1.0
    if(n2.lt.0) then
      nchoos = 0
      return
    endif
    if(n2.eq.0) then
      nchoos = 1
      return
    endif
    do i = 1, n2
      xh = xh * float(n1+1-i)/float(i)
    enddo
    nchoos = int(xh + 0.5)
end

subroutine computematrix(ioffset,nblock,ntot,b,ns,abas,nsp1,  &
         nmax,nleng,imat,idmat,iddomat,ntable,idouble,myid)
!     this is meant to be a subroutine that generates
!     1) the blocks
!     2) the non-diagonal matrix elements within each block (imat)
!     3) the matrix elements of d^{dagger} -> idmat
    implicit none
    integer, intent(in) :: nblock,ntot,ns,nmax,nsp1,myid,b(ns)
    integer, intent(out) :: ioffset(nblock+1),imat(nmax,nmax,nblock)
    integer, intent(out) :: idmat(ntot),iddomat(ntot),idouble(ntot)
    integer, intent(out) :: abas(ns),nleng(nblock),ntable(ntot)
    integer i,j,nup,ndo
    
!     compute ioffsets
    ioffset(1)=0
    ioffset(2)=1
    do i=2,nblock-1
       nup=mod(i-1,nsp1)
       ndo=(i-1)/nsp1
       ioffset(i+1)=ioffset(i)+nchoos(ns,nup)*nchoos(ns,ndo)
    enddo
    ioffset(nblock+1)=ntot
    call buildblocks(ns,ntot,nblock,ntable,ioffset,nleng,  &
         b,abas,nsp1)
!     loop over blocks -> compute the Hamiltonian
    idouble=0
    do i=1,nblock
       do j=1,nleng(i)
          call findnupndo(ntable(ioffset(i)+j)-1,nup,ndo,b,  &
               abas,ns)
          if (abas(1).eq.2) idouble(ioffset(i)+j)=1
       enddo
       call hamilt(ntable(ioffset(i)+1),nleng(i),  &
            imat(1,1,i),b,ns,abas,nmax)
       if (mod(i-1,nsp1).ne.ns) then
          call ddag(ntable(ioffset(i)+1),ntable(ioffset(i+1)+1),  &
        nleng(i),nleng(i+1),idmat(ioffset(i)+1),b,ns,abas)
       endif
       if(i.gt.(ns+1)) then
          call ddown(ntable(ioffset(i)+1),ntable(ioffset(i-ns-1)+1),  &
       nleng(i),nleng(i-ns-1),iddomat(ioffset(i)+1),b,ns,abas)
       endif
    enddo
    if (myid.eq.0) then
       write(*,*)'# of independent blocks:',nblock
    endif
    return
end

subroutine hamilt(ntablein,nlen,imatin,b,ns,abas,nmax)
    implicit none
    integer, intent(in)  :: ntablein(nlen),nlen
    integer, intent(out) :: imatin(nmax,nmax)
    integer, intent(in)  :: b(ns),ns,nmax
    integer abas(ns),i,j,nup,ndo,ntnew,inew,k,isegno

    imatin = 0
!     diagonal part is already coded in ntab!!!!
!     non-diagonal part
    do i=1,nlen
       call findnupndo(ntablein(i)-1,nup,ndo,b,abas,ns)
       if (abas(1).gt.0) then
!     loop over final sites
          do j=2,ns
             if (abas(j).le.0) then
      ntnew=ntablein(i)+b(j)*(1-2*abas(j))+b(1)*(1-2*abas(1))
                call findinew(ntnew,ntablein,nlen,inew)
                nup=0
                do k=2,j-1
                   if (abas(k).gt.0) nup=nup+1
                enddo
                isegno= 1-2*mod(nup,2)
                imatin(i,inew)=j*isegno
                imatin(inew,i)=j*isegno
             endif
          enddo
       endif
!     if down is on site 1, it can hop
    if (abas(1).eq.2.or.abas(1).eq.-1) then
       do j=2,ns
          if (abas(j).eq.0.or.abas(j).eq.1) then
      ntnew=ntablein(i)+b(j)*(2*abas(j)-1)-b(1)*abas(1)/abs(abas(1))
             call findinew(ntnew,ntablein,nlen,inew)
             ndo=0
             do k=2,j-1
                if (abas(k).eq.2.or.abas(k).eq.-1) ndo=ndo+1
             enddo
             isegno=1-2*mod(ndo,2)
             imatin(i,inew)=j*isegno
             imatin(inew,i)=j*isegno                  
          endif
       enddo
    endif
    enddo
    return
end

subroutine findinew(ntnew,ntab,nlen,inew)
    implicit none
    integer nlen,i,ntnew,inew
    integer ntab(nlen)
    do i=1,nlen
       if (ntab(i).eq.ntnew) goto 777
    enddo
 777  inew=i
    return
end

subroutine buildblocks(ns,ntot,nblock,ntable,ioffset,nfound,  &
         b,abas,nsp1)
    implicit none
    integer,intent(in) ::ns,ntot,nblock,ioffset(nblock+1),b(ns),nsp1
    integer,intent(out) :: ntable(ntot),nfound(nblock),abas(ns)
    integer nup,ndo,indblock,i

    nfound = 0
    do i=1,ntot
       call findnupndo(i-1,nup,ndo,b,abas,ns)
       indblock=nup+nsp1*ndo+1
       nfound(indblock)=nfound(indblock)+1

       if (nfound(indblock).gt.(ioffset(indblock+1)-  &
            ioffset(indblock)))  &
            then
          write(*,*)'Il blocco ',indblock,  &
               ' ', nfound(indblock), ' should be < ',   &
              ioffset(indblock+1),'-',ioffset(indblock)
          write(*,*)'nup, ndo =',nup,ndo
          error stop
       endif
       ntable(nfound(indblock)+ioffset(indblock))=i
    enddo
    return
end

subroutine findnupndo(nstato,nup,ndo,b,abas,ns)
    implicit none
    integer,intent(in)  :: ns, b(ns),nstato
    integer,intent(out) :: nup,ndo,abas(ns)
    integer nb,i,is

    is = nstato
    nup=0
    ndo=0
    do 1 i = ns,1,-1
      nb = is/b(i)
      abas(i) = nb - 1
      is = is - nb * b(i)
      if (abas(i).gt.0) nup=nup+1
      if (abas(i).eq.2.or.abas(i).eq.-1) ndo=ndo+1
  1 continue
    
    return
end

subroutine ddag(ntab0,ntab1,nlen0,nlen1,idmatin,b,ns,abas)
    implicit none

    integer nlen0,nlen1,ns,j,iold,nup,ndo,inew,k
    integer ntab0(nlen0)
    integer ntab1(nlen1)
    integer b(ns)
    integer abas(ns)
    integer idmatin(nlen0)
    do j=1,nlen0
       iold=ntab0(j)
       call findnupndo(iold-1,nup,ndo,b,abas,ns)
       if (abas(1).le.0) then
          inew=iold+(1-2*abas(1))
          do k=1,nlen1
             if (ntab1(k).eq.inew) then
                idmatin(j)=k
                goto 555
             endif
          enddo
 555        continue
       else
          idmatin(j)=0
       endif
    enddo
    return
end

subroutine ddown(ntab0,ntab1,nlen0,nlen1,iddomatin,b,ns,abas)
    implicit none
    integer nlen0,nlen1,ns,j,iold,nup,ndo,inew,k,isegno
    integer ntab0(nlen0)
    integer ntab1(nlen1)
    integer b(ns)
    integer abas(ns)
    integer iddomatin(nlen0)
    do j=1,nlen0
       iold=ntab0(j)
       call findnupndo(iold-1,nup,ndo,b,abas,ns)

!ccc     If down is on site 1
       if (abas(1).eq.-1.or.abas(1).eq.2) then
          inew=iold + (1-2*abas(1))/3
!c          ddown operator has to hop nup spin to reach his spin down
          isegno=1-2*mod(nup,2)
     
          do k=1,nlen1
             if (ntab1(k).eq.inew) then
                iddomatin(j)=k*isegno
                goto 555
             endif
          enddo
 555        continue
       else
          iddomatin(j)=0
       endif
    enddo
    return
end


! =============== DMFT Routines =================

subroutine difference(nbparm,x,f,epsk,tpar,ns,piob,Iwmax,G0w,G0wand,om,symm)
    implicit none
    logical,    intent(in) :: symm
    integer,    intent(in) :: ns,Iwmax,nbparm
    real(dp),   intent(in) :: piob,om(0:Iwmax)
    complex(dp),intent(in) :: G0w(0:Iwmax)
    real(dp),   intent(inout) :: x(nbparm)
    real(dp),   intent(out):: tpar(ns),epsk(ns),f
    complex(dp),intent(out):: G0wand(0:Iwmax)
    integer                :: icount,i,nsym
    real(dp)               :: diff,norm
    complex(dp)            :: cdummy1

    icount=1

    if(symm) then
        nsym=ns/2
        do i=1,nsym
            icount=icount+1
            Epsk(icount)=x(i)
        end do
        icount=0
        do i=nsym+1,2*nsym
            icount=icount+1
            tpar(icount)=x(i)
        end do
        do i=nsym+2,ns
            Epsk(i)=-Epsk(i-nsym)
            tpar(i-1)=tpar(i-nsym-1)
        enddo 
        Epsk(nsym+1)=Epsk(nsym+1)*dfloat(ns-2*(ns/2))  ! Sets epsilon(#bathsites/2+1)=0 for odd number of bath sites
        norm=0.d0
    else
        do i=1,ns-1
            icount=icount+1
            Epsk(icount)=x(i)
        end do
        icount=0
        do i=ns,2*(ns-1)
            icount=icount+1
            tpar(icount)=x(i)
        end do
        ! norm=0.d0
        ! do i=1,ns-1
        !     norm=norm+tpar(i)**2
        ! enddo
        ! norm=dsqrt(0.25d0/norm)
        ! do i=1,ns-1
        !     tpar(i)=norm*tpar(i)
        ! enddo
    endif

    diff=0.d0
    do i=0,Iwmax
        call calcg0(om(i),cdummy1,tpar,epsk,ns,Iwmax)
        g0wand(i)=cdummy1
        diff = diff + abs(g0w(i)-g0wand(i))/dfloat(i+1)
    end do
    f=diff/dfloat(Iwmax+1)
       
    if(symm) then
        icount=1
        do i=1,nsym
            icount=icount+1
            x(i)=Epsk(icount)
        end do
        icount=0
        do i=nsym+1,2*nsym
            icount=icount+1
            x(i)=tpar(icount)
        end do
    else
        icount=1
        do i=1,ns-1
            icount=icount+1
            x(i)=Epsk(icount)
        enddo
        icount=0
        do i=ns,2*(ns-1)
            icount=icount+1
            x(i)=tpar(icount)
        enddo
    endif
end

end
