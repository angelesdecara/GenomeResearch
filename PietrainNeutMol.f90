! This program is to manage Pietrain using coancestry measured as IBS averaged over individuals
!
!
program ManagePietrain
  implicit none
   integer :: i,j,k,l, icrom, igen, irep, n, iguy,igamete,gametetemp,idummy,io,imarker
   integer, parameter :: dp=KIND(1.d0)
   integer, parameter :: nindbase=156,nmarkers=49803,ncrom=18
   integer, parameter :: nfem=23, nmale=23, nind=nfem+nmale, ngen=10, nrep=1000
   integer :: genotbase(nindbase,nmarkers,2),genotBaseHere(nind,nmarkers,2)
   integer :: genot(nind,nmarkers,2),genottemp(nind,nmarkers,2),chr(nmarkers)
   integer::  pos(nmarkers), ncross, geneal(nind,2),crosses(100)
   integer :: sol(nind,2),dummy, NumberMarkers(ncrom),nbins
   real(dp):: fmol(nind,nind), recombrate(ncrom), poisson, meandummy
   real(dp):: recomb(ncrom)!!=.362726d0 ! cM per Mb
   real (dp):: meanOH(0:ngen),varOH(0:ngen),meanGD(0:ngen),varGD(0:ngen)
   real(dp):: meanmalecontrib(0:ngen),meanfemalecontrib(0:ngen),varfemalecontrib(0:ngen)
   real(dp):: varmalecontrib(0:ngen), meanfixedpos(0:ngen),varfixedpos(0:ngen)
   real(dp) :: meanmaxcontrib(0:ngen),meanmincontrib(0:ngen)
   real(dp)::oh,gd,fixpos, x, sumx
   integer:: aux(nind)
   character(len=6)::a
   character(len=2)::base
   character(len=50)::fileHaps,markername, indName,sampleName
   real:: ranmar
   call rmarin(5879,6958)
   genot=99; genottemp=99

!
! read phased genotypes data from chromosome files
!
   imarker=0; NumberMarkers=0
   do icrom=1,ncrom
      io=0
      write(fileHaps,*)icrom
      fileHaps='PietrainNeut/chr'//trim(adjustl(fileHaps))//'.phased.haps'
      open(53,file=trim(adjustl(fileHaps)))
      do while(io.eq.0)
         imarker=imarker+1
!         write(6,*)io
         read(53,*,IOSTAT=io)idummy,markername,pos(imarker),base,base,((genotbase(iguy,imarker,igamete),igamete=1,2),iguy=1,nindbase)
         chr(imarker)=icrom
      enddo
      close(53)
      imarker=imarker-1
      if(icrom.gt.1)then
         NumberMarkers(icrom)=imarker-sum(NumberMarkers(1:icrom-1))
      else if(icrom.eq.1)then
         NumberMarkers(icrom)=imarker
      endif
!      write(6,*)'icrom,nmarkers',icrom,imarker,NumberMarkers(icrom)
   enddo

!
! take the individuals I want from these genotypes
!
   open(21,file='PietrainNeut/Individuals.dat')

   do i=1,nind
      read(21,*)a,indName,j,j,j,j,j
      open(22,file='PietrainNeut/chr1.phased.sample')
      read(22,*)
      read(22,*)
      do iguy=1,nindbase
         read(22,*)a,sampleName,j,j,j,j,j
         if(sampleName.eq.indName)then
!            write(6,*)'sampling guys',sampleName,indName
            genotBaseHere(i,:,:)=genotbase(iguy,:,:)
            close(22)
            exit
         endif
      enddo
   enddo
   close(21)

!
! recomb rate per chromosome
   do icrom=1,ncrom
      sumx=0.d0
      open(23,file='rec_rate.csv')
      read(23,*)
      io=0
      do while(io.eq.0)
         read(23,*,IOSTAT=io)i,n,x
         if(i.eq.icrom)then
            sumx=sumx+x
            nbins=n
     !       write(6,*)icrom,i,n,x,sumx
         endif
      enddo
      close(23)
      recomb(icrom)=sumx/dble(nbins+1)
      recombrate(icrom)=recomb(icrom)*(pos(sum(NumberMarkers(1:icrom)))-pos(1+sum(NumberMarkers(1:icrom-1))))/1.d8
!      write(6,*)icrom,recombrate(icrom),sumx
   enddo

!   recombrate=recomb*dble(pos(nmarkers)-pos(1))/1.d8 ! divided by 10^6 for Mb and 10^2 for cM to Morgan

! everything 0.d0
   meanOH=0.d0;varOH=0.d0; meanGD=0.d0; varGD=0.d0
   meanmalecontrib=0.d0;meanfemalecontrib=0.d0;
   varmalecontrib=0.d0;varfemalecontrib=0.d0;
   meanfixedpos =0.d0; varfixedpos=0.d0
   meanmaxcontrib=0.d0; meanmincontrib=0.d0


   do irep=1,nrep

      genot=genotBaseHere
!
! at t=0 meanOH and so on
! 
      call heter(genot,nind,nmarkers,oh,gd,fixpos)
      meanOH(0)=meanOH(0)+oh
      varOH(0)=varOH(0)+oh**2
      meanGD(0)=meanGD(0)+gd
      varGD(0)=varGD(0)+gd**2
      meanfixedpos(0)=meanfixedpos(0)+fixpos
      varfixedpos(0)=varfixedpos(0)+fixpos**2

      write(6,*)irep,meanGD(0)/irep,gd
!
! management during ngen generations
!
      do igen=1,ngen
!
! we need the coancestries to later calculate how many offspring each guy leaves
! ibdsize is the minimum size for a segment
!
         call molcoanc(nind,nmarkers,genot,fmol)
         write(6,*)'averfmol',igen,sum(fmol(1:nind,1:nind))/dble(nind**2)
!
! we use fmol in the annealing
!
 !        write(6,*)'before annealing',sum(fibd(1:nind,1:nind))/dble(nind**2)
         call modannealing(fmol,nfem,nmale,sol)
!
! sol(nind,2), sol(i,1) is the number of sons that individual i leaves
!         sol(i,2) is the number of daughters that individual i leaves
!
! variance in contribs each gen

!         write(6,*)'sum contribs',sum(sol(1+nmale:nind,1)),sum(sol(1+nmale:nind,2))

         meandummy=0.d0
         do idummy=1,nmale
            meandummy=sol(idummy,1)+sol(idummy,2)
         enddo
         meandummy=meandummy/dble(nmale)
         meanmalecontrib(igen)=meanmalecontrib(igen)+meandummy
         varmalecontrib(igen)=varmalecontrib(igen)+meandummy**2
!
         meandummy=0.d0
         do idummy=1+nmale,nind
            meandummy=sol(idummy,1)+sol(idummy,2)
!            write(6,*)idummy,meandummy
         enddo
         meandummy=meandummy/dble(nfem)
         meanfemalecontrib(igen)=meanfemalecontrib(igen)+meandummy
         varfemalecontrib(igen)=varfemalecontrib(igen)+meandummy**2

  ! fathers' genealogies
         dummy=0; aux=0
         do n=1,2
            do i=1,nmale
               if(sol(i,n).gt.0)then
                  do j=1,sol(i,n)
                     dummy=dummy+1
                     aux(dummy)=i
!                     geneal(dummy,1)=i ! geneal(i,1) is father of i
                  enddo
               endif
            enddo
         enddo
         call desordena(aux(1:nmale),nmale)
         call desordena(aux(1+nmale:nind),nfem)
         do i=1,nind
            geneal(i,1)=aux(i)
         enddo
  ! mothers' genealogies
         dummy=0; aux = 0
         do n=1,2
            do i=nmale+1,nind
               if(sol(i,n).gt.0)then
                  do j=1,sol(i,n)
                     dummy=dummy+1
                     aux(dummy)=i
!                     geneal(dummy,2)=i ! geneal(i,2) is mother of i
                  enddo
               endif
            enddo
         enddo
         call desordena(aux(1:nmale),nmale)
         call desordena(aux(1+nmale:nind),nfem)
         do i=1,nind
            geneal(i,2)=aux(i)
         enddo
!
! generate new pop
!
         do i=1,nind
            do icrom=1,ncrom

               do k=1,2
!
! crossovers for this guy and this gamete
!
                  ncross=poisson(recombrate(icrom))
                  do l=1,ncross
                     crosses(l)=ranmar()*(sum(NumberMarkers(1:icrom))-sum(NumberMarkers(1:icrom-1)))+sum(NumberMarkers(1:icrom-1))+1
!*(pos(sum(NumberMarkers(1:icrom)))-pos(1+sum(NumberMarkers(1:icrom-1))))+1
!(pos(nmarkers)-pos(1))+1
                  enddo
                  call sortincrease(crosses,ncross)
                  crosses(ncross+1)=sum(NumberMarkers(1:icrom))
!
!
                  gametetemp=ranmar()*2+1
                  genottemp(i,1+sum(NumberMarkers(1:icrom-1)):sum(NumberMarkers(1:icrom)),k)=&
                       genot(geneal(i,k),1+sum(NumberMarkers(1:icrom-1)):sum(NumberMarkers(1:icrom)),gametetemp)

                  do n=1,ncross+1
                     gametetemp=3-gametetemp
                     genottemp(i,crosses(n)+1:crosses(n+1),k)=genot(geneal(i,k),crosses(n)+1:crosses(n+1),gametetemp)
                  enddo

               enddo
            enddo
         enddo
         genot=genottemp
!
         oh=0.d0;gd=0.d0;fixpos=0.d0
         call heter(genot,nind,nmarkers,oh,gd,fixpos)
         meanOH(igen)=meanOH(igen)+oh
         varOH(igen)=varOH(igen)+oh**2
         meanGD(igen)=meanGD(igen)+gd
         varGD(igen)=varGD(igen)+gd**2
         meanfixedpos(igen)=meanfixedpos(igen)+fixpos
         varfixedpos(igen)=varfixedpos(igen)+fixpos**2

      enddo
!
!
!
      open(54,file='PietrainNeut/manMolPI.dat')
      write(54,'(a120,i4)')'#igen,meanOH,varOH,meanGD,varGD,meanfixedpos,varfixedpos,&
           meanfemalecontrib,varfemalecontrib,meanmalecontrib,varmalecontrib',irep
      do igen=0,ngen
         write(54,'(i4,10f18.14)')igen,meanOH(igen)/dble(irep),varOH(igen)/dble(irep)-(meanOH(igen)/dble(irep))**2,&
              meanGD(igen)/dble(irep),varGD(igen)/dble(irep)-(meanGD(igen)/dble(irep))**2,&
              meanfixedpos(igen)/dble(irep),varfixedpos(igen)/dble(irep)-(meanfixedpos(igen)/dble(irep))**2,&
              meanfemalecontrib(igen)/dble(irep),varfemalecontrib(igen)/dble(irep)-(meanfemalecontrib(igen)/dble(irep))**2,&
              meanmalecontrib(igen)/dble(irep),varmalecontrib(igen)/dble(irep)-(meanmalecontrib(igen)/dble(irep))**2
      enddo
      close(54)
   enddo

 end program ManagePietrain
!
! sorting in increasing order the recomb points
! poisson
! ranmar
! fmol for nind,nind 
! annealing for different number of males and females
!------
subroutine desordena(ix,nx)
  implicit none
  integer::i,n,nx,ix(nx),ic
  real::ranmar
  do i=1,nx
     n=ranmar()*nx+1
     ic=ix(i)
     ix(i)=ix(n)
     ix(n)=ic
  enddo
end subroutine desordena
!---------- oh, gd, 
 subroutine heter(genot,nind,nmarkers,oh,gd,fixpos)
   implicit none
   integer, parameter :: dp=KIND(1.d0)
   integer :: nind,nmarkers,genot(nind,nmarkers,2)
   real(dp)::oh,gd,fixpos,p,q
   integer :: i,j,k,l

   p=0.d0; q=0.d0; gd=0.d0;fixpos=0.d0
   do i=1,nmarkers
      p=count(genot(:,i,:)==0)/dble(nind*2)
      q=count(genot(:,i,:)==1)/dble(nind*2)
      gd=gd+(1.d0-p**2-q**2)         
      if((p.eq.1).or.(q.eq.1.d0))fixpos=fixpos+1.d0
   enddo
   gd=gd/dble(nmarkers)
   fixpos=fixpos/dble(nmarkers)

   oh=0.d0
   do i=1,nind
      do l=1,nmarkers
         if(genot(i,l,1).ne.genot(i,l,2))oh=oh+1.d0
      enddo
   enddo
   oh=oh/dble(nind*nmarkers)

 end subroutine heter
!------ subroutine to order array in increasing order --
subroutine sortincrease(vectortosort,ncross)
  implicit none
  integer :: ncross, i,j, temp
  integer:: vectortosort(ncross), vectorlength

  vectorlength=ncross
  do i=1,vectorlength-1
     do j=i,vectorlength
        if (vectortosort(j).lt.vectortosort(i))then
           temp=vectortosort(i)
           vectortosort(i)=vectortosort(j)
           vectortosort(j)=temp
        endif
     enddo
  enddo
end subroutine sortincrease
!-----------------------------------------------
!--------------
function poisson(media)
  implicit none
  integer :: k, poisson, sem1,sem2
  real :: ranmar
  double precision :: L, p, media

  L=exp(-media)
  k=0
  p=1
  do while(p>L)
     k=k+1
     p=p*dble(ranmar())
  enddo

  poisson=k-1
end function poisson
!====================
!
subroutine molcoanc(nind,nmarkers,genot,fmol)
  implicit none
  integer, parameter :: dp=KIND(1.d0)
  integer :: i,j,k,l,nind,nmarkers,i1,i2,j1,j2, cont
  integer :: genot(nind,nmarkers,2)
  real (dp) :: fmol(nind,nind)

  fmol = 0.d0
  do i=1,nind
     do j=1,nind
        cont=0
        do i1=1,2
           do j1=1,2
              do l=1,nmarkers
                 if (genot(i,l,i1).eq.genot(j,l,j1))then
                    cont=cont+1
                 endif
!
              enddo
           enddo
!
!        write(66,*)frohman(i,j)
        enddo
        fmol(i,j)=dble(cont)/dble(4*nmarkers)
     enddo
  enddo

end subroutine molcoanc
!-----------------------------------
! annealing!
!
subroutine modannealing(fibd,nfem,nmale,sol)
  implicit none
  integer, parameter :: dp=KIND(1.d0)
  integer :: nfem, nmale, m, n,nind, nchanges
  integer :: sol(nfem+nmale,2),sola(nfem+nmale,2),solopt(nfem+nmale,2)
  real (dp) :: kt,fibd(nfem+nmale,nfem+nmale)
  real (dp) :: camb, eneal, omega, eneopt, eneac, delta,ene
  integer :: ch,pos,anc, niv,ivig,rep,l,i,j
  real ranmar
  real (dp), parameter::t=.01d0,k=.9d0
  character(len=10)::dia,hora,zona
  integer ::valores(8)
  call date_and_time(dia,hora,zona,valores)
  call rmarin(valores(5)*11+43,valores(8)+11)

!
! start
!
  nind=nfem+nmale
  sol=0; sola=0;solopt=0
! random solution
  do m=1,nmale
     n=ranmar()*nmale+1
     sol(n,1)=sol(n,1)+1
     n=ranmar()*nfem+nmale+1
     sol(n,1)=sol(n,1)+1
  enddo
  do m=1,nfem
     n=ranmar()*nmale+1
     sol(n,2)=sol(n,2)+1
     n=ranmar()*nfem+nmale+1
     sol(n,2)=sol(n,2)+1
  enddo
!
  kt=t/k

!first test
!  write(6,*)'in annealing',sum(sol(1:nmale,1)),sum(sol(1:nmale,2)),sum(sol(nmale+1:nind,1))

  call energ(sol,fibd,nfem,nmale,ene)
!  write(6,*)'in annealing sol 0',ene
  solopt=sol
  eneopt=ene
  ivig=5000
! Levels loop
!  
  nchanges = 0
  do niv=1,200
     camb=1+(4*ivig/5000)
     kt=kt*k
     ivig=0
! try many changes
     
     do rep=1,5000
! alternative solution
! and copy the current solution into the alternative
        sola=sol
        do l=1,int(camb)
! A guy loses a contribution
           i=int(ranmar()*2)+1
20         m=int(ranmar()*nind)+1
          !! write(6,*)'l287',i,m,sola(m,i)
           if(sola(m,i).eq.0) goto 20
           sola(m,i)=sola(m,i)-1
! Someone else "gains" it
           pos=0
           anc=nmale
           if(m.gt.nmale) then
              pos=nmale
              anc=nfem
           endif
30         n=int(ranmar()*anc)+1+pos
           if(n.eq.m) goto 30
           sola(n,i)=sola(n,i)+1
        enddo
!
! check energy of this alternative state
!
        call energ(sola,fibd,nfem,nmale,eneal)
!        write(6,*)'altern',eneal,ene
! change?
        ch=0
        delta=dmax1(eneal-ene,0.)
        omega=exp(-delta/kt)
        if(omega.ge.1.d0) then
           ch=1
           nchanges=nchanges+1
           if(eneal.lt.eneopt) then
              eneopt=eneal
              solopt=sola
           endif
        else
           if(dble(ranmar()).lt.omega) then
              ch=1
              nchanges=nchanges+1
           endif
        endif
! if ch=1, change
        if (ch.eq.1)then
           sol=sola
           ene=eneal
           ivig=ivig+1
        endif
     enddo
!
! stops changing , go back to main program
     if(ivig.eq.0)then
        return
     endif
  enddo
  write(6,*)'all loops done'!,nchanges,ch,delta,ene,eneal !! check how many times it enters the change to see whether it actually changes or not

end subroutine modannealing
!
!--- energy
subroutine energ(sol,fibd,nfem,nmale,ene)
  implicit none
  integer, parameter :: dp=KIND(1.d0)
  integer ::i,k,j,nfem,nmale,nind
  integer :: sol(nfem+nmale,2), weight(2)
  real(dp):: ene, fibd(nfem+nmale,nfem+nmale), prodfsol(nfem+nmale)


  weight(1)=4*(nmale**2)
  weight(2)=4*(nfem**2)

  ene =0.d0; prodfsol=0.d0
  do i=1,2
     prodfsol=matmul(dble(sol(:,i)),fibd)
     ene=ene+dot_product(dble(sol(:,i)),prodfsol)/weight(i)
  enddo
  prodfsol=matmul(dble(sol(:,1)),fibd)
  ene=ene+dot_product(dble(sol(:,2)),prodfsol)/(2*nfem*nmale)
  ene=ene/4.d0

end subroutine energ
!-------------------
!----------------------------------------------
! ranmar y rmarin
!---------------------------------------------
! I received this from STUART@ADS.COM,
! "Toward a Universal Random Number Generator" by George Marsaglia and Arif
! Zaman. Florida State University Report: FSU-SCRI-87-50 (1987).
! It was later modified by F. James and  published in "A Review of
! Pseudo-random Number Generators"
! Stuart says this is the BEST KNOWN RANDOM NUMBER GENERATOR available.
! It passes all tests for a random number generator, and has a period
! of 2^144, and will give bit-identical results on all machines with
! at least 24-bit manitissas in the floating point representation.
! The algorithm is a combination of Fibonacci sequence (with lags of 97
! and 33, and operation "subtraction plus one, modulo one") and an
! "arithmetic sequence" (using subtraction.
! There are three routines contained herein. I have separated them and the
! notes with double lines of asterisks (just delete all between them and the
! lines of asterisks).
! The first program, called TstRAN is a test for correct operation.
! The second part, RANMAR, is the function itself, which must be included
! in any program where you wish to call it from.
! The third part is the initialization routine, which asks for the number seeds.
! Note the way the function is called, for anyone not familiar. The value
! returned by the RANMAR function is fully scalable: I use it between 0 & 100
! by saying
!           x=100*RANMAR()
! Now, RMARIN and RANMAR share initialization variables in a common block,
! so they can be separated between subroutines, etc (again, for anyone
! unfamiliar, I call RMARIN as a subroutine from my main, and then include
! RANMAR as a function in some of my subroutines where I need to generate
! random numbers).
! **************************************************************************
! **************************************************************************
!      PROGRAM tstran
!      INTEGER ij, kl, i
!C These are the seeds needed to produce the test case results
!      ij = 1802
!      kl = 9373
!C Do the initialization
!      CALL rmarin(ij,kl)
!C Generate 20000 random numbers
!      DO 10 i = 1, 20000
!        x = ranmar()
!   10 CONTINUE
!      PRINT *,'GOT THIS FAR!'
!c If the random number generator is working properly, the next six random
!C numbers should be:
!C           6533892.0  14220222.0  7275067.0
!C           6172232.0  8354498.0   10633180.0
!      WRITE(*,*)'Correct:     6533892.0  14220222.0  7275067.0 '
!      WRITE(*,*)'Correct:     6172232.0  8354498.0   10633180.0'
!      WRITE(*,20) (4096.0*4096.0*ranmar(), i=1,6)
!   20 FORMAT (3f12.1)
!      END
! **************************************************************************
! **************************************************************************
      SUBROUTINE rmarin(ij,kl)
! This is the initialization routine for the random number generator
!     RANMAR()
! NOTE: The seed variables can have values between:    0 <= IJ <= 31328
!                                                      0 <= KL <= 30081
      REAL u(97), c, cd, cm
      INTEGER i97, j97
      LOGICAL test
      COMMON /raset1/ u, c, cd, cm, i97, j97, test
      DATA test /.false./
      IF( ij .LT. 0  .OR.  ij .GT. 31328  .OR.&
         kl .LT. 0  .OR.  kl .GT. 30081 ) THEN
        PRINT '(A)', ' The first random number seed must have a value &
             &between 0 and 31328'
        PRINT '(A)',' The second seed must have a value between 0 and &
             &30081'
        STOP
      ENDIF
      i = mod(ij/177, 177) + 2
      j = mod(ij    , 177) + 2
      k = mod(kl/169, 178) + 1
      l = mod(kl,     169)
      DO 2 ii = 1, 97
        s = 0.0
        t = 0.5
        DO 3 jj = 1, 24
          m = mod(mod(i*j, 179)*k, 179)
          i = j
          j = k
          k = m
          l = mod(53*l+1, 169)
          IF (mod(l*m, 64) .ge. 32) then
            s = s + t
          ENDIF
          t = 0.5 * t
    3   CONTINUE
        u(ii) = s
    2 CONTINUE
      c = 362436.0 / 16777216.0
      cd = 7654321.0 / 16777216.0
      cm = 16777213.0 /16777216.0
      i97 = 97
      j97 = 33
      test = .true.
      RETURN
      END
! **************************************************************************
! **************************************************************************
      real FUNCTION ranmar()
      REAL u(97), c, cd, cm
      INTEGER i97, j97
      LOGICAL test
      COMMON /raset1/ u, c, cd, cm, i97, j97, test
      IF( .NOT. test ) THEN
        PRINT '(A)',' Call the init routine (RMARIN) before calling RANMAR'
        STOP
      ENDIF
      uni = u(i97) - u(j97)
      IF( uni .LT. 0.0 ) uni = uni + 1.0
      u(i97) = uni
      i97 = i97 - 1
      IF(i97 .EQ. 0) i97 = 97
      j97 = j97 - 1
      IF(j97 .EQ. 0) j97 = 97
      c = c - cd
      IF( c .LT. 0.0 ) c = c + cm
      uni = uni - c
      IF( uni .LT. 0.0 ) uni = uni + 1.0
      ranmar = uni
      RETURN
      END
! **************************************************************************

