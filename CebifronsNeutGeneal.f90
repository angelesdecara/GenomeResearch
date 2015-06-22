!
! program to manage with all chromosomes
! with fped from genealogies
program ManageCebifrons
  implicit none
   integer :: i,j,k,l, icrom, igen, irep, n, iguy,igamete,gametetemp,idummy,imarker, io, m
   integer, parameter :: dp=KIND(1.d0), nindbase=7, nbinsdens=10000, ncrom=18
   integer, parameter :: nfem0=4,nmale0=1, nind0=nfem0+nmale0,nmarkers=104035
   integer, parameter :: nind=10, ngen=10, nrep=1000,ibdsize=2000000
   integer :: genot(nind,nmarkers,2),genottemp(nind,nmarkers,2),NumberMarkers(ncrom),nbins
   integer::  pos(nmarkers), ncross, geneal(nind,2),nfem,nmale,crosses(100),chr(nmarkers)
   integer :: sol(nind,2),dummy, genotbase(nindbase,nmarkers,2),genotBaseHere(nind0,nmarkers,2)
   real(dp):: fpedbase(nind,nind),fped(nind,nind), recombrate(ncrom), poisson, meandummy, sumx,x
   real(dp):: recomb(ncrom),fpedtemp(nind,nind)!=.362726d0 ! cM per Mb
   real (dp):: meanOH(0:ngen),varOH(0:ngen),meanGD(0:ngen),varGD(0:ngen)
   real(dp):: meanmalecontrib(0:ngen),meanfemalecontrib(0:ngen),varfemalecontrib(0:ngen)
   real(dp):: varmalecontrib(0:ngen), meanfixedpos(0:ngen),varfixedpos(0:ngen)
   real(dp) :: meanmaxcontrib(0:ngen),meanmincontrib(0:ngen),meanfped(0:ngen)
   real(dp)::oh,gd,fixpos, frecs(nmarkers,2),gd2,diff, f
   integer :: aux(nind)
   character(len=6)::a
   character(len=2)::base
   character(len=50)::fileHaps,markername, indName,sampleName
   character (len=50):: sizeforfile,formGenCrom
   real:: ranmar
   call rmarin(5879,6958)
   genot=99; genottemp=99

! everything 0.d0
   meanOH=0.d0;varOH=0.d0; meanGD=0.d0; varGD=0.d0
   meanmalecontrib=0.d0;meanfemalecontrib=0.d0;
   varmalecontrib=0.d0;varfemalecontrib=0.d0;
   meanfixedpos =0.d0; varfixedpos=0.d0; 
   meanmaxcontrib=0.d0; meanmincontrib=0.d0; meanfped=0.d0

   write(sizeforfile,*)ibdsize/1000
   sizeforfile=trim(adjustl(sizeforfile))
   write(formGenCrom,*)(ngen+1)*ncrom
   formGenCrom='('//trim(adjustl(formGenCrom))//'f12.8)'

!
! read phased genotypes data from chromosome files
!
   imarker=0; NumberMarkers=0
   do icrom=1,ncrom
      io=0
      write(fileHaps,*)icrom
      fileHaps='CebifronsNeut/chr'//trim(adjustl(fileHaps))//'.phased.haps'
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
   if(sum(NumberMarkers(1:ncrom)).ne.nmarkers)write(6,*)'NumberMarkers',sum(NumberMarkers(1:ncrom)),nmarkers
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

!
! 7 guys, 2 first are the males from negros, then 5 from Panay: female, female,male, female, female
!
   genotBaseHere(1,:,:)=genotbase(5,:,:)
   genotBaseHere(2,:,:)=genotbase(3,:,:)
   genotBaseHere(3,:,:)=genotbase(4,:,:)
   genotBaseHere(4,:,:)=genotbase(6,:,:)
   genotBaseHere(5,:,:)=genotbase(7,:,:)


   do irep=1,nrep

      ! start from real base
!      genot=genotBaseHere
     fped=0
     do i=1,nind
        fped(i,i)=0.5
     enddo

! open coancestries
     open(33,file='coancestries.txt')
     do i=1,400
        read(33,*)j,k,l,m,f
!
! using their pedigree
!
! 125 -> 1, male
! 103 -> 2
! 12 -> 3
! 120 -> 4
! 132 -> 5
        if((l.eq.125))then
           if(m.eq.125)fped(1,1)=f
           if(m.eq.103)fped(1,2)=f
           if(m.eq.12)fped(1,3)=f
           if(m.eq.120)fped(1,4)=f
           if(m.eq.132)fped(1,5)=f
        endif
        if((l.eq.103))then
           if(m.eq.103)fped(2,2)=f
           if(m.eq.125)fped(2,1)=f
           if(m.eq.12)fped(2,3)=f
           if(m.eq.120)fped(2,4)=f
           if(m.eq.132)fped(2,5)=f
        endif
        if((l.eq.12))then
           if(m.eq.12)fped(3,3)=f
           if(m.eq.103)fped(3,2)=f
           if(m.eq.125)fped(3,1)=f
           if(m.eq.120)fped(3,4)=f
           if(m.eq.132)fped(3,5)=f
        endif
        if((l.eq.120))then
           if(m.eq.120)fped(4,4)=f
           if(m.eq.103)fped(4,2)=f
           if(m.eq.12)fped(4,3)=f
           if(m.eq.125)fped(4,1)=f
           if(m.eq.132)fped(4,5)=f
        endif
        if((l.eq.132))then
           if(m.eq.132)fped(5,5)=f
           if(m.eq.103)fped(5,2)=f
           if(m.eq.12)fped(5,3)=f
           if(m.eq.120)fped(5,4)=f
           if(m.eq.125)fped(5,1)=f
        endif
     enddo
     close(33)
     write(6,*)'meanfped bef expansion',sum(fped(1:5,1:5))/25.d0
     fpedbase=fped

!population expansion to nind 
      do i=1,nind

         fped=fpedbase
         geneal(i,1)=1
         geneal(i,2)=ranmar()*nfem0+1+nmale0

         do icrom=1,ncrom
            do k=1,2
!
! crossovers for this guy and this gamete
!
               ncross=poisson(recombrate(icrom))
               do l=1,ncross
                  crosses(l)=ranmar()*(sum(NumberMarkers(1:icrom))-sum(NumberMarkers(1:icrom-1)))+sum(NumberMarkers(1:icrom-1))+1
               enddo
               call sortincrease(crosses,ncross)
               crosses(ncross+1)=sum(NumberMarkers(1:icrom))
!
!
               gametetemp=ranmar()*2+1
               genottemp(i,1+sum(NumberMarkers(1:icrom-1)):sum(NumberMarkers(1:icrom)),k)&
                    =genotBaseHere(geneal(i,k),1+sum(NumberMarkers(1:icrom-1)):sum(NumberMarkers(1:icrom)),gametetemp)
               do n=1,ncross+1
                  gametetemp=3-gametetemp
                  genottemp(i,crosses(n)+1:crosses(n+1),k)=genotBaseHere(geneal(i,k),crosses(n)+1:crosses(n+1),gametetemp)
               enddo
            enddo

         enddo
      enddo
      genot=genottemp

      do i=1,nind

         fpedtemp(i,i)=.5d0*(1+fped(geneal(i,1),geneal(i,2)))
         do j=i+1,nind
            fpedtemp(i,j)=.25d0*(fped(geneal(i,1),geneal(j,1))+&
                 fped(geneal(i,1),geneal(j,2))+&
                 fped(geneal(i,2),geneal(j,1))+&
                 fped(geneal(i,2),geneal(j,2)))
            fpedtemp(j,i)=fpedtemp(i,j)
         enddo
      enddo
      fped=fpedtemp
!
! the first two individuals 1,2 are males and the rest females 3:10
!
   nmale=2;  nfem=8
!
! at t=0 meanOH and so on
   write(6,*)'averfped 0',sum(fped(1:nind,1:nind))/dble(nind**2)
   write(6,*)'0',count(genot(:,:,:).gt.1)/dble(2*nind*nmarkers),sum(NumberMarkers(1:ncrom)),nmarkers


! 
!   meanOH=0.d0;varOH=0.d0; meanGD=0.d0; varGD=0.d0
!   meanfixedpos =0.d0; varfixedpos=0.d0
         call heter(genot,nind,nmarkers,oh,gd,fixpos)
         meanOH(0)=meanOH(0)+oh
         varOH(0)=varOH(0)+oh**2
         meanGD(0)=meanGD(0)+gd
         varGD(0)=varGD(0)+gd**2
         meanfixedpos(0)=meanfixedpos(0)+fixpos
         varfixedpos(0)=varfixedpos(0)+fixpos**2

         call freqs(genot,nind,nmarkers,frecs)
         gd2=0.d0; gd=0.d0
         do l=1,nmarkers
            gd2=gd2+(frecs(l,1))**2+(frecs(l,2))**2
            gd=gd+(count(genot(:,l,:)==0)/dble(2*nind))**2+(count(genot(:,l,:)==1)/dble(2*nind))**2
            diff=abs((dble(count(genot(:,l,:)==0))/dble(2*nind))**2+(dble(count(genot(:,l,:)==1))/dble(2*nind))**2-&
                 ((frecs(l,1))**2+(frecs(l,2))**2))
         enddo
         gd2=1.d0-gd2/dble(nmarkers)
         gd=1.d0-gd/dble(nmarkers)
         write(6,*)irep,meanGD(0)/irep,gd,gd2
   
!
! management during ngen generations
!
      do igen=1,ngen
!
! we need the coancestries to later calculate how many offspring each guy leaves
!
         write(6,*)'averfped',igen,sum(fped(1:nind,1:nind))/dble(nind**2)
!
! we use fibd in the annealing
!
 !        write(6,*)'before annealing',sum(fibd(1:nind,1:nind))/dble(nind**2)
         call modannealing(fped,nfem,nmale,sol)
!
! sol(nind,2), sol(i,1) is the number of sons that individual i leaves
!         sol(i,2) is the number of daughters that individual i leaves
!
! variance in contribs each gen

!         write(6,*)'sum contribs',sum(sol(1+nmale:nind,1)),sum(sol(1+nmale:nind,2))

         meandummy=0.d0
         do idummy=1,nmale
            if(sol(idummy,1)+sol(idummy,2).ge.1)meandummy=meandummy+1.d0
         enddo
         meandummy=meandummy/dble(nmale)
         meanmalecontrib(igen)=meanmalecontrib(igen)+meandummy
         varmalecontrib(igen)=varmalecontrib(igen)+meandummy**2
!
         meandummy=0.d0
         do idummy=1+nmale,nind
            if(sol(idummy,1)+sol(idummy,2).ge.1)meandummy=meandummy+1.d0
!            meandummy=sol(idummy,1)+sol(idummy,2)
!            write(6,*)idummy,meandummy
         enddo
         meandummy=meandummy/dble(nfem)
         meanfemalecontrib(igen)=meanfemalecontrib(igen)+meandummy
         varfemalecontrib(igen)=varfemalecontrib(igen)+meandummy**2

  ! fathers' genealogies
         dummy=0
         do n=1,2
            do i=1,nmale
               if(sol(i,n).gt.0)then
                  do j=1,sol(i,n)
                     dummy=dummy+1
                     aux(dummy)=i
!!                     geneal(dummy,1)=i ! geneal(i,1) is father of i
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
         dummy=0
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
                     crosses(l)=ranmar()*(sum(NumberMarkers(1:icrom))-sum(NumberMarkers(1:icrom-1)))+sum(NumberMarkers(1:icrom-1))
                  enddo
                  call sortincrease(crosses,ncross)
                  crosses(ncross+1)=sum(NumberMarkers(1:icrom))
!
!
                  gametetemp=ranmar()*2+1

                  genottemp(i,1+sum(NumberMarkers(1:icrom-1)):sum(NumberMarkers(1:icrom)),k)&
                       =genot(geneal(i,k),1+sum(NumberMarkers(1:icrom-1)):sum(NumberMarkers(1:icrom)),gametetemp)

                  do n=1,ncross+1
                     gametetemp=3-gametetemp
                     genottemp(i,crosses(n):crosses(n+1),k)=genot(geneal(i,k),crosses(n):crosses(n+1),gametetemp)
                  enddo
               enddo

            enddo
         enddo
         genot=genottemp
!
!
         do i=1,nind
            fpedtemp(i,i)=.5d0*(1+fped(geneal(i,1),geneal(i,2)))
            do j=i+1,nind
               fpedtemp(i,j)=.25d0*(fped(geneal(i,1),geneal(j,1))+&
                    fped(geneal(i,1),geneal(j,2))+&
                    fped(geneal(i,2),geneal(j,1))+&
                    fped(geneal(i,2),geneal(j,2)))
               fpedtemp(j,i)=fpedtemp(i,j)
            enddo
         enddo
         fped=fpedtemp


         oh=0.d0;gd=0.d0;fixpos=0.d0
         call heter(genot,nind,nmarkers,oh,gd,fixpos)
         meanOH(igen)=meanOH(igen)+oh
         varOH(igen)=varOH(igen)+oh**2
         meanGD(igen)=meanGD(igen)+gd
         varGD(igen)=varGD(igen)+gd**2
         meanfixedpos(igen)=meanfixedpos(igen)+fixpos
         varfixedpos(igen)=varfixedpos(igen)+fixpos**2

         meanfped(igen)=meanfped(igen)+sum(fped(1:nind,1:nind))/dble(nind**2)
      enddo
!
!
!
      open(54,file='CebifronsNeut/manGen.dat')
      write(54,'(a130,i4)')'#meanOH,varOH,meanGD,varGD,meanfixedpos,varfixedpos,meanfemalecontrib,varfemalecontrib,meanmalecontrib,varmalecontrib,meanfped',irep
      do igen=0,ngen
         write(54,'(i4,11f18.14)')igen,meanOH(igen)/dble(irep),varOH(igen)/dble(irep)-(meanOH(igen)/dble(irep))**2,&
              meanGD(igen)/dble(irep),varGD(igen)/dble(irep)-(meanGD(igen)/dble(irep))**2,&
              meanfixedpos(igen)/dble(irep),varfixedpos(igen)/dble(irep)-(meanfixedpos(igen)/dble(irep))**2,&
              meanfemalecontrib(igen)/dble(irep),varfemalecontrib(igen)/dble(irep)-(meanfemalecontrib(igen)/dble(irep))**2,&
              meanmalecontrib(igen)/dble(irep),varmalecontrib(igen)/dble(irep)-(meanmalecontrib(igen)/dble(irep))**2,&
              meanfped(igen)/dble(irep)
      enddo
      close(54)


   enddo

 end program ManageCebifrons
!
! sorting in increasing order the recomb points
! poisson
! ranmar
! fibd for nind,nind with ibdsize
! annealing for different number of males and females
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
!--- freqs alels
subroutine freqs(genot,nind,nmarkers,frec)
  implicit none
  integer, parameter :: dp=KIND(1.d0)
  integer :: i,j,k,l, nind,nmarkers,genot(nind,nmarkers,2),al
  real(dp):: frec(nmarkers,2), sumfrecs

  frec=0.d0
  do l=1,nmarkers
     do i=1,nind
        do k=1,2
           al=genot(i,l,k)+1
          ! if(al.gt.2)then
         !     write(6,*)'error',i,l,k,genot(i,l,k)
         !     pause
         !  endif
           frec(l,al)=frec(l,al)+1.d0
        enddo
     enddo
  enddo
!
  do l=1,nmarkers
     sumfrecs=0.d0
     do al=1,2
        sumfrecs=sumfrecs+frec(l,al)
     enddo
     do al=1,2
        frec(l,al)=frec(l,al)/sumfrecs
     enddo

!     if(frec(l,1)+frec(l,2).ne.1.d0)then
 !       write(6,*)'freq prob',frec(l,1)+frec(l,2)
 !    endif
  enddo

end subroutine freqs
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
! need to change to include pos
!
subroutine ibdcoanc(nind,ncrom,NumberMarkers,nmarkers,genot,pos,windowsize,fibd)
  implicit none
  integer, parameter :: dp=KIND(1.d0)
  integer :: i,j,k,l,nind,nmarkers,i1,i2,j1,j2, length0,windowsize
  integer :: genot(nind,nmarkers,2),lengthroh,ibdsize, pos(nmarkers)
  integer :: ncrom, NumberMarkers(ncrom), icrom
  real (dp) :: fibd(nind,nind),fibdpercrom

  fibd = 0.d0; fibdpercrom=0.d0
  do i=1,nind
     do j=i,nind
        do icrom=1,ncrom
           fibdpercrom=0.d0
           do i1=1,2
              do j1=1,2
                 if((j.gt.i).or.((i.eq.j).and.(j1.gt.i1)))then
                    l=1+sum(NumberMarkers(1:icrom-1));lengthroh=0;
                    if(icrom.eq.1)then
                       length0=pos(1)
                       l=1
                    else
                       length0=pos(1+sum(NumberMarkers(1:icrom-1)))
                    endif
                    do while (l.le.sum(NumberMarkers(1:icrom)))
                       if (genot(i,l,i1).eq.genot(j,l,j1))then
                          lengthroh = pos(l)-length0
                          l = l+1
                       else
                          if(lengthroh.ge.windowsize)then
                             fibdpercrom=fibdpercrom+lengthroh
                          endif
                          lengthroh=0
                          length0=pos(l)
                          l=l+1
                       endif
                    enddo
! now let's include the last ROH if homozygous at nloci
                    if (genot(i,sum(NumberMarkers(1:icrom)),i1).eq.genot(j,sum(NumberMarkers(1:icrom)),j1))then
                       if(lengthroh.ge.windowsize)then
                          fibdpercrom=fibdpercrom+lengthroh
                       endif
                    endif
                 endif
!
              enddo
           enddo
!
           fibdpercrom=fibdpercrom/dble(4*(pos(sum(NumberMarkers(1:icrom)))-pos(1+sum(NumberMarkers(1:icrom-1)))))


           fibd(i,j)=fibd(i,j)+fibdpercrom
        enddo
        fibd(i,j)=fibd(i,j)/dble(ncrom)
        fibd(j,i)=fibd(i,j)


        if (fibd(i,j).gt.1.d0)then
           write(6,*)'error fibd',fibd(i,j),i,j
           stop
        endif

     enddo
  enddo

!  write(6,*)'aver fibd in subr',sum(fibd(1:nind,1:nind))/dble(nind*nind)

end subroutine ibdcoanc
!====================
!-----------------------------------
!====================
subroutine ibddistrib(nind,ncrom,NumberMarkers,nmarkers,genot,pos,nbins,densibd)
  implicit none
  integer, parameter :: dp=KIND(1.d0)
  integer :: i,j,k,l,nind,nmarkers,i1,i2,j1,j2, length0,windowsize
  integer :: genot(nind,nmarkers,2),lengthroh,ibdsize, pos(nmarkers)
  integer :: ncrom, NumberMarkers(ncrom), icrom, nbins
  real (dp) :: densibd(ncrom,nbins+1),interv(ncrom), sumdens
!
! first interval per chromosome
!
  do icrom=1,ncrom
     interv(icrom)=(pos(sum(NumberMarkers(1:icrom)))-pos(sum(1+NumberMarkers(1:icrom-1))))/dble(nbins)
  enddo

  densibd=0.d0
  do i=1,nind
     do j=i,nind
        do icrom=1,ncrom

           do i1=1,2
              do j1=1,2
                 if((j.gt.i).or.((i.eq.j).and.(j1.gt.i1)))then
!
                    l=1+sum(NumberMarkers(1:icrom-1));lengthroh=0;
                    if(icrom.eq.1)then
                       length0=pos(1)
                       l=1
                    else
                       length0=pos(1+sum(NumberMarkers(1:icrom-1)))
                    endif
                    do while (l.le.sum(NumberMarkers(1:icrom)))
                       if (genot(i,l,i1).eq.genot(j,l,j1))then
                          lengthroh = pos(l)-length0
                          l = l+1
                       else
                          if(lengthroh.ge.2)then
                             densibd(icrom,int(lengthroh/interv(icrom))+1)=densibd(icrom,int(lengthroh/interv(icrom))+1)+1.d0
                          endif
                          lengthroh=0
                          length0=pos(l)
                          l=l+1
                       endif
                    enddo
! now let's include the last ROH if homozygous at nloci
                    if (genot(i,sum(NumberMarkers(1:icrom)),i1).eq.genot(j,sum(NumberMarkers(1:icrom)),j1))then
                       if(lengthroh.ge.2)then
                          densibd(icrom,int(lengthroh/interv(icrom))+1)=densibd(icrom,int(lengthroh/interv(icrom))+1)+1.d0
                       endif
                    endif
                    !
                 endif
              enddo
           enddo
        enddo
!
     enddo
  enddo

!
!normalise
!
  do icrom=1,ncrom
     sumdens=0.d0
     do i=1,nbins+1
        sumdens=sumdens+densibd(icrom,i)
     enddo
     densibd(icrom,:)=densibd(icrom,:)/sumdens
  enddo

end subroutine ibddistrib
!-------------------------------------------
! annealing!
!
subroutine modannealing(fibd,nfem,nmale,sol)
  implicit none
  integer, parameter :: dp=KIND(1.d0)
  integer :: nfem, nmale, m, n,nind
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
  solopt=sol
  eneopt=ene
  ivig=5000
! Levels loop
!  
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
! change?
        ch=0
        delta=dmax1(eneal-ene,0.)
        omega=exp(-delta/kt)
        if(omega.ge.1.d0) then
           ch=1
           if(eneal.lt.eneopt) then
              eneopt=eneal
              solopt=sola
           endif
        else
           if(dble(ranmar()).lt.omega) then
              ch=1
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
  write(6,*)'all loops done'

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

