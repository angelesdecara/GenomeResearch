!
! program to manage Cebifrons with segment based coancestry
! and including deleterious variants
! but not in the calculation of coancestry
! fibd
program ManageCebifrons
  implicit none
   integer :: i,j,k,l, icrom, igen, irep, n, iguy,igamete,gametetemp,idummy,imarker, io,jmarker
   integer, parameter :: dp=KIND(1.d0), nindbase=7, nbinsdens=10000, ncrom=18, ibdsize=100000
   integer, parameter :: nfem0=4,nmale0=1, nind0=nfem0+nmale0,nmarkers=107129
   integer, parameter :: nind=10, ngen=10, nrep=1000
   real (dp), parameter :: lambda=0.5d0, beta=1.d0,smean=0.005d0, hmean=.35d0
   integer, parameter :: ndelmarkers=3129
   integer :: genot(nind,nmarkers,2),genottemp(nind,nmarkers,2),NumberMarkers(ncrom),nbins
   integer :: genotdummy(nind,nmarkers,2)
   integer::  pos(nmarkers), ncross, geneal(nind,2),nfem,nmale,crosses(100),chr(nmarkers)
   integer :: sol(nind,2),dummy, genotbase(nindbase,nmarkers,2),genotBaseHere(nind0,nmarkers,2)
   real(dp):: fibd(nind,nind), recombrate(ncrom), poisson, meandummy, sumx,x
   real(dp):: recomb(ncrom)!=.362726d0 ! cM per Mb
   integer :: posdel(ndelmarkers), NumberDelMarkers(ncrom), chrdel(ndelmarkers)
   real (dp):: selec(ndelmarkers), heter(ndelmarkers), kheter
   real (dp):: meanOH(0:ngen),varOH(0:ngen),meanGD(0:ngen),varGD(0:ngen)
   real(dp):: meanmalecontrib(0:ngen),meanfemalecontrib(0:ngen),varfemalecontrib(0:ngen)
   real(dp):: varmalecontrib(0:ngen), meanfixedpos(0:ngen),varfixedpos(0:ngen)
   real(dp) :: meanmaxcontrib(0:ngen),meanmincontrib(0:ngen),meanfibd(0:ngen)
   real(dp):: densibd(ncrom,nbinsdens+1), meandensibd(0:ngen,ncrom,nbinsdens+1)
   real(dp)::oh,gd,fixpos, frecs(nmarkers,2),gd2,diff, meanfitness(0:ngen), fitness(nind)
   real(dp)::meandelfreq(0:ngen),meanhomdel(0:ngen),meanhetdel(0:ngen)
   integer :: delposinall(ndelmarkers), cont, aux(nind), locusstart
   logical:: selsite(nmarkers), first
   character(len=6)::a
   character(len=2)::base,baseref(nmarkers),basegood(ndelmarkers)
   character(len=50)::fileHaps,markername, indName,sampleName
   character (len=50):: sizeforfile,formGenCrom,filename
   character(len=10)::schar,hchar
   real:: ranmar
   real (dp):: shape, scale

   INTERFACE
      FUNCTION random_gamma(s, b, first) RESULT(fn_val)
        IMPLICIT NONE
        INTEGER, PARAMETER  :: dp = KIND(1.d0) !SELECTED_REAL_KIND(12, 60)
        REAL (dp), INTENT(IN)  :: s, b
        LOGICAL, INTENT(IN)    :: first
        REAL (dp)              :: fn_val
      END FUNCTION random_gamma
   end INTERFACE

   call rmarin(5879,6958)
   genot=99; genottemp=99

   write(sizeforfile,*)ibdsize/1000
   sizeforfile=trim(adjustl(sizeforfile))
   write(formGenCrom,*)(ngen+1)*ncrom
   formGenCrom='('//trim(adjustl(formGenCrom))//'f12.8)'


   write(schar,'(g0.2)')smean
   write(hchar,'(g0.2)')hmean
   filename='s'//trim(adjustl(schar))//'_h'//trim(adjustl(hchar))
! everything 0.d0
   meanOH=0.d0;varOH=0.d0; meanGD=0.d0; varGD=0.d0
   meanmalecontrib=0.d0;meanfemalecontrib=0.d0;
   varmalecontrib=0.d0;varfemalecontrib=0.d0;
   meanfixedpos =0.d0; varfixedpos=0.d0; 
   meanmaxcontrib=0.d0; meanmincontrib=0.d0; meanfibd=0.d0; meanfitness=0.d0; fitness=0.d0;
   meandelfreq=0.d0;meanhomdel=0.d0; meanhetdel=0.d0
!
! read phased genotypes data from chromosome files
!
   imarker=0; NumberMarkers=0
   do icrom=1,ncrom
      io=0
      write(fileHaps,*)icrom
      fileHaps='CebifronsNonNeut/chr'//trim(adjustl(fileHaps))//'.phased_D.haps'
      open(53,file=trim(adjustl(fileHaps)))
      do while((io.eq.0))
         imarker=imarker+1
!         write(6,*)io
         read(53,*,IOSTAT=io)idummy,markername,pos(imarker),baseref(imarker),base,((genotbase(iguy,imarker,igamete),igamete=1,2),iguy=1,nindbase)
!         write(6,*)idummy,markername,pos(imarker),baseref(imarker),base,((genotbase(iguy,imarker,igamete),igamete=1,2),iguy=1,nindbase)
!         pause
         chr(imarker)=icrom
      enddo
      close(53)
      imarker=imarker-1
      if(icrom.gt.1)then
         NumberMarkers(icrom)=imarker-sum(NumberMarkers(1:icrom-1))
      else if(icrom.eq.1)then
         NumberMarkers(icrom)=imarker
      endif
!      write(6,*)'icrom,nmarkers',icrom,imarker,NumberMarkers(icrom),nmarkers
   enddo
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
! deleterious variants
!
   selsite=.false.
   open(23,file='CebifronsNonNeut/delsites')
   do imarker=1,ndelmarkers
      read(23,*)icrom,posdel(imarker),basegood(imarker),base !! 1 is the deleterious
 !!     write(6,*)icrom,posdel(imarker),basegood(imarker),base 
!      pause
      chrdel(imarker)=icrom
      if(icrom.gt.1)then
         NumberDelMarkers(icrom)=imarker-sum(NumberDelMarkers(1:icrom-1))
      else if(icrom.eq.1)then
         NumberDelMarkers(icrom)=imarker
      endif
!
! find which marker this is
!
      do jmarker=1+sum(NumberMarkers(1:icrom-1)),sum(NumberMarkers(1:icrom))
         if(posdel(imarker).eq.pos(jmarker))then
            delposinall(imarker)=jmarker
            selsite(jmarker)=.true.
         endif
      enddo

   enddo
   close(23)

!
! freq of 1 del
!
   first=.true.
   shape = beta
   scale = smean/shape
  kheter = (beta/smean) * ( (2*hmean)**(-1./beta) -1.)

   do imarker=1,ndelmarkers
      selec(imarker)=random_gamma(shape,scale,first)
      heter(imarker)=ranmar()*exp(-kheter*selec(imarker))
      first=.false.

      write(29,*)selec(imarker),heter(imarker)
      if(baseref(delposinall(imarker)).ne.basegood(imarker))then
         do i=1,nindbase
            do k=1,2
               genotbase(i,delposinall(imarker),k)=1- genotbase(i,delposinall(imarker),k)
            enddo
         enddo
      endif
!
! initial freqs 
!
   enddo
   write(6,*)sum(selec)/dble(ndelmarkers),sum(heter)/dble(ndelmarkers)

!
! 7 guys, 2 first are the males from negros, then 5 from Panay: female, female,male, female, female
!
   genotBaseHere(1,:,:)=genotbase(5,:,:)
   genotBaseHere(2,:,:)=genotbase(3,:,:)
   genotBaseHere(3,:,:)=genotbase(4,:,:)
   genotBaseHere(4,:,:)=genotbase(6,:,:)
   genotBaseHere(5,:,:)=genotbase(7,:,:)


!
! before expansion, check allelic freqs
!
   write(6,*)'divers',dble(count(genotBaseHere==1))/dble(2*nmarkers*5)
   do i=1,ndelmarkers
      meandelfreq(0)=meandelfreq(0)+count(genotBaseHere(:,delposinall(i),:)==1)
   enddo
   write(6,*)'del divers',meandelfreq(0)/(2*ndelmarkers*5)

   meandelfreq(0)=0.d0

   do irep=1,nrep

      ! start from real base
!      genot=genotBaseHere


!population expansion to nind 
      do i=1,nind
         geneal(i,1)=1
         geneal(i,2)=ranmar()*nfem0+1+nmale0

         do icrom=1,ncrom
            do k=1,2
!
! crossovers for this guy and this gamete
!
               ncross=poisson(recombrate(icrom))
               locusstart=1+sum(NumberMarkers(1:icrom-1))
               do l=1,ncross
                  crosses(l)=ranmar()*(sum(NumberMarkers(1:icrom))-sum(NumberMarkers(1:icrom-1)))+sum(NumberMarkers(1:icrom-1))+1
               !   crosses(l)=ranmar()*nmarkers+1!(pos(nmarkers)-pos(1))+1
               enddo
               call sortincrease(crosses,ncross)
               crosses(ncross+1)=sum(NumberMarkers(1:icrom))
!               crosses(ncross+1)=nmarkers
!
!
               gametetemp=ranmar()*2+1
               genottemp(i,1+sum(NumberMarkers(1:icrom-1)):sum(NumberMarkers(1:icrom)),k)&
                    =genotBaseHere(geneal(i,k),1+sum(NumberMarkers(1:icrom-1)):sum(NumberMarkers(1:icrom)),gametetemp)
               do n=1,ncross
                  gametetemp=3-gametetemp
                  locusstart=crosses(n)+1
                  genottemp(i,crosses(n)+1:crosses(n+1),k)=genotBaseHere(geneal(i,k),crosses(n)+1:crosses(n+1),gametetemp)
               enddo
!               gametetemp=3-gametetemp
               genottemp(i,crosses(ncross)+1:crosses(ncross+1),k)=genotBaseHere(geneal(i,k),crosses(ncross)+1:crosses(ncross+1),gametetemp)

            enddo
         enddo
      enddo
      genot=genottemp
!
! the first two individuals 1,2 are males and the rest females 3:10
!
      nmale=2;  nfem=8
!
! at t=0 meanOH and so on
      do i=1,nind
         fitness(i)=1.d0
         cont=0
         do imarker=1,ndelmarkers
            fitness(i)=fitness(i)*(1.d0-&
                 heter(imarker)*selec(imarker)*dble( (1-genot(i,delposinall(imarker),1))*genot(i,delposinall(imarker),2)+(1-genot(i,delposinall(imarker),2))*genot(i,delposinall(imarker),1) )-&
                 selec(imarker)* dble( (genot(i,delposinall(imarker),1))*(genot(i,delposinall(imarker),2)) ) )
            cont=cont+genot(i,delposinall(imarker),1)+genot(i,delposinall(imarker),2)
!
            meanhomdel(0)=meanhomdel(0)+genot(i,delposinall(imarker),1)*genot(i,delposinall(imarker),2)/dble(nind*ndelmarkers)
            meanhetdel(0)=meanhetdel(0)+(genot(i,delposinall(imarker),1)*(1-genot(i,delposinall(imarker),2))+(1-genot(i,delposinall(imarker),1))*genot(i,delposinall(imarker),2))/dble(nind*ndelmarkers)
         enddo
      enddo
      do imarker=1,ndelmarkers
         meandelfreq(0)=meandelfreq(0)+count(genot(:,delposinall(imarker),:)==1)/dble(2*nind*ndelmarkers)
      enddo

         fibd=0.d0
         call ibdcoancNoDel(nind,ncrom,NumberMarkers,nmarkers,selsite,genot,pos,ibdsize,fibd)
!         write(6,*)'averfibd 0',sum(fibd(1:nind,1:nind))/dble(nind**2)
!         write(6,*)'freq',meandelfreq(0)/dble(irep),meanhomdel(0)/dble(irep),meanhetdel(0)/dble(irep)

         meanfibd(0)=meanfibd(0)+sum(fibd(1:nind,1:nind))/dble(nind**2)

! 
!   meanOH=0.d0;varOH=0.d0; meanGD=0.d0; varGD=0.d0
!   meanfixedpos =0.d0; varfixedpos=0.d0
         call heteroz(genot,nind,nmarkers,oh,gd,fixpos)
         meanOH(0)=meanOH(0)+oh
         varOH(0)=varOH(0)+oh**2
         meanGD(0)=meanGD(0)+gd
         varGD(0)=varGD(0)+gd**2
         meanfixedpos(0)=meanfixedpos(0)+fixpos
         varfixedpos(0)=varfixedpos(0)+fixpos**2

         meanfitness(0)=meanfitness(0)+sum(fitness)/dble(nind)

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
!         write(6,*)irep,meanGD(0)/irep,gd,gd2
   
!
! management during ngen generations
!
      do igen=1,ngen
!
! we need the coancestries to later calculate how many offspring each guy leaves
!
         fibd=0.d0
         call ibdcoancNoDel(nind,ncrom,NumberMarkers,nmarkers,selsite,genot,pos,ibdsize,fibd)
!         call ibdcoanc(nind,ncrom,NumberMarkers,nmarkers,genot,pos,ibdsize,fibd)
!         write(6,*)'averfibd',igen,sum(fibd(1:nind,1:nind))/dble(nind**2)
!
! we use fibd in the annealing
!
 !        write(6,*)'before annealing',sum(fibd(1:nind,1:nind))/dble(nind**2)
         call modannealing(fibd,nfem,nmale,sol)
!
! sol(nind,2), sol(i,1) is the number of sons that individual i leaves
!         sol(i,2) is the number of daughters that individual i leaves
!
! variance in contribs each gen

!         write(6,*)'sum contribs',sum(sol(1+nmale:nind,1)),sum(sol(1+nmale:nind,2))

         meandummy=0.d0
         do idummy=1,nmale
            if(sol(idummy,1)+sol(idummy,2).ge.1) meandummy=meandummy+1
         enddo
         meandummy=meandummy/dble(nmale)
         meanmalecontrib(igen)=meanmalecontrib(igen)+meandummy
         varmalecontrib(igen)=varmalecontrib(igen)+meandummy**2
!
         meandummy=0.d0
         do idummy=1+nmale,nind
            if(sol(idummy,1)+sol(idummy,2).ge.1) meandummy=meandummy+1
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
!!                     geneal(dummy,2)=i ! geneal(i,2) is mother of i
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
         i=1
         do idummy=1,10000*nind
            if(i.gt.nind)exit
            do icrom=1,ncrom
               do k=1,2
!
! crossovers for this guy and this gamete
!
                  ncross=poisson(recombrate(icrom))
                  do l=1,ncross
                     crosses(l)=ranmar()*(sum(NumberMarkers(1:icrom))-sum(NumberMarkers(1:icrom-1)))+sum(NumberMarkers(1:icrom-1))+1
                  !crosses(l)=ranmar()*(pos(nmarkers)-pos(1))+1
                  enddo
                  call sortincrease(crosses,ncross)
                  crosses(0)=1+sum(NumberMarkers(1:icrom-1))
                  crosses(ncross+1)=sum(NumberMarkers(1:icrom))
!               crosses(ncross+1)=pos(nmarkers)
!
!
                  locusstart=1+sum(NumberMarkers(1:icrom-1))
                  gametetemp=int(ranmar()*2)+1
                  genottemp(i,1+sum(NumberMarkers(1:icrom-1)):sum(NumberMarkers(1:icrom)),k)&
                       =genot(geneal(i,k),1+sum(NumberMarkers(1:icrom-1)):sum(NumberMarkers(1:icrom)),gametetemp)
                  do n=1,ncross
                     gametetemp=3-gametetemp
                     locusstart=crosses(n)+1
                     genottemp(i,crosses(n)+1:crosses(n+1),k)=genot(geneal(i,k),crosses(n)+1:crosses(n+1),gametetemp)
                  enddo
                  genottemp(i,crosses(ncross)+1:crosses(ncross+1),k)=genot(geneal(i,k),crosses(ncross)+1:crosses(ncross+1),gametetemp)
               enddo

            enddo
!
! evaluate fitness to decide whether to keep it
!
            fitness(i)=1.d0
            do imarker=1,ndelmarkers

               fitness(i)=fitness(i)*(1.d0-&
                    heter(imarker)*selec(imarker)*dble( (1-genot(i,delposinall(imarker),1))*genot(i,delposinall(imarker),2)+(1-genot(i,delposinall(imarker),2))*genot(i,delposinall(imarker),1) )-&
                    selec(imarker)* dble( (genot(i,delposinall(imarker),1))*(genot(i,delposinall(imarker),2)) ) )

            enddo
            if(dble(ranmar()).le.fitness(i))then !! keep
!               write(6,*)'here',idummy,i
               genotdummy(i,:,:)=genottemp(i,:,:)
               i=i+1
            endif


         enddo
         genot=genotdummy
         write(6,*)'attempts segman',idummy,i
!
         oh=0.d0;gd=0.d0;fixpos=0.d0
         call heteroz(genot,nind,nmarkers,oh,gd,fixpos)
         meanOH(igen)=meanOH(igen)+oh
         varOH(igen)=varOH(igen)+oh**2
         meanGD(igen)=meanGD(igen)+gd
         varGD(igen)=varGD(igen)+gd**2
         meanfixedpos(igen)=meanfixedpos(igen)+fixpos
         varfixedpos(igen)=varfixedpos(igen)+fixpos**2

         fibd=0.d0

         call ibdcoancNoDel(nind,ncrom,NumberMarkers,nmarkers,selsite,genot,pos,ibdsize,fibd)

!         call ibdcoanc(nind,ncrom,NumberMarkers,nmarkers,genot,pos,ibdsize,fibd)
         meanfibd(igen)=meanfibd(igen)+sum(fibd(1:nind,1:nind))/dble(nind**2)

         call ibddistrib(nind,ncrom,NumberMarkers,nmarkers,genot,pos,nbinsdens,densibd)
         do icrom=1,ncrom
            do i=1,nbinsdens+1
               meandensibd(igen,icrom,i)=meandensibd(igen,icrom,i)+densibd(icrom,i)
            enddo
         enddo

         do imarker=1,ndelmarkers
            do i=1,nind
               meanhomdel(igen)=meanhomdel(igen)+genot(i,delposinall(imarker),1)*genot(i,delposinall(imarker),2)/dble(nind*ndelmarkers)
               meanhetdel(igen)=meanhetdel(igen)+(genot(i,delposinall(imarker),1)*(1-genot(i,delposinall(imarker),2))+(1-genot(i,delposinall(imarker),1))*genot(i,delposinall(imarker),2))/dble(nind*ndelmarkers)
            enddo
            meandelfreq(igen)=meandelfreq(igen)+count(genot(:,delposinall(imarker),:)==1)/dble(2*nind*ndelmarkers)
         enddo
!
         do i=1,nind
            fitness(i)=1.d0
            do imarker=1,ndelmarkers
               fitness(i)=fitness(i)*(1.d0-&
                    heter(imarker)*selec(imarker)*dble( (1-genot(i,delposinall(imarker),1))*genot(i,delposinall(imarker),2)+(1-genot(i,delposinall(imarker),2))*genot(i,delposinall(imarker),1) )-&
                    selec(imarker)* dble( (genot(i,delposinall(imarker),1))*(genot(i,delposinall(imarker),2)) ) )

            enddo            
         enddo
         meanfitness(igen)=meanfitness(igen)+sum(fitness)/dble(nind)

!         write(6,*)' irep,igen',irep,igen,meanOH(igen)/irep,meanfitness(igen)/dble(irep)
!         write(6,*)'freq',meandelfreq(igen)/dble(irep),meanhomdel(igen)/dble(irep),meanhetdel(igen)/dble(irep)

      enddo
!
!
!
      open(54,file='CebifronsNonNeut/manSegsNoDel_muk'//trim(adjustl(filename))//'_'//trim(adjustl(sizeforfile))//'kb.dat')
      write(54,'(a175,i4)')'#meanOH,varOH,meanGD,varGD,meanfixedpos,varfixedpos,meanfemalecontrib,varfemalecontrib,meanmalecontrib,varmalecontrib,meanfibd,meanfitness,meandelfreq,meanhomdel,meanhetdel',irep
      do igen=0,ngen
         write(54,'(i4,15f18.14)')igen,meanOH(igen)/dble(irep),varOH(igen)/dble(irep)-(meanOH(igen)/dble(irep))**2,&
              meanGD(igen)/dble(irep),varGD(igen)/dble(irep)-(meanGD(igen)/dble(irep))**2,&
              meanfixedpos(igen)/dble(irep),varfixedpos(igen)/dble(irep)-(meanfixedpos(igen)/dble(irep))**2,&
              meanfemalecontrib(igen)/dble(irep),varfemalecontrib(igen)/dble(irep)-(meanfemalecontrib(igen)/dble(irep))**2,&
              meanmalecontrib(igen)/dble(irep),varmalecontrib(igen)/dble(irep)-(meanmalecontrib(igen)/dble(irep))**2,&
              meanfibd(igen)/dble(irep),meanfitness(igen)/dble(irep),&
              meandelfreq(igen)/dble(irep),meanhomdel(igen)/dble(irep),meanhetdel(igen)/dble(irep)
      enddo
      close(54)

      open(55,file='distIBDsegsCeb_muk005.dat')
!      write(6,*)'densibd ', meandensibd(1:1,1,1)
      do i=1,nbinsdens+1
         write(55,trim(adjustl(formGenCrom)))((meandensibd(igen,icrom,i)/dble(irep),icrom=1,ncrom),igen=0,ngen)
      enddo !!  meandensibd(0:ngen,ncrom,nbinsdens+1)
      close(55)


   enddo

 end program ManageCebifrons
!
! sorting in increasing order the recomb points
! poisson
! ranmar
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
!---------- oh, gd
 subroutine heteroz(genot,nind,nmarkers,oh,gd,fixpos)
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

 end subroutine heteroz
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
!====================
subroutine ibdcoancNoDel(nind,ncrom,NumberMarkers,nmarkers,selsite,genot,pos,windowsize,fibd)
!
!ibdcoanc(nind,ncrom,NumberMarkers,nmarkers,genot,pos,windowsize,fibd)
  implicit none
  integer, parameter :: dp=KIND(1.d0)
  integer :: i,j,k,l,nind,nmarkers,i1,i2,j1,j2, length0,windowsize
  integer :: genot(nind,nmarkers,2),lengthroh, pos(nmarkers)
  integer :: ncrom, NumberMarkers(ncrom), icrom, cont(nind,nind)
  real (dp) :: fibd(nind,nind),fibdpercrom
  logical :: selsite(nmarkers)

  fibd = 0.d0; fibdpercrom=0.d0
  do i=1,nind
     do j=i,nind
        do icrom=1,ncrom
           fibdpercrom=0.d0; cont=0
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
!
!in case the first one is a deleterious variant
!
                    do while (.not.selsite(l))
                       l=l+1
                       length0=pos(l)
                    enddo
                    if(.not.selsite(l))l=l-1

                    do while (l.le.sum(NumberMarkers(1:icrom)))
                       if(.not.selsite(l))then
                          if (genot(i,l,i1).eq.genot(j,l,j1))then
                             lengthroh = pos(l)-length0
                             l = l+1
                             cont(i,j)=cont(i,j)+1
                          else
                             if(lengthroh.ge.windowsize)then
                                fibdpercrom=fibdpercrom+lengthroh
                             endif
                             lengthroh=0
                             length0=pos(l)
                             l=l+1
                          endif
                       else
                          l=l+1
                       endif
                    enddo
! now let's include the last ROH if homozygous at nloci
                    if (genot(i,sum(NumberMarkers(1:icrom)),i1).eq.genot(j,sum(NumberMarkers(1:icrom)),j1))then
                       if(lengthroh.ge.windowsize)then
                          fibdpercrom=fibdpercrom+lengthroh
                       endif
                    endif
!
                 endif
!
              enddo
           enddo
!
           if(i.ne.j)then
              fibdpercrom=fibdpercrom/dble(4*(pos(sum(NumberMarkers(1:icrom)))-pos(1+sum(NumberMarkers(1:icrom-1)))))
           else
              fibdpercrom=fibdpercrom/dble((pos(sum(NumberMarkers(1:icrom)))-pos(1+sum(NumberMarkers(1:icrom-1)))))
           endif

!           if((i.eq.1).and.((j.eq.1)))then
!              write(6,*)'----',i,j,icrom,'---'
!              write(6,*)fibdpercrom,cont(i,j)/dble(NumberMarkers(icrom))
!          else if((i.eq.1).and.((j.eq.2)))then
!              write(6,*)'----',i,j,icrom,'---'
!              write(6,*)fibdpercrom,cont(i,j)/dble(NumberMarkers(icrom))
!           endif

           fibd(i,j)=fibd(i,j)+fibdpercrom
        enddo
        fibd(i,j)=fibd(i,j)/dble(ncrom)
        if(i.eq.j)fibd(i,j)=(1.d0+fibd(i,j))/2.d0
        fibd(j,i)=fibd(i,j)
        cont(j,i)=cont(i,j)
     enddo
  enddo


  !write(6,*)windowsize,fibd(1,1),fibd(1,2),pos(1)
 ! write(6,*)cont(1,1)/dble(nmarkers),cont(1,2)/dble(4*nmarkers)
!  write(6,*)(dble((pos(sum(NumberMarkers(1:icrom)))-pos(1+sum(NumberMarkers(1:icrom-1))))),icrom=1,2)
  !stop

end subroutine ibdcoancNoDel
!-----------------------------------------------------------------------
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
     interv(icrom)=(pos(sum(NumberMarkers(1:icrom)))-pos(1+sum(NumberMarkers(1:icrom-1))))/dble(nbins)
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
 !                            write(6,*)'in ibddist',lengthroh,int(lengthroh/interv(icrom))
                             densibd(icrom,int(lengthroh/interv(icrom))+1)=densibd(icrom,int(lengthroh/interv(icrom))+1)+1.d0
!                             stop
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
     if(sumdens.eq.0.d0)write(6,*)'no rohs?',icrom,sumdens,sum(NumberMarkers(1:icrom-1)),pos(sum(NumberMarkers(1:icrom-1))+1)
     do i=1,nbins+1
        if(densibd(icrom,i).gt.0.d0)densibd(icrom,i)=densibd(icrom,i)/sumdens
     enddo
  enddo

end subroutine ibddistrib
!-----------------------------------------------------------------------
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
!-------------------
! shape k = alpha = s
! scale theta = 1/beta = b
FUNCTION random_gamma(s, b, first) RESULT(fn_val)
!
! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

!     N.B. This version is in `double precision' and includes scaling

!     FUNCTION GENERATES A RANDOM GAMMA VARIATE.
!     CALLS EITHER random_gamma1 (S > 1.0)
!     OR random_exponential (S = 1.0)
!     OR random_gamma2 (S < 1.0).

!     S = SHAPE PARAMETER OF DISTRIBUTION (0 < REAL).
!     B = Scale parameter

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), INTENT(IN)  :: s, b
LOGICAL, INTENT(IN)    :: first
REAL (dp)              :: fn_val

! Local parameters
REAL (dp), PARAMETER  :: one = 1.0_dp, zero = 0.0_dp

IF (s <= zero) THEN
  WRITE(*, *) 'SHAPE PARAMETER VALUE MUST BE POSITIVE'
  STOP
END IF

IF (s > one) THEN
  fn_val = random_gamma1(s, first)
ELSE IF (s < one) THEN
  fn_val = random_gamma2(s, first)
ELSE IF (s == one) THEN
  fn_val = dble(random_exponential())
END IF

! Now scale the random variable
fn_val = b * fn_val
RETURN

CONTAINS

FUNCTION random_exponential() RESULT(ran_exp)
! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A NEGATIVE EXPONENTIAL DlSTRIBUTION WlTH DENSITY PROPORTIONAL
! TO EXP(-random_exponential), USING INVERSION.

IMPLICIT NONE
REAL  :: ran_exp
!     Local variable
REAL  :: r

DO
  CALL RANDOM_NUMBER(r)
  IF (r > zero) EXIT
END DO

ran_exp = -LOG(r)
RETURN

END FUNCTION random_exponential


FUNCTION random_gamma1(s, first) RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO GAMMA**(S-1)*EXP(-GAMMA),
! BASED UPON BEST'S T DISTRIBUTION METHOD

!     S = SHAPE PARAMETER OF DISTRIBUTION
!          (1.0 < REAL)

REAL (dp), INTENT(IN)  :: s
LOGICAL, INTENT(IN)    :: first
REAL (dp)              :: fn_val

!     Local variables
REAL (dp)             :: d, r, g, f, x
REAL (dp), SAVE       :: b, h
REAL (dp), PARAMETER  :: sixty4 = 64.0_dp, three = 3.0_dp, pt75 = 0.75_dp,  &
                         two = 2.0_dp, half = 0.5_dp

IF (s <= one) THEN
  WRITE(*, *) 'IMPERMISSIBLE SHAPE PARAMETER VALUE'
  STOP
END IF

IF (first) THEN                        ! Initialization, if necessary
  b = s - one
  h = SQRT(three*s - pt75)
END IF

DO
  CALL RANDOM_NUMBER(r)
  g = r - r*r
  IF (g <= zero) CYCLE
  f = (r - half)*h/SQRT(g)
  x = b + f
  IF (x <= zero) CYCLE
  CALL RANDOM_NUMBER(r)
  d = sixty4*g*(r*g)**2
  IF (d <= zero) EXIT
  IF (d*x < x - two*f*f) EXIT
  IF (LOG(d) < two*(b*LOG(x/b) - f)) EXIT
END DO
fn_val = x

RETURN
END FUNCTION random_gamma1



FUNCTION random_gamma2(s, first) RESULT(fn_val)

! Adapted from Fortran 77 code from the book:
!     Dagpunar, J. 'Principles of random variate generation'
!     Clarendon Press, Oxford, 1988.   ISBN 0-19-852202-9

! FUNCTION GENERATES A RANDOM VARIATE IN [0,INFINITY) FROM
! A GAMMA DISTRIBUTION WITH DENSITY PROPORTIONAL TO
! GAMMA2**(S-1) * EXP(-GAMMA2),
! USING A SWITCHING METHOD.

!    S = SHAPE PARAMETER OF DISTRIBUTION
!          (REAL < 1.0)

REAL (dp), INTENT(IN)  :: s
LOGICAL, INTENT(IN)    :: first
REAL (dp)              :: fn_val

!     Local variables
REAL (dp)             :: r, x, w
REAL (dp), SAVE       :: a, p, c, uf, vr, d
REAL (dp), PARAMETER  :: vsmall = EPSILON(one)

IF (s <= zero .OR. s >= one) THEN
  WRITE(*, *) 'SHAPE PARAMETER VALUE OUTSIDE PERMITTED RANGE'
  STOP
END IF

IF (first) THEN                        ! Initialization, if necessary
  a = one - s
  p = a/(a + s*EXP(-a))
  IF (s < vsmall) THEN
    WRITE(*, *) 'SHAPE PARAMETER VALUE TOO SMALL'
    STOP
  END IF
  c = one/s
  uf = p*(vsmall/a)**s
  vr = one - vsmall
  d = a*LOG(a)
END IF

DO
  CALL RANDOM_NUMBER(r)
  IF (r >= vr) THEN
    CYCLE
  ELSE IF (r > p) THEN
    x = a - LOG((one - r)/(one - p))
    w = a*LOG(x)-d
  ELSE IF (r > uf) THEN
    x = a*(r/p)**c
    w = x
  ELSE
    fn_val = zero
    RETURN
  END IF

  CALL RANDOM_NUMBER(r)
  IF (one-r <= w .AND. r > zero) THEN
    IF (r*(w + one) >= one) CYCLE
    IF (-LOG(r) <= w) CYCLE
  END IF
  EXIT
END DO

fn_val = x
RETURN

END FUNCTION random_gamma2

END FUNCTION random_gamma


RECURSIVE SUBROUTINE quick_sort(list, order)

! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.

IMPLICIT NONE
INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(12, 60)

REAL (dp), DIMENSION (:), INTENT(IN OUT)  :: list
INTEGER, DIMENSION (:), INTENT(OUT)       :: order

! Local variable
INTEGER :: i

DO i = 1, SIZE(list)
  order(i) = i
END DO

CALL quick_sort_1(1, SIZE(list))
RETURN


CONTAINS


RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
REAL (dp)           :: reference, temp
INTEGER, PARAMETER  :: max_simple_sort_size = 6

IF (right_end < left_end + max_simple_sort_size) THEN
  ! Use interchange sort for small lists
  CALL interchange_sort(left_end, right_end)

ELSE
  ! Use partition ("quick") sort
  reference = list((left_end + right_end)/2)
  i = left_end - 1
  j = right_end + 1

  DO
    ! Scan list from left end until element >= reference is found
    DO
      i = i + 1
      IF (list(i) >= reference) EXIT
    END DO
    ! Scan list from right end until element <= reference is found
    DO
      j = j - 1
      IF (list(j) <= reference) EXIT
    END DO


    IF (i < j) THEN
      ! Swap two out-of-order elements
      temp = list(i)
      list(i) = list(j)
      list(j) = temp
      itemp = order(i)
      order(i) = order(j)
      order(j) = itemp
    ELSE IF (i == j) THEN
      i = i + 1
      EXIT
    ELSE
      EXIT
    END IF
  END DO

  IF (left_end < j) CALL quick_sort_1(left_end, j)
  IF (i < right_end) CALL quick_sort_1(i, right_end)
END IF

RETURN
END SUBROUTINE quick_sort_1


SUBROUTINE interchange_sort(left_end, right_end)

INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
REAL (dp)           :: temp

DO i = left_end, right_end - 1
  DO j = i+1, right_end
    IF (list(i) > list(j)) THEN
      temp = list(i)
      list(i) = list(j)
      list(j) = temp
      itemp = order(i)
      order(i) = order(j)
      order(j) = itemp
    END IF
  END DO
END DO

RETURN
END SUBROUTINE interchange_sort

END SUBROUTINE quick_sort
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

