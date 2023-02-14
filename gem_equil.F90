MODULE gem_equil
  IMPLICIT NONE
  integer :: itube,iperi,iperidf,ibunit,icandy=1,isprime=1,ildu=1,eldu=1,isym=0,igmrkr=0,igmrkre=0
  real :: mimp=2,mcmp=12,chgi=1,chgc=6
  real :: elon0=1.0,tria0=0.0,rmaj0=1000.0,r0,a, selon0=0.0,&
             stria0=0.0,rmaj0p=-0.0,q0p=0.006,q0=1.4, elonp0=0.,triap0=0.,erp=0.01,er0=0.,q0abs
  real :: beta,Rovera,shat0,teti,tcti,rhoia,Rovlni,Rovlti,Rovlne,Rovlte,Rovlnc,Rovltc,ncne,nuacs
  real :: rhoi !yjhu added
  real :: gamma_E,mach
  real :: f0, f0p,candyf0p,bunit,debye
  real :: rin,rout,dr,dth,delz,jacmax
  real :: cn0e,cn0i,cn0b,cn0c,n0emax,n0imax,n0bmax,n0cmax
  real :: r0a,lxa,lymult,lxmult,delra,delri,delre,delrn,rina,routa,betai,betae, &
               tir0,xnir0

  integer :: nr=400,nr2=400,ntheta=1024,idiag=0,isgnf=1,isgnq=-1
  real,dimension(:,:),allocatable :: bfld,qhat,radius,gr,gth,grdgt,grcgt, &
                                        gxdgy,dydr,dbdr,dbdth,jacob, &
                                        yfn,hght,thflx
  real,dimension(:),allocatable :: psi,f,psip,sf,jacoba,jfn,zfnth,thfnz,&
                                      t0i,t0e,t0b,t0c,t0ip,t0ep,t0bp,t0cp,&
                                      xn0i,xn0e,xn0c,xn0b,xn0ip,xn0ep,xn0bp,&
                                      xn0cp,vpari,vparc,vparb,&
                                      vparip,vparcp,vparbp, &
                                      capti,capte,captb,captc,capni,capne,&
                                      capnb,capnc,zeff,nue0,phinc,phincp,&
                                      er,upari,dipdr

  !arrays for computing (b. grad rho) effect
  real,dimension(:,:), allocatable :: e1gx,e1gy,e2gx,e2gy,bdge1gx,bdge1gy,bdge2gx,bdge2gy

  real,dimension(:,:),allocatable :: t0s,xn0s,capts,capns,vpars,vparsp
  real,dimension(:),allocatable :: cn0s,n0smax,tgis
  real :: tge
  real:: bu,tu,nu,xu,frequ,vu,eru

!for including bstar effects
  real,dimension(:),allocatable :: psip2
  real,dimension(:,:),allocatable :: curvbz,srbr,srbz,thbr,thbz,prsrbr,prsrbz,pthsrbr,pthsrbz,bdcrvb

!  real,external :: erf

contains
  subroutine new_equil()
    use gem_com,only: myid
    implicit none
    real(8) :: pi,pi2,r,th,s,c1,c2,lti,s0,delsi,lte,delse,teti,lnh,delsh,sh
    real(8) :: kappat,kappan,wt,wn
    parameter(kappat=6.9,kappan=2.2,wt=0.3,wn=0.3)
    parameter(c1=0.43236,c2=2.33528,lti=1200.0,s0=0.5,delsi=0.9, &
         lte=1200.0,delse=0.9,lnh=600,delsh=0.2,sh=0.5)
    integer :: i,j,k,m,i1,j1,j2
    real(8) :: dum,x,denom,dum1,dum2,dum3,dum4,xp,xpp,cutoff,delrai=0.01,delrao=0.01,psitmax,psin
    real(8) :: vpar,vperpc,vperp,v
    real(8) :: e,proton,omegau
    real(8) :: rdum,rdum1,thdum,thdum1,wx0,wx1,wy0,wy1
    integer :: jleft, jright
    real :: drgem,dthgem
    integer :: nrgem=501,nthgem=1025,nrxgc=501
    real, dimension(:,:), allocatable :: dqhdr    
    real, dimension(:,:), allocatable :: rdata,zdata
    real, dimension(:), allocatable :: rgem,psigem,qgem,fgem,tigem,tegem,nigem,negem,psi_nor_output

    !xgc profile data for n T
    integer, dimension(:), allocatable :: itoxgc
    real, dimension(:), allocatable :: psinxgc,nexgc,texgc,tixgc
    
    !arays for computing (b.grad rho) effect in vorticity
    real(8) :: e1r(0:nr,0:ntheta),e1z(0:nr,0:ntheta),e2r(0:nr,0:ntheta),e2z(0:nr,0:ntheta),e2zet(0:nr,0:ntheta)
    real(8) :: bdge1r(0:nr,0:ntheta),bdge1z(0:nr,0:ntheta),bdge2r(0:nr,0:ntheta),bdge2z(0:nr,0:ntheta),bdge2zet(0:nr,0:ntheta)      

    !global equilibrium data
    allocate(bfld(0:nr,0:ntheta),qhat(0:nr,0:ntheta),radius(0:nr,0:ntheta), &
         gr(0:nr,0:ntheta),gth(0:nr,0:ntheta),grdgt(0:nr,0:ntheta), &
         grcgt(0:nr,0:ntheta),gxdgy(0:nr,0:ntheta),dydr(0:nr,0:ntheta),&
         dbdr(0:nr,0:ntheta),dbdth(0:nr,0:ntheta),dqhdr(0:nr,0:ntheta),&
         jacob(0:nr,0:ntheta),jfn(0:ntheta), zfnth(0:ntheta),thfnz(0:ntheta),&
         yfn(0:nr,0:ntheta),hght(0:nr,0:ntheta),thflx(0:nr,0:ntheta))

    allocate(sf(0:nr), psi(0:nr),&
         psip(0:nr),&
         f(0:nr),jacoba(0:nr),t0i(0:nr),t0c(0:nr),t0e(0:nr),t0ip(0:nr),&
         t0ep(0:nr),capti(0:nr),captc(0:nr),capte(0:nr),&
         xn0e(0:nr),xn0i(0:nr),xn0c(0:nr),capne(0:nr),capni(0:nr),capnc(0:nr),zeff(0:nr),nue0(0:nr),&
         vpari(0:nr),vparc(0:nr),vparb(0:nr),phinc(0:nr), &
         vparip(0:nr),vparcp(0:nr),vparbp(0:nr),phincp(0:nr),er(0:nr), &
         upari(0:nr),dipdr(0:nr))

    allocate(psip2(0:nr), curvbz(0:nr,0:ntheta),srbr(0:nr,0:ntheta),srbz(0:nr,0:ntheta),&
         thbr(0:nr,0:ntheta),thbz(0:nr,0:ntheta),bdcrvb(0:nr,0:ntheta), &
         prsrbr(0:nr,0:ntheta),prsrbz(0:nr,0:ntheta), &
         pthsrbr(0:nr,0:ntheta),pthsrbz(0:nr,0:ntheta))

    allocate(cn0s(1:5),n0smax(1:5),t0s(1:5,0:nr),xn0s(1:5,0:nr),&
         capts(1:5,0:nr),capns(1:5,0:nr),vpars(1:5,0:nr),&
         vparsp(1:5,0:nr),tgis(1:5))

    !interface with experiment
    allocate(rdata(nrgem,nthgem),zdata(nrgem,nthgem),rgem(nrgem),psigem(nrgem),qgem(nrgem),fgem(nrgem), &
         tigem(nrgem),tegem(nrgem),nigem(nrgem),negem(nrgem),itoxgc(nrgem),psi_nor_output(nrgem))
    
    allocate(psinxgc(nrxgc),nexgc(nrxgc),texgc(nrxgc),tixgc(nrxgc))
    
    allocate(e1gx(0:nr,0:ntheta),e1gy(0:nr,0:ntheta),e2gx(0:nr,0:ntheta),e2gy(0:nr,0:ntheta), &
         bdge1gx(0:nr,0:ntheta),bdge1gy(0:nr,0:ntheta),bdge2gx(0:nr,0:ntheta),bdge2gy(0:nr,0:ntheta))
#ifdef GPU
!$acc enter data create( bfld,qhat,radius,gr,gth,grdgt,grcgt)
!$acc enter data create( gxdgy,dydr,dbdr,dbdth,jacob)
!$acc enter data create( yfn,hght,thflx) 
!$acc enter data create( psi)
!$acc enter data create( f,psip,sf,jacoba,jfn,zfnth,thfnz)
!$acc enter data create( t0i,t0e,t0b,t0c,t0ip,t0ep,t0bp,t0cp)
!$acc enter data create( xn0i,xn0e,xn0c,xn0b,xn0ip,xn0ep,xn0bp)
!$acc enter data create( xn0cp,vpari,vparc,vparb)
!$acc enter data create( vparip,vparcp,vparbp)
!$acc enter data create( capti,capte,captb,captc,capni,capne)
!$acc enter data create( capnb,capnc,zeff,nue0,phinc,phincp)
!$acc enter data create( er,upari,dipdr)
!$acc enter data create( e1gx,e1gy,e2gx,e2gy,bdge1gx,bdge1gy,bdge2gx,bdge2gy)
!$acc enter data create( psip2)
!$acc enter data create( curvbz,srbr,srbz,thbr,thbz,prsrbr,prsrbz,pthsrbr,pthsrbz,bdcrvb)
!$acc enter data create( t0s,xn0s,capts,capns,vpars,vparsp)
!$acc enter data create( cn0s,n0smax,tgis)
#endif


    !Normalization
    e = 1.6022e-19
    proton = 1.6720e-27
    Bu = 5.3
    Tu = 2.845593e+04*e
    omegau = e*Bu/proton
    frequ = omegau
    rhoi=sqrt(tu/(mimp*proton))/(abs(chgi)*e*Bu/(mimp*proton)) !yjhu added
    vu = sqrt(Tu/proton)
    xu = proton*vu/(e*Bu)
    nu = 8.441875e19
    beta = 4*3.14159*1e-7*nu*Tu/Bu**2 
    debye = (8.85e-12*Tu/(nu*e**2)) / xu**2
    !write(*,*) 'betaU=', beta

    open(113,file='profiles-1d.dat',status='old',action='read')
    do i = 1,nrgem
       read(113,*)rgem(i),qgem(i),fgem(i) ,negem(i),nigem(i),tegem(i),tigem(i),dum,dum,dum,dum,psi_nor_output(i)
    end do
    close(113)
    open(116,file='rdata.dat',status='old',action='read')
    do i = 1,nrgem
       read(116,*)(rdata(i,j),j=1,nthgem)
    end do
    close(116)
    open(117,file='zdata.dat',status='old',action='read')
    do i = 1,nrgem
       read(117,*)(zdata(i,j),j=1,nthgem)
    end do
    close(117)
201 format(2e16.9)    
    
    dum=negem(nrgem)
    nigem=nigem/dum
    negem=negem/dum
    dum=tegem(nrgem)
    tegem=tegem/dum
    tigem=tigem/dum
   
   
    !normalize psigem
    psigem=rgem
    dum = psigem(nrgem)
    psigem = psigem/dum
    psinxgc=psigem

    nexgc=negem
    texgc=tegem
    tixgc=tigem
    !goto 100
    !open(113,file='ECBC_den_bet00p500_corrected_gem.dat',status='old',action='read')
    !read(113,*)k
    !do i = 1,nrxgc
    !   read(113,204)psinxgc(i),nexgc(i)
    !end do
    !close(113)
    !open(113,file='ECBC_te_corrected_gem.dat',status='old',action='read')
    !read(113,*)k    
    !do i = 1,nrxgc
    !   read(113,204)psinxgc(i),texgc(i)
    !end do
    !close(113)
    !open(113,file='ECBC_te_corrected_gem.dat',status='old',action='read')
    !read(113,*)k    
    !do i = 1,nrxgc
    !   read(113,204)psinxgc(i),tixgc(i)
    !end do
    !!close(113)

     if(myid==0)then
       open(113,file='1d.dat',status='replace')
       write(113,*)k
       do i=1,nrxgc
          write(113,204)psigem(i),psinxgc(i),tixgc(i),texgc(i),nexgc(i),qgem(i),fgem(i)
       enddo
       close(113)
    endif
    texgc = texgc      !input data in eV. For tu=1000 eV
    tixgc = tixgc
    nexgc = nexgc
   
  

    !from nexgc(1:nrxgc) to negem(1:nrgem)
    do i = 1,nrgem
       do j = 1,nrxgc-1
          if(psinxgc(j)<=psigem(i) .and. psigem(i)<=psinxgc(j+1)) itoxgc(i)=j
       end do
    end do
    do i = 1,nrgem
       j = itoxgc(i)
       dum = psinxgc(j+1)-psinxgc(j) 
       wx0 = (psinxgc(j+1)-psigem(i))/dum
       wx1 = 1-wx0
       tigem(i) = wx0*tixgc(j)+wx1*tixgc(j+1)
       tegem(i) = wx0*texgc(j)+wx1*texgc(j+1)
       negem(i) = wx0*nexgc(j)+wx1*nexgc(j+1)            
    end do

!100 continue
    
    if(myid==0)then
       open(113,file='rz.dat',status='unknown')
       do i = 1,nrgem
                 write(113,205)(rdata(i,j),j=1,nthgem)
                 write(113,205)(zdata(i,j),j=1,nthgem)
       end do
       close(113)
    end if
205 format(1025(1x,e16.9))

    !if(myid==0)then
    !   open(113,file='1d.dat',status='replace')
    !   do i=1,nrxgc
    !      write(113,204)psinxgc(i),tixgc(i!)
    !   enddo
    !    close(113)
    !endif
204 format(5e16.9) 
    rdata = rdata/xu
    zdata = zdata/xu
    rgem = rgem/xu
    fgem = fgem/(xu*bu)
    
    rmaj0 = rdata(1,nthgem/2)
    a = rgem(nrgem)

    routa = r0a+lxa/2
    rina = r0a-lxa/2
    r0 = r0a*a
    rin = rina*a
    rout = routa*a

    pi = atan(1.0)*4
    drgem = rgem(nrgem)/(nrgem-1)
    dthgem = pi*2/(nthgem-1)
    dr = (rout-rin)/nr
    dth = pi*2/ntheta

    do i = 0,nr
       r = rin+i*dr
       j = int(r/drgem)+1
       j = min(j,nrgem-1)
       j = max(j,0)
       if(j>(nrgem-1) .or. j<0)write(*,*)'r,drgem,j', r,drgem,j
       wx0 = (j*drgem-r)/drgem
       wx1 = 1-wx0
       f(i) = wx0*fgem(j)+wx1*fgem(j+1)
       sf(i) = wx0*qgem(j)+wx1*qgem(j+1)
    end do

    f = f*isgnf
    sf = sf*isgnq
    f0 = f(nr/2)
    q0 = sf(nr/2)
    q0abs = abs(q0)
    q0p = (sf(nr/2+1)-sf(nr/2-1))/(2*dr)
    shat0 = r0/q0*q0p
    
    do i = 0,nr
       r = rin+i*dr
       i1 = int(r/drgem)+1
       i1 = min(i1,nrgem-1)
       wx0 = (i1*drgem-r)/drgem
       wx1 = 1-wx0
       do j = 0,ntheta
          th = -pi+j*dth
          if(th<0.)th = th+pi*2
          j1 = int(th/dthgem)+1
          j1 = min(j1,nthgem-1)
          wy0 = (j1*dthgem-th)/dthgem
          wy1 = 1-wy0
          radius(i,j) = wx0*wy0*rdata(i1,j1)+wx0*wy1*rdata(i1,j1+1) &
                       +wx1*wy0*rdata(i1+1,j1)+wx1*wy1*rdata(i1+1,j1+1)
          hght(i,j) = wx0*wy0*zdata(i1,j1)+wx0*wy1*zdata(i1,j1+1) &
                       +wx1*wy0*zdata(i1+1,j1)+wx1*wy1*zdata(i1+1,j1+1)
       end do
    end do

    !make up-down symmettry
    if(isym==1)then   !set upper half to be equal to the lower half
       do i = 0,nr
          do j = 1,ntheta/2
             radius(i,ntheta/2+j) = radius(i,ntheta/2-j)
             hght(i,ntheta/2+j) =2*hght(i,ntheta/2)-hght(i,ntheta/2-j)
          end do
       end do
    end if
    if(isym==2)then   !set lower half to be equal to the upper half
       do i = 0,nr
          do j = 1,ntheta/2
             radius(i,ntheta/2-j) = radius(i,ntheta/2+j)
             hght(i,ntheta/2-j) = 2*hght(i,ntheta/2)-hght(i,ntheta/2+j)
          end do
       end do
    end if
    

    ! compute grad r,grad theta, grdgt,grcgt
    do i = 0,nr
       do j = 0,ntheta
          jleft=j-1
          jright=j+1
          if(j==0) jleft=ntheta-1
          if(j==ntheta) jright=1

          dum1=(hght(i,jright)-hght(i,jleft))/(2*dth)
          dum2=(radius(i,jright)-radius(i,jleft))/(2*dth)
          if(i==0) then !one-side finite difference
             dum3=(hght(i+1,j)-hght(i,j))/(dr)
             dum4=(radius(i+1,j)-radius(i,j))/(dr)
          elseif(i==nr) then !one-side finite difference
             dum3=(hght(i,j)-hght(i-1,j))/(dr)
             dum4=(radius(i,j)-radius(i-1,j))/(dr)
          else !centered finite difference
             dum3=(hght(i+1,j)-hght(i-1,j))/(2*dr)
             dum4=(radius(i+1,j)-radius(i-1,j))/(2*dr)
          endif

          denom = dum4*dum1-dum3*dum2

          srbr(i,j) = dum1/denom
          srbz(i,j) = -dum2/denom
          thbr(i,j) = -dum3/denom
          thbz(i,j) = dum4/denom
          if(denom<0)write(*,*)'denom<0 !', denom
          gr(i,j) = sqrt(dum1**2+dum2**2)/denom !formula checked yjhu
          gth(i,j) = sqrt(dum3**2+dum4**2)/denom !formula checked yjhu
          grdgt(i,j) = (-dum1*dum3-dum2*dum4)/denom**2 !formula checked yjhu
          grcgt(i,j) = 1/denom !formula checked yjhu
       end do
    end do

    gr(0,:) =   2*gr(1,:)-gr(2,:)
    gr(nr,:)=   2*gr(nr-1,:)-gr(nr-2,:)
    gth(0,:) =   2*gth(1,:)-gth(2,:)
    gth(nr,:)=   2*gth(nr-1,:)-gth(nr-2,:)
    grdgt(0,:) =   2*grdgt(1,:)-grdgt(2,:)
    grdgt(nr,:)=   2*grdgt(nr-1,:)-grdgt(nr-2,:)

    !assign T, n profiles 
    do i = 0,nr     
       r = rin+i*dr
       !s = r/a
       j = int(r/drgem)+1
       j = min(j,nrgem-1)
       if(j>(nrgem-1) .or. j<0)write(*,*)'r,drgem,j', r,drgem,j
       wx0 = (j*drgem-r)/drgem
       wx1 = 1-wx0
       !psin = wx0*psigem(j)+wx1*psigem(j+1)
       t0i(i) = wx0*tixgc(j)+wx1*tixgc(j+1) 
       t0e(i) = wx0*texgc(j)+wx1*texgc(j+1)

       xn0e(i) = wx0*nexgc(j)+wx1*nexgc(j+1)
       xn0i(i) = xn0e(i) !wx0*nigem(j)+wx1*nigem(j+1)
       phincp(i) = 0.
       nue0(i) = 1.
       zeff(i) = 1.
    end do


    !Gradients of T and n profiles for global simulation. Flux-tube parameters done at the end
    capti = 0.
    capte = 0.
    captc = 0.
    capni = 0.
    capne = 0.
    capnc = 0.
    do i = 1,nr-1
       r = rin+i*dr
       s = r/a
       !cutoff=1.-exp(-((routa-s)/delrao)**2)-exp(-((s-rina)/delrai)**2) !yjhu: reduce the gradient drive near the radial boundaries to avoid possbile boundary effects
       cutoff=1. !yjhu
       t0ip(i) = (t0i(i+1)-t0i(i-1))/(2*dr)
       t0ep(i) = (t0e(i+1)-t0e(i-1))/(2*dr)
       capti(i) = -t0ip(i)/t0i(i)*cutoff
       capte(i) = -t0ep(i)/t0e(i)*cutoff

       !       dum = (t0c(i+1)-t0c(i-1))/(2*dr)
       !       captc(i)=  -dum/t0c(i)*cutoff
       dum = (xn0e(i+1)-xn0e(i-1))/(2*dr)
       capne(i) = -dum/xn0e(i)*cutoff
       dum = (xn0i(i+1)-xn0i(i-1))/(2*dr)
       capni(i) = -dum/xn0i(i)*cutoff 

       !       dum = (xn0c(i+1)-xn0c(i-1))/(2*dr)
       !       capnc(i) = -dum/xn0c(i)*cutoff 
    end do

    ! compute psip(r), from f and q
    do i = 0,nr
       dum = 0.
       do j = 0,ntheta-1
          dum = dum+dth/radius(i,j)/grcgt(i,j)
       end do
       psip(i) = f(i)/2/pi/sf(i)*dum !yjhu, formula checked, correct
    end do
    bunit = q0/r0*psip(nr/2) !used in GYRO to define gyro's rhos=mi*sqrt(Te/mi)/(qi*bunit)
#ifndef __GEM_XGC_COUPLING 
    !compute psi(r) from psip. psi(:) is used only in the toroidal canonical momentum pzeta, the initial condition psi(0) not needed
    dum = 0.
    psi(0) = 0.
    do i = 1,nr
       psi(i) = psi(i-1)+dr*(psip(i-1)+psip(i))/2
    end do
#else
    do i = 0,nr                               
       psi(i)=psi_nor_output(31+i)                         
    enddo  
#endif

    do i = 0,nr
       do j = 0,ntheta
          bfld(i,j) = sqrt((f(i)/radius(i,j))**2+ &
               (psip(i)/radius(i,j)*gr(i,j))**2) !formula checked yjhu

          qhat(i,j) = f(i)/radius(i,j)/(psip(i)*grcgt(i,j))  !yjhu, general formula, checked
       end do
    end do

    do i = 1,nr-1
       do j = 1,ntheta-1
          dbdr(i,j) = (bfld(i+1,j)-bfld(i-1,j))/(2*dr)
          dbdth(i,j) = (bfld(i,j+1)-bfld(i,j-1))/(2*dth)
       end do
       dbdr(i,0) = (bfld(i+1,0)-bfld(i-1,0))/(2*dr)
       dbdr(i,ntheta) = dbdr(i,0)
       dbdth(i,0) = (bfld(i,1)-bfld(i,ntheta-1))/(2*dth)
       dbdth(i,ntheta) = dbdth(i,0)
    end do
    do j = 0,ntheta
       dbdr(0,j) = dbdr(1,j)
       dbdr(nr,j) = dbdr(nr-1,j)
       dbdth(0,j) = dbdth(1,j)
       dbdth(nr,j) = dbdth(nr-1,j)
    end do

    !compute dydr(r,theta)
    do i = 1,nr-1
       do j = 0,ntheta
          dqhdr(i,j) = (qhat(i+1,j)-qhat(i-1,j))/(2*dr)
       end do
    end do
    do j = 0,ntheta
       dqhdr(0,j) = (qhat(1,j)-qhat(0,j))/dr
       dqhdr(nr,j) = (qhat(nr,j)-qhat(nr-1,j))/dr
    end do


    do i = 0,nr
       yfn(i,ntheta/2) = 0.
       dydr(i,ntheta/2) = 0.
       dum = 0.
       dum1 = 0.
       do j = ntheta/2+1,ntheta
          dum = dum+(dqhdr(i,j-1)+dqhdr(i,j))*dth/2
          dydr(i,j) = r0/q0*dum
          dum1 = dum1+r0/q0*(qhat(i,j-1)+qhat(i,j))*dth/2
          yfn(i,j) = dum1
       end do
       dum = 0.
       dum1 = 0.
       do j = ntheta/2-1, 0, -1 !yjhu
          dum = dum-(dqhdr(i,j+1)+dqhdr(i,j))*dth/2
          dydr(i,j) = r0/q0*dum
          dum1 = dum1-r0/q0*(qhat(i,j+1)+qhat(i,j))*dth/2
          yfn(i,j) = dum1
       end do
    end do
    
    !compute the flux coordinate theta
    do i = 0,nr
       thflx(i,ntheta/2) = 0.
       dum = 0.
       do j = ntheta/2+1,ntheta !ntheta/2 corresponds to theta=0 in gcrz.f90
          dum = dum+(qhat(i,j-1)+qhat(i,j))*dth/2
          thflx(i,j) = dum/sf(i) !correct, yjhu checked
       end do
       dum = 0.
       do j = ntheta/2-1, 0, -1 !yjhu
          dum = dum-(qhat(i,j+1)+qhat(i,j))*dth/2
          thflx(i,j) = dum/sf(i) !correct, yjhu checked
       end do
    end do

    ! compute gxdgy
    jacmax = 0.
    do i = 0,nr
       dum = 0.
       do j = 0,ntheta
          gxdgy(i,j) = dydr(i,j)*gr(i,j)**2+r0/q0*qhat(i,j)*grdgt(i,j)
          jacob(i,j) = 1./(r0*rmaj0/radius(i,j)*grcgt(i,j)) !yjhu, formula checked
          if(jacob(i,j)>jacmax)jacmax = jacob(i,j)
          if(j<ntheta)dum = dum+jacob(i,j)
       end do
       jacoba(i) = dum/ntheta
    end do
    
    !compute cn0e,cn0i,cn0b,cn0e
    cn0e = 0.
    cn0i = 0.
    cn0c = 0.
    dum = 0.
    do i = 0,nr-1
       dum = dum+(jacoba(i)+jacoba(i+1))/2.0
    end do
    do i = 0,nr-1
       cn0e = cn0e+(xn0e(i)+xn0e(i+1))/2*(jacoba(i)+jacoba(i+1))/2.0
       cn0i = cn0i+(xn0i(i)+xn0i(i+1))/2*(jacoba(i)+jacoba(i+1))/2.0
       cn0c = cn0c+(xn0c(i)+xn0c(i+1))/2*(jacoba(i)+jacoba(i+1))/2.0
    end do
    cn0e = cn0e/dum
    cn0i = cn0i/dum
    cn0c = cn0c/dum
    do i = 0,nr
       xn0e(i) = xn0e(i)/cn0e
       xn0i(i) = xn0i(i)/cn0i 
       xn0c(i) = xn0c(i)/cn0c
    end do
    !yjhu: cn0i is (approximate) volume average of equilibrium density xn0i (xn0i is in unit of nu and so is cn0i).
    !cn0i is used to normalize the density profile (originally in unit of nu). The pertubed density is therefore also
    !in unit of cn0i. And later in grid1() subroutine, cn0i is used to multify the perturbed density,
    !and thus the resulting perturbed density is in unit of nu.

    !assign value to xn0s...
    xn0s(1,:) = xn0i(:)
    xn0s(2,:) = xn0c(:)
    !      xn0s(3,:) = xn0b(:)

    t0s(1,:) = t0i(:)
    t0s(2,:) = t0c(:)
    !      t0s(3,:) = t0b(:)

    tgis(1) = t0s(1,1)
    tgis(2) = t0s(2,1)
    tgis(3) = t0s(3,1)
    tge = t0e(1)

    capts(1,:) = capti(:)
    capts(2,:) = captc(:)
    !      capts(3,:) = captb(:)

    capns(1,:) = capni(:)
    capns(2,:) = capnc(:)
    !      capns(3,:) = capnb(:)

    cn0s(1) = cn0i
    cn0s(2) = cn0c
    cn0s(3) = cn0b

    vpars(1,:) = 0. !vpari(:)
    vpars(2,:) = 0. !vparc(:)
    !      vpars(3,:) = vparb(:)

    vparsp(1,:) = 0. !vparip(:)
    vparsp(2,:) = 0. !vparcp(:)
    !      vparsp(3,:) = vparbp(:)


    !compute jfn(theta)
    do j = 0,ntheta
       dum = 0.
       do i = 0,nr-1
          dum = dum+(jacob(i,j)+jacob(i+1,j))/2.0
       end do
       jfn(j) = dum/nr
    end do
    dum = 0.
    do j = 0,ntheta-1
       dum = dum+(jfn(j)+jfn(j+1))/2
    end do
    dum = dum/ntheta
    do j = 0,ntheta
       jfn(j) = dum/jfn(j)
    end do
    jfn = 1.

    !bstar effects
    ! compute psip2(r), from psip(r)
    do i = 1,nr-1
       psip2(i) = (psip(i+1)-psip(i-1))/(2*dr)
    end do
    psip2(0) = psip2(1)
    psip2(nr) = psip2(nr-1)

    ! compute prsrbr, prsrbz, pthsrbr,pthsrbz
    do j = 0,ntheta
       do i = 1,nr-1
          prsrbr(i,j) = (srbr(i+1,j)-srbr(i-1,j))/(2.*dr)
          prsrbz(i,j) = (srbz(i+1,j)-srbz(i-1,j))/(2.*dr)
       end do
       prsrbr(0,j) = (srbr(1,j)-srbr(0,j))/dr
       prsrbr(nr,j) = (srbr(nr,j)-srbr(nr-1,j))/dr
       prsrbz(0,j) = (srbz(1,j)-srbz(0,j))/dr
       prsrbz(nr,j) = (srbz(nr,j)-srbz(nr-1,j))/dr         
    end do
    do i = 0,nr
       do j = 1,ntheta-1
          pthsrbr(i,j) = (srbr(i,j+1)-srbr(i,j-1))/(2.*dth)
          pthsrbz(i,j) = (srbz(i,j+1)-srbz(i,j-1))/(2.*dth)
       end do
       pthsrbr(i,0) = (srbr(i,1)-srbr(i,0))/dth
       pthsrbr(i,ntheta) = (srbr(i,ntheta)-srbr(i,ntheta-1))/dth
       pthsrbz(i,0) = (srbz(i,1)-srbz(i,0))/dth
       pthsrbz(i,ntheta) = (srbz(i,ntheta)-srbz(i,ntheta-1))/dth
    end do

    ! compute curvbz
    do i = 0,nr
       do j = 0,ntheta
          dum1 = prsrbz(i,j)*srbz(i,j)+pthsrbz(i,j)*thbz(i,j)
          dum2 = prsrbr(i,j)*srbr(i,j)+pthsrbr(i,j)*thbr(i,j)
          curvbz(i,j) = psip(i)*(dum1/radius(i,j)-srbr(i,j)/radius(i,j)**2 &
               +dum2/radius(i,j))
          bdcrvb(i,j) = f(i)/(bfld(i,j)**2*radius(i,j))*(psip2(i)*(gr(i,j))**2/radius(i,j)+curvbz(i,j))                       
       end do
    end do

    do i = 0,nr
       do j = 0,ntheta
          dum=sqrt(srbr(i,j)**2+srbz(i,j)**2)
          e1r(i,j) = srbr(i,j)/dum
          e1z(i,j) = srbz(i,j)/dum
          e2r(i,j) = -e1z(i,j)*f(i)/(radius(i,j)*bfld(i,j))
          e2z(i,j) = e1r(i,j)*f(i)/(radius(i,j)*bfld(i,j))
          e2zet(i,j) = (-e1r(i,j)*srbr(i,j)-e1z(i,j)*srbz(i,j))*psip(i)/(radius(i,j)*bfld(i,j))
       end do
    end do
    call bdgrad(e1r,bdge1r)
    call bdgrad(e1z,bdge1z)
    call bdgrad(e2r,bdge2r)
    call bdgrad(e2z,bdge2z)
    call bdgrad(e2zet,bdge2zet)      

    do i = 0,nr
       do j = 0,ntheta
          e1gx(i,j) = sqrt(srbr(i,j)**2 + srbz(i,j)**2)
          e1gy(i,j) = dydr(i,j)*e1gx(i,j) + (r0/q0)*qhat(i,j)*grdgt(i,j)/e1gx(i,j)
          e2gx(i,j) = e2r(i,j)*srbr(i,j) + e2z(i,j)*srbz(i,j)
          dum = e2r(i,j)*thbr(i,j)+e2z(i,j)*thbz(i,j)
          e2gy(i,j) = dydr(i,j)*e2gx(i,j)+(r0/q0)*qhat(i,j)*dum-e2zet(i,j)*r0/(q0*radius(i,j))
          bdge1gx(i,j) = bdge1r(i,j)*srbr(i,j)+bdge1z(i,j)*srbz(i,j)
          dum = bdge1r(i,j)*thbr(i,j)+bdge1z(i,j)*thbz(i,j)
          bdge1gy(i,j) = dydr(i,j)*bdge1gx(i,j)+(r0/q0)*qhat(i,j)*dum-(r0/q0)*e1r(i,j)*f(i)/(bfld(i,j)*radius(i,j)**3)
          bdge2gx(i,j) = (bdge2r(i,j)-e2zet(i,j)*f(i)/(bfld(i,j)*radius(i,j)**2))*srbr(i,j)+bdge2z(i,j)*srbz(i,j)
          dum = (bdge2r(i,j)-e2zet(i,j)*f(i)/(bfld(i,j)*radius(i,j)**2))*thbr(i,j)+bdge2z(i,j)*thbz(i,j)
          bdge2gy(i,j) = dydr(i,j)*bdge2gx(i,j)+(r0/q0)*qhat(i,j)*dum-(bdge2zet(i,j)+e2r(i,j)*f(i)/(bfld(i,j)*radius(i,j)**2))*r0/(q0*radius(i,j))
       end do
    end do
    
    
    dipdr = 0.
    do i = 1,nr-1
       dipdr(i) = (f(i+1)-f(i-1))/(2*dr)
    end do
!    bdcrvb = 0.
!    curvbz = 0.
!    psip2 = 0.

    if(myid==0)then
       open(11,file='xpp',status='replace')
       write(11,*)'xu,omegau,vu,debye = ', xu, omegau, vu, debye
       write(11,*)'q0,q0p,shat0,rmaj0,dr,a,rin,rout=',q0,q0p,shat0,rmaj0*xu,dr*xu,a*xu,rin*xu,rout*xu
       write(11,*)'capte(nr/2),capti(nr/2),capne(nr/2),capni(nr/2)',capte(nr/2)*rmaj0,capti(nr/2)*rmaj0,capne(nr/2)*rmaj0,capni(nr/2)*rmaj0
       write(11,*)'bdcrvb=', maxval(abs(bdcrvb)),minval(abs(bdcrvb))
       write(11,*)'curvbz=', maxval(abs(curvbz)),minval(abs(curvbz))
       write(11,*)'psip2=', maxval(abs(psip2)),minval(abs(psip2))
       write(11,*)'dipdr=', maxval(abs(dipdr)),minval(abs(dipdr))
       write(11,*)'gr=', maxval(abs(gr)),minval(abs(gr))
       write(11,*)'gth=', maxval(abs(gth)),minval(abs(gth))
       write(11,*)'grcgt=', maxval(abs(grcgt)),minval(abs(grcgt))                     
       do j = 0,ntheta
          !            write(11,10)j,jacob(nr/2,j),bfld(nr/2,j),dbdr(nr/2,j),gr(nr/2,j),dydr(nr/2,j),gxdgy(nr/2,j),radius(nr/2,j),hght(nr/2,j)
       end do
10     format(1x,i5,10(1x,e16.9))
       do i = 0,nr
          write(11,10)i,sf(i),f(i),t0i(i),t0e(i),xn0i(i),xn0e(i),capti(i)*rmaj0,capte(i)*rmaj0,capni(i)*rmaj0,capne(i)*rmaj0
       end do
       call flush(11)
       close(11)
    end if
      
  end subroutine new_equil

  subroutine bdgrad(u,v)
    implicit none
    real(8) :: pi,pi2,r,th
    integer :: i,j,k,m
    real(8) :: dum,x,denom,dum1,dum2,dum3,dum4
    real(8) :: u(0:nr,0:ntheta),v(0:nr,0:ntheta),w(0:nr,0:ntheta)


    do i = 0,nr
       do j = 1,ntheta-1
          w(i,j) = (u(i,j+1)-u(i,j-1))/(2*dth)
       end do
       w(i,0) = (u(i,1)-u(i,ntheta-1))/(2*dth)
       w(i,ntheta) = w(i,0)
    end do

    do i = 0,nr
       do j = 0,ntheta
          v(i,j) = psip(i)/(bfld(i,j)*radius(i,j))*w(i,j)*grcgt(i,j)
       end do
    end do
  end subroutine bdgrad
  
END MODULE gem_equil
