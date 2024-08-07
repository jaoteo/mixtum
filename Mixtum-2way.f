!	2-way admixture problem: AA=alpha*C1+(1-alpha)*C2.
!	AA: hybrid, C1 & C2: contributors, alpha: mixing proportion.
!	Computes alpha estimate.
!	Computes pre-JL and post-JL estimates of alpha, ancestral mixing proportion.
!	Computes pre-JL and post-JL estimates of angle, admixture tests.
!	Input: 'selectedPopulations.dat', whose first line reads XYZ (the name of the 
!	Allele Frequencies dataset, 'XYZ-Allele-Frequencies.dat'). The next lines are 
!	the population labels and are read.
!	Header in 'XYZ-Allele-Frequencies.dat' is population names in each column:
!	Column 1=hybrid, Column 2=Parental-1, Column 3=Parental-2,
!	Columns 4...kpop = auxiliary populations.
!	kpopmax: max number of populations in allele freqs dataset (=max number of columns).
!	kpop: real number of populations in the dataset (=number of data columns).
!	kgenomax: max number of SNP to read.
!	Outcome: 'Mixtum-f4-points.dat': f_4's outcomes in columns ready to be plotted, and f4-ratios,
!	also. Summary of results at the bottom of the file.
!	'Mixtum-f4-alphas-OK.dat': only f4-ratios in the interval [0,1] with 
!	the corresponding aux pops pair.
!	'Mixtum-out-Brief.dat': summary (incremental, every new run is added).
! 	Note: removing 'implicit double precision (a-h,o-z)' allows to increase array dimensions.
!	Compile (Intel, Windows):  ifort  MixtuM.f
!	Compile (gfortran, Windows): gfortran -o mixtum mixtum.f
!	Execute: Mixtum 
!	CPU time: few secs.
!	****************************************************************************************
!	oteo@uv.es, Departament de Fisica Teorica, Universitat de Valencia, Spain.
!	July, 2024
!	****************************************************************************************
	implicit double precision(a-h,o-z)
	parameter(kgenomax=597562,kpopmax=42,kx=(kpopmax-3)*(kpopmax-4)/2)
cc	parameter(kgenomax=2000,kpopmax=43,kx=(kpopmax-3)*(kpopmax-4)/2)	!for testing
	dimension f(kgenomax,kpopmax), AA(kgenomax)
	dimension C1(kgenomax),C2(kgenomax)
	dimension auxi(kgenomax),auxj(kgenomax)
	dimension vx(kx),vy(kx),vw(kx)
	dimension vx0(kx),vy0(kx),vw0(kx)
	dimension f4r(kx)
	character*100  pops(kpopmax),res,dataset
	character*30 auxpair
	character*150 outfreqs
c
	open(unit=9, file='selectedPopulations.dat',status='old')	!driver
	read(9,*) dataset
	outfreqs=
     &   adjustl(adjustr(dataset)//'-Allele-Frequencies.dat')
	do k=1,10*kpopmax
	read(9,*,end=555) res
	enddo
555	kpop=k-1
	if(kpop.ge.kpopmax) then 
	print *, 
     &   'Warning: number of pops read has reached the max (= kpopmax)'
	print *, 'allowed in the code. Possible missing pops.'
	endif
	open(unit=10,file=outfreqs) 		!allele freq	
	open(unit=40,file='Mixtum-f4-points.dat')	!outcome for 2D plot
	open(unit=50,file='Mixtum-out-Brief.dat',access='append')	!outcome, incremental
	open(unit=60,file='Mixtum-f4-alphas-OK.dat')	!outcome, if  0<alpha<1
c
	rad=180./(4*atan(1.))   !Pi radians in degrees
	write(40,987)		!header
987	format('       f4primeAB       f4primeXB',
     &   '              f4AB            f4XB      f4-ratio'
     &   '          Aux1   Aux2')
c
	read(10,*) pops(1:kpop)		!reads header (populations) in the dataset
	do k=1,kgenomax
	read(10,*,end=777) f(k,1:kpop)
	enddo
	kgeno=kgenomax
	goto 778
777	kgeno=k-1			!number of SNPs
778	continue
	AA(1:kgeno)=f(1:kgeno,1)	!hybrid
	C1(1:kgeno)=f(1:kgeno,2)	!contributor 1
	C2(1:kgeno)=f(1:kgeno,3)	!contributor 2
c
	write(50,*) 'Admixture model:'
	write(50,756) 
     &   adjustl(adjustr(pops(1))//'='//adjustl(adjustr(pops(2))
     &   //'+'//adjustl(pops(3))))
 756	format(1x,A90)
	write(50,*) 'Auxiliary populations:'
	write(50,766) pops(4:kpop)
766	format(200(3(A30,1x),/))
	print *, 'Admixture model: ',
     &   adjustl(trim(adjustr(pops(1))//'='//adjustl(adjustr(pops(2))
     &   //'+'//adjustl(pops(3)))))
	print *, 'Auxiliary populations: '
	print 505,pops(4:kpop)
505	format(5(A20,2x))
c
	!  f-statistics are computed as dot-products
	icomb=0			!aux pops pairs counter
	krat=0			![0,1] f4-ratio counter
	tstudent=1.960		!t-Student factor 95% CI
	do 1 i=4,kpop-1
	auxi(1:kgeno)=f(1:kgeno,i)
	do 2 j=i+1,kpop
	icomb=icomb+1
	auxj(1:kgeno)=f(1:kgeno,j)		!auxiliary populations
	x0=dot_product(C1-C2,auxi-auxj)		!f_4's
	y0=dot_product(AA-C2,auxi-auxj)
	w0=dot_product(AA-C1,auxi-auxj)		!used in post-JL angle
	f2=dot_product(auxi-auxj,auxi-auxj)	!f_2's
	rf2=sqrt(f2)
	x=x0/rf2				!f_4^prime's (renormalised)
	y=y0/rf2
	w=w0/rf2
	vx(icomb)=x
	vy(icomb)=y
	vw(icomb)=w
	x0=x0/kgeno				!f_4 (standard)
	y0=y0/kgeno
	w0=w0/kgeno
	vx0(icomb)=x0
	vy0(icomb)=y0
	vw0(icomb)=w0
	auxpair=adjustl(adjustr(pops(i))//'  '//adjustl(pops(j)))	!aux pops label
	aaa=y0/x0				!f4-ratio 
	write(40,123) x,y,x0,y0,aaa,auxpair	!data to 2D plot
	if(aaa.le.1) then
	if(aaa.gt.0) then 
	write(60,*) real(aaa),'   ',auxpair !fair f4-ratios
	krat=krat+1
	f4r(krat)=aaa		!fair f4-ratios statistics
	endif
	endif
2	continue
1	continue
123	format(2E16.6,2x,2E16.6,4x,f10.4,a40)
	!f4-ratio statistical estimate
	f4rMean=9.9999		!when f4-ratio fails
	f4rSD=9.9999
	if(krat.le.1) goto 8765
	f4rmean=sum(f4r(1:krat))/krat		!mean
	f4r(1:krat)=f4r(1:krat)-f4rMean
	f4rSD=sqrt(dot_product(f4r,f4r)/krat)	!standard deviation
8765	print *, 'Auxiliary pairs:',icomb
	print *, 'SNPs read: ', kgeno
	write(40,*)
	write(40,*) 'Auxiliary pairs:',icomb
	write(40,*) 'SNPs read: ', kgeno
	write(50,*)
	write(50,*) 'Auxiliary pairs:',icomb
	write(50,*) 'SNPs read: ', kgeno
	print *
	write(50,*)
	write(40,*)
	res=
     &   adjustr(pops(1))//'='//adjustl(adjustr(pops(2))
     &   //'+'//adjustl(pops(3)))
	write(40,606) adjustl(adjustr(pops(1))//'='
     &   //adjustl(adjustr(pops(2))//'+'//adjustl(pops(3))))
606	format('Admixture model: ',A100)
c
	! admixture angles: phipre & phipost JL
	t1=dot_product(AA-C1,AA-C2)
	f3test=t1/kgeno		!f_3(c1,c2;aa)<0? standard admixture test
	t2=dot_product(AA-C1,AA-C1)
	t3=dot_product(AA-C2,AA-C2)
	phipre=t1/sqrt(t2*t3)		!cos angle pre-JL
	aphipre=rad*acos(phipre)	!angle pre-JL
	zphipre=aphipre*100/180.	!% of angle /180
	t1=dot_product(vy,vw)
	t2=dot_product(vy,vy)
	t3=dot_product(vw,vw)
	phipost=t1/sqrt(t2*t3)		!cos angle post-JL
	aphipost=rad*acos(phipost)	!angle post-JL
	zphipost=aphipost*100./180.	!% of angle /180
	print 400, phipre,rad*acos(phipre)
	write(40,400)  phipre,aphipre,zphipre
	write(50,400)  phipre,aphipre,zphipre
400	format('cos pre-JL: ',F8.4,'  --->  ', 
     &   'angle pre-JL: ',F8.2,' deg',4x,'vs. 180 deg:',f5.1,'%')
	print 401, phipost,rad*acos(phipost)
	write(40,401)  phipost,aphipost,zphipost
	write(50,401)  phipost,aphipost,zphipost
401	format('cos post-JL:',F8.4,'  --->  ',
     &   'angle post-JL:',F8.2,' deg',4x,'vs. 180 deg:',f5.1,'%')
c
	! alpha estimates: post JL (f4' & f4 slopes)
	aPre=dot_product(C1-C2,AA-C2)/dot_product(C1-C2,C1-C2) !alpha pre-JL
	aPost=dot_product(vy,vx)/dot_product(vx,vx)  !alpha post-JL f_4^prime
	aPostNR=dot_product(vy0,vx0)/dot_product(vx0,vx0) !alpha post-JL f_4 standard
	! alpha C.I.: pre & post JL (slope CI's)
	t1=dot_product(aPost*vx-vy,aPost*vx-vy)/(krat-2)
	vxmean=sum(vx(1:krat))/krat
	t2=dot_product(vx-vxmean,vx-vxmean)
	sig=sqrt(t1/t2)*tstudent	!CI
	t1=dot_product(aPost*vx0-vy0,aPost*vx0-vy0)/(krat-2)
	vxmean=sum(vx0(1:krat))/krat
	t2=dot_product(vx0-vxmean,vx0-vxmean)
	sigNR=sqrt(t1/t2)*tstudent	!CI
c	!writes summary
	print 402, aPre,aPost,sig,aPostNR,sigNR
	write(40,402) aPre,aPost,sig,aPostNR,sigNR
	write(50,402) aPre,aPost,sig,aPostNR,sigNR
402	format('alpha-pre-JL:     ',F8.4,/,'alpha-post-JL:    ',F8.4,
     &   2x,F8.4,' (95% CI)',4x,
     &   '(f_4-prime,renormalised)',/,'alphaNR-post-JL:',F10.4,
     &   2x,F8.4,' (95% CI)',4x,'(f_4, standard)')
	write(40,403) f4rmean,f4rSD*tstudent,krat,f3test
	write(50,403) f4rmean,f4rSD*tstudent,krat,f3test
	print 403, f4rmean,f4rSD,krat,f3test
403	format('f4-ratio mean if [0,1]: ',f8.4,' +/-',f8.4,2x,'(95% CI)',
     &     2x,i6,1x,'cases',/,
     &    'standard admixture test, f_3(c1,c2;x)<0?:', f10.6)
	write(50,*)
c
	print *
	print *, 'done'
	CALL CPU_TIME(TIME)
	print *,'Total CPU= ',real(TIME),'s.'
	end
	