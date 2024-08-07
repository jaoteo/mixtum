! 	ifort/heap-arrays Mixtum-3way.f
c	implicit double precision (a-h,o-z)
	parameter(kgenoz=597562,kpop=43)
c	parameter(kgenoz=200000,kpop=43)
	dimension  freqs(kpop,kgenoz)
	dimension AA(kgenoz),C1(kgenoz),C2(kgenoz),C3(kgenoz)
	dimension A(kgenoz),B(kgenoz),C(kgenoz),D(kgenoz)
	dimension u(kgenoz),v(kgenoz),w(kgenoz),vol(3,3)
	dimension uv(kgenoz),uw(kgenoz),vw(kgenoz)
	dimension uJL((kpop-4)*(kpop-5)/2),vJL((kpop-4)*(kpop-5)/2)
	dimension wJL((kpop-4)*(kpop-5)/2)
	dimension uvJL((kpop-4)*(kpop-5)/2),uwJL((kpop-4)*(kpop-5)/2)
	dimension vwJL((kpop-4)*(kpop-5)/2)
	character*15  pops(kpop)
c
c	open(unit=10,file='v37.2.1240K_HumanOrigins_unzip.ind')
c	open(unit=10,file='v37.2.1240K_HumanOrigins_unzip.geno')
c	open(unit=10,file='v37.2.1240K_HOP-Freqs-11pops.dat')
	open(unit=10,file=
     &   'v37.2.1240K_HOP-Allele-Freqs-CLM=IBS+YRI+PEL.dat')
	open(unit=50,file=
     &  'v37.2.1240K_HOP-3way-MixtuM.dat',access='append')
	open(unit=51,file=
     &   'v37.2.1240K_HOP-3way-Mixtum-Points.dat')
c
	!Model
	!hybrid: 	pops(1)
	!contributor 1: pops(1)
	!contributor 2: pops(2)
	!contributor 3: pops(3)
c
	write(51,*) ' f4Pacij        f4Pbcij        f4Pxcij'
	pi=4.d0*atan(1.d0)
	uno=1.d0
	cero=0.d0
	read(10,*) pops(1:kpop)	!population labels
	print *, pops(1:kpop)
	write(50,*) pops(1:kpop)
	do 44 k=1,kgenoz
	read(10,654,end=777) freqs(1:kpop,k)
44	continue
654	format(10000e14.6)
	kgeno=k+1
	goto 778
777	print *, 'Freqs. reading error at SNP:', k
	stop
778	continue
	kgeno=k-1	!#SNPs
	print *, 'SNPs read:',kgeno
	write(50,*) 'SNPs read:',kgeno
	sqgeno=sqrt(kgeno*1.)
	do k=1,kgeno		!admixture model
	AA(k)=freqs(1,k)
	C1(k)=freqs(2,k)
	C2(k)=freqs(3,k)
	C3(k)=freqs(4,k)
	enddo
c
	print *, 'MODEL: ',adjustr(pops(1)),' = ',pops(2:4)
	write(50,*) 'MODEL: ',adjustr(pops(1)),' = ',pops(2:4)
c	! 3-way admixture coefficients in full phase space
	call coeffs(kgeno,AA(1:kgeno),C1(1:kgeno),C2(1:kgeno),
     &    C3(1:kgeno),aFull,bFull,cFull)	!pre-JL
	print *, '   alpha     beta      gamma'
	write(50,*) '   alpha     beta      gamma'
	print 987, aFull,bFull,cFull
	write(50,987) aFull,bFull,cFull
987	format(3f10.4,2x,':(',2x,'Full Phase Space')
	call coeffsJL(kgeno,kpop,AA(1:kgeno),C1(1:kgeno),
     &  C2(1:kgeno),C3(1:kgeno),freqs(1:kpop,1:kgeno),aJL,bJL,gJL)	!post-JL
	print 986, aJL,bJL,gJL
	write(50,986)  aJL,bJL,gJL
986	format(3f10.4,2x,':)',2x,'JL subspace')
c
c	! tetrahedron volume in Phase Space
	u=AA-C1		!edges
	v=AA-C2
	w=AA-C3
	uv=C1-C2
	uw=C1-C3
	vw=C2-C3
	vol(1,1)=dot_product(u,u)	!squared edges
	vol(2,2)=dot_product(v,v)
	vol(3,3)=dot_product(w,w)
	vol(1,2)=dot_product(uv,uv)
	vol(2,1)=vol(1,2)
	vol(1,3)=dot_product(uw,uw)
	vol(3,1)=vol(1,3)
	vol(2,3)=dot_product(vw,vw)
	vol(3,2)=vol(2,3)
c		print *, vol(1:3,1:3)
c	! Cayley
	d12=vol(1,2)	! squared distances
	d13=vol(1,3)
	d23=vol(2,3)
	d14=vol(1,1)
	d24=vol(2,2)
	d34=vol(3,3)
! Cayley–Menger determinant: 
! https://en.wikipedia.org/wiki/Cayley%E2%80%93Menger_determinant
! 288*vol**2=determinant
!	expression exported from Maple code
      detCM=- 2 * d12 ** 2 * d34 - 2 * d12 * d13 * d23 + 2*d12*d13* 
     &d24 + 2 * d12 * d13 * d34 + 2 * d12 * d14 * d23 - 2 *d12*d14* 
     &d24 + 2 * d12 * d14 * d34 + 2 * d12 * d23 * d34 + 2 *d12*d24*
     &d34 - 2 * d12 * d34 ** 2 - 2 * d13 ** 2 * d24 + 2 *d13*d14*d2
     &3 + 2 * d13 * d14 * d24 - 2 * d13 * d14 * d34 + 2 *d13*d23*d2
     &4 - 2 * d13 * d24 ** 2 + 2 * d24 * d13 * d34 - 2 *d14**2*d23 
     &- 2 * d14 * d23 ** 2 + 2 * d14 * d23 * d24 + 2 * d23*d14*d34-
     & 2 * d23 * d24 * d34	
c	
	if(detCM.lt.cero) print *, vol(1:3,1:3)
	if(detCM.lt.cero) print *, detCM
	volFullCM=sqrt(detCM/288)
c	! Heron formula
        heronCM =d12**2-2*d13 * d12 - 2 * d12 * d23 +d13**2-2*d13
     & * d23 + d23 ** 2
	areabaseFull=sqrt(-heronCM/16)
c	
c	! JL subspace
c
	icomb=1
	do 6660 k=6,kpop
	do 666 j=5,k-1
c
	A(1:kgeno)=freqs(k,1:kgeno)
	B(1:kgeno)=freqs(j,1:kgeno)
c
c	! renormalization
	fab=sqrt(dot_product(A-B,A-B))
c
c	! f4-prime = (AA-C1)*(A-B)/|A-B|/sqrt(kgeno)
	uJL(icomb)=dot_product(u,A-B)/fab	!(AA-C1)*(A-B)/|A-B|
	vJL(icomb)=dot_product(v,A-B)/fab
	wJL(icomb)=dot_product(w,A-B)/fab
	uvJL(icomb)=dot_product(uv,A-B)/fab
	uwJL(icomb)=dot_product(uw,A-B)/fab
	vwJL(icomb)=dot_product(vw,A-B)/fab
	icomb=icomb+1
	write(51,*) real(w(icomb)/sqgeno),real(uw(icomb)/sqgeno),
     &    real(vw(icomb)/sqgeno)
c
666	enddo		!combinat
6660	enddo		!combinat
c
	icomb=icomb-1
c	! tetrahedron volume in JL
	z1=dot_product(ujl,ujl)	!edges
	z2=dot_product(vjl,vjl)
	z3=dot_product(wjl,wjl)
	x1=dot_product(uvjl,uvjl)
	x2=dot_product(uwjl,uwjl)
	x3=dot_product(vwjl,vwjl)
c
	d12=x1		!edge lengths
	d13=x2
	d23=x3
	d14=z1
	d24=z2
	d34=z3
c	! Cayley–Menger determinant: https://en.wikipedia.org/wiki/Cayley%E2%80%93Menger_determinant
c	! 288*vol**2=determinant
	!expression exported from Maple code
      volCM= -2 * d12 ** 2 * d34 - 2 * d12 * d13 * d23+2*d12*d13* 
     &d24 + 2 * d12 * d13 * d34 + 2 * d12 * d14 * d23 -2*d12*d14* 
     &d24 + 2 * d12 * d14 * d34 + 2 * d12 * d23 * d34 +2*d12*d24* 
     &d34 - 2 * d12 * d34 ** 2 - 2 * d13 ** 2 * d24 +2*d13*d14*d2
     &3 + 2 * d13 * d14 * d24 - 2 * d13 * d14 * d34 +2*d13*d23*d2
     &4 - 2 * d13 * d24 ** 2 + 2 * d24 * d13 * d34 - 2*d14**2*d23 
     &- 2 * d14 * d23 ** 2 + 2 * d14 * d23 * d24 + 2*d23*d14*d34-
     & 2 * d23 * d24 * d34	
c	! Heron formula
        heronJL =d12**2-2*d13 * d12 - 2 * d12 * d23 + d13 ** 2 - 2 * d13
     & * d23 + d23 ** 2
	areabaseJL=sqrt(-heronJL/16)
c
	volJLCM=sqrt(volCM/288)
	print *, '                    pre-JL         post-JL'
	print *, 'Volum:           ',real(volFullCM),real(volJLCM)
	write(50,*) '                    pre-JL         post-JL'
	write(50,*) 'Volum:           ',
     &   real(volFullCM),real(volJLCM)
c	! Vol tetrahedron: base*h/3
	print *, 'Area(base):      ',
     &   real(areabaseFull), real(areabaseJL)
	write(50,*) 'Area(base):      ',
     &   real(areabaseFull), real(areabaseJL)
c	! Heights
	hFull=3*volFullCM/areabaseFull
	hJL=3*volJLCM/areabaseJL
	print *, 'Height:          ', real(hFull), real(hJL)
	write(50,*) 'Height:          ', real(hFull), real(hJL)
c		
	if(volFullCM.lt.0.) print *,'warning, volume pre-JL < 0'
	if(volJLCM.lt.0.) print *,'warning, volume post-JL < 0'
	! flatness index = 3*V/A**(3/2)
	flat=3*volFullCM/areabaseFull**(1.5)	!pre-JL
	flatJL=3*volJLCM/areabaseJL**(1.5)	!post-JL
	print *, 'Flatness:         ', real(flat),real(flatJL)
	Write(50,*) 'Flatness:         ', real(flat),real(flatJL)
c
c	! Tetrahedron solid angle 
c	! Edge vectors u,v,w
	call solidangle(kgeno,u,v,w,solid)	! full phase space
	call solidangle(icomb,uJL,vJL,wJL,solidJL)	! JL  subspace
	print *, 'Solid angle:      ',real(solid/pi/2),
     &   real(solidJL/pi/2),'    (ideally: 1)'	
	write(50,*) 'Solid angle:      ',real(solid/pi/2),
     &   real(solidJL/pi/2),'    (ideally: 1)'	
	write(51,*)
	write(51,*) 'MODEL: ',adjustr(pops(1)),' = ',pops(2:4)
	write(50,*)
c	
	CALL CPU_TIME(TIME)
	print *,'Total CPU= ',real(TIME),'s.'
	end
c	!-----------------------------------------------------------------
	subroutine solidangle(n,vxa,vxb,vxc,solid)
c	implicit double precision(a-h,o-z)
	dimension vxa(n),vxb(n),vxc(n)
	real*4 s1,s2,s3
c	! tetrahedron solid angle subtended in radians: 'solid'
c	! vxa, vxb, vxc vectors from vertex to a b c tips.
c
	pi=4.d0*atan(1.d0)
!	 solid angle
!	https://math.stackexchange.com/questions/3342761/octahedral-facet-solid-angle
	xa=dot_product(vxa,vxa)
	xb=dot_product(vxb,vxb)
	xc=dot_product(vxc,vxc)
	xaxb=dot_product(vxa,vxb)	!cos a
	xaxc=dot_product(vxa,vxc)	!cos b
	xbxc=dot_product(vxb,vxc)	!cos c
	cosa=xaxb/sqrt(xa*xb)
	cosb=xaxc/sqrt(xa*xc)
	cosc=xbxc/sqrt(xb*xc)
	anga=acos(cosa)
	angb=acos(cosb)
	angc=acos(cosc)
	sina=sqrt(1-cosa**2)
	sinb=sqrt(1-cosb**2)
	sinc=sqrt(1-cosc**2)
	s1=(cosa-cosb*cosc)/sinb/sinc
	s2=(cosb-cosa*cosc)/sina/sinc
	s3=(cosc-cosa*cosb)/sina/sinb
c
	solid=acos(s1)+acos(s2)+acos(s3)-pi
	return
	end
c	!-----------------------------------------------------------------
	subroutine coeffs(kgeno,x,a,b,c,alpha,beta,gamma)
c	implicit double precision(a-h,o-z)
	dimension x(kgeno),a(kgeno),b(kgeno),c(kgeno)
!	3-way mixture coeffcients and hybrid angle in full phase space
!	x=alpha*a+beta*b+(1-alpha-beta)*c
	ac=dot_product(a-c,a-c)
	bc=dot_product(b-c,b-c)
	xcac=dot_product(x-c,a-c)
	xcbc=dot_product(x-c,b-c)
	acbc=dot_product(a-c,b-c)
	denom=ac*bc-acbc**2
	alpha=(bc*xcac-acbc*xcbc)/denom
	beta=(ac*xcbc-acbc*xcac)/denom
	gamma=1.d0-alpha-beta
	return
	end
!	-----------------------------------------------------------------
	subroutine coeffsJL(kgeno,kpop,x,a,b,c,freqs,
     &   alpha,beta,gamma)
c	implicit double precision(a-h,o-z)
	dimension x(kgeno),a(kgeno),b(kgeno),c(kgeno)
	dimension freqs(kpop,kgeno)
	dimension u(kgeno),v(kgeno)
	dimension ac((kpop-4)*(kpop-5)/2),bc((kpop-4)*(kpop-5)/2)
	dimension xc((kpop-4)*(kpop-5)/2)
c	!3-way mixture coeffcients and hybrid angle in JL subspace
c	! x=alpha*a+beta*b+(1-alpha-beta)*c
c	! kombi: #pairs of aux pops
c	! 4 first columns of freqs() are replaced with x,a,b,c 
c	! and columns 5 to kpop are auxiliary populations: u and v
c
	ac=0.d0
	bc=0.d0
	xc=0.d0
	icomb=1
	do 6660 k=6,kpop
	do 666 j=5,k-1
	u(1:kgeno)=freqs(k,1:kgeno)
	v(1:kgeno)=freqs(j,1:kgeno)
!	renomalization factor
	fab=sqrt(dot_product(u-v,u-v))
! 	scalar products in JL
	ac(icomb)=dot_product(a-c,u-v)/fab
	bc(icomb)=dot_product(b-c,u-v)/fab
	xc(icomb)=dot_product(x-c,u-v)/fab
	icomb=icomb+1
666	enddo		!combinat
6660	enddo		!combinat
	icomb=icomb-1
c
	d1=dot_product(ac,ac)*dot_product(bc,bc)
	d2=(dot_product(ac,bc))**2
	delta=d1-d2
	a1=dot_product(bc,bc)*dot_product(ac,xc)
	a2=dot_product(ac,bc)*dot_product(bc,xc)
	alpha=(a1-a2)/delta
	b1=dot_product(ac,ac)*dot_product(bc,xc)
	b2=dot_product(ac,bc)*dot_product(ac,xc)
	beta=(b1-b2)/delta
	gamma=1.d0-alpha-beta	
	return
	end	
	
	
	
	
	
	
	
	