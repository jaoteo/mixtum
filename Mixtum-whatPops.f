!	reads file 'selectedPopulations.dat', with the XYZ label of the dataset 
!	in the first line and sort 
!	population labels and number of individuals in each of them.
!	Outcome in file: 'XYZ-popula-number.dat'.
!	****************************************************************************************
!	oteo@uv.es, Departament de Fisica Teorica, Universitat de Valencia, Spain.
!	July, 2024
!	****************************************************************************************
c
	parameter(indivz=50000,kpop=10000)	!upper bounds
	character*100, popindiv(indivz),pops(kpop),dicco(kpop)
	character*100 res,dataset,popula,snp
	character*1 resmf
	dimension ndicco(kpop)
	
	open(unit=9, file='selectedPopulations.dat',status='old')	!driver
	read(9,*) dataset
	res=adjustl(adjustr(dataset)//'.ind')
	open(unit=10,file=res)
	popula=adjustl(adjustr(dataset)//'-PopulaNumber.dat')
	open(unit=20,file=popula)
	snp=adjustl(adjustr(dataset)//'.snp')
	open(unit=30,file=snp)
c
	do k=1,indivz
	read(10,*,end=555) res,resmf,popindiv(k)
	enddo
555	indiv=k-1
	type *, 'read done',indiv
	pops(1)=popindiv(1)
	np=1
	do 10 k=2,indiv
	do j=1,np
	if(pops(j).eq.popindiv(k)) goto 10 
	enddo
	np=np+1
	pops(np)=popindiv(k)
10	continue
c
	ndicco=0
	dicco=pops
	do 100 j=1,np
	do k=1,indiv
	if(dicco(j).eq.popindiv(k)) ndicco(j)=ndicco(j)+1
	enddo
	type *, j,ndicco(j)
100	continue
!	 writes population and number of individuals  !!!!!!!!!!!!!!!!!!!!
	do k=1,np
	write(20,888)  dicco(k),ndicco(k)
	enddo
888	format(A100,2x,I4)
c	!SNPs
	do k=1,10**7
	read(30,2020,end=777) res
	enddo
777	k=k-1
	write(20,*) 'Number of SNPs: ', k
2020	format(A30)
	end
	
	
	