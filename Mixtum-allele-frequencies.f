!	Computes Allele Frequencies of a set of populations in the dataset XYZ,
!	given by the three files in eigenformat: XYZ.geno  & XYZ.ind  & XYZ.snp.
!	The label XYZ (less that 100 charcaters) is given in the first line  
!	of the driver file 'selectedPopulations.dat'.
!	Each succesive line is a population label of the selected populations.
!	kpopmax: max number of selected populations allowed. Default: kpopmax=400 (may be changed).
!	max1pop: max number of individual in one population. Default: max1pop=1000 (may be changed).
!	The first label read is intended to be the hybrid of an admixture model whereas the 2dn & the 3rd 
!	labels read are the parental-1 and the parental-2.  The remaining labels read below these three ones
!	are considered as the Auxiliary Populations of the f-statistics formalism.
!	Output: 'XYZ-Allele-Frequencies.dat', has the allele freqs 
!	in kpop columns (populations) and kgeno rows (SNPs). See 'write(30,*)', The
!	first line contains the population labels (text). 
!	Output: 'XYZ-AlleleFreqsSNP.dat', has the labels
!	of the effective SNPs kept in the allele freqs final dataset.
!	Caveat: input data in XYZ.geno are (0,1,2 or 9). No separator, e.g., 022110091....
!	Values: 9 = N/A; 0,1,2 = number of mutations.
!	Caveat: data in 'XYZ.ind' are of type 'character', in 3 columns &
!	some special characters may not be read or may fool the reading.
!	Only the 3rd column of 'XYZ.ind' is used.
!	indiv: number of individuals = number of columns in datafile.
!	indivmax: max number of individuals (columns). May be changed.
!	kpop: total number of populations selected. Default k=40. May be changed.
!	freqs(*): allele frequencies as a function of population.
!	pops(*): population labels.
!	kgeno: max number of SNPs to be read. May be increased if allowed in compilation.
!	kgenoeff: number of effective SNPs. Same number of SNPs for all populations. 
!	A SNP is dropped from the dataset whenenever '9' in every individual of one population.
! 	Removing 'implicit double precision (a-h,o-z)' allows for higher dimension arrays in compilation.
!	Compile (Intel, Windows): ...\workdir> ifort make-allele-frequencies-Mixtum.f    
!	Compile (GNU, Windows): ...\workdir> gfortran -o make-allele-frequencies-Mixtum make-allele-frequencies-Mixtum.f
!	Then execute: make-allele-frequencies-Mixtum
!	CPU time depends on the number of populations and SNPs, and
!	on the location of the rightmost individual to be read in *.geno. The further to the right, the longer the CPU time.
!	****************************************************************************************
!	oteo@uv.es, Departament de Fisica Teorica, Universitat de Valencia, Spain.
!	July, 2024
!	****************************************************************************************
c
	implicit double precision (a-h,o-z)
c	parameter(kgeno=597568,indivmax=7744,kpopmax=1320)	! TOP values for v37.2.1240K_HumanOrigins
	parameter(kgeno=600000,indivmax=10000,kpopmax=400)
c	parameter(kgeno=500,indivmax=10000,kpopmax=400)	!for testing, trials
	parameter(max1pop=1000)	! max number of individuals in one population
	dimension kdicco(kpopmax,max1pop), ndicco(max1pop)
	dimension  freqs(kpopmax), ng(indivmax)
	character*100  pops(kpopmax)
	character*30, res,popula(indivmax)
	character*1 resmf,base1,base2
	character*20 snplabel
	character*100 dataset
	character*150 labelind,labelgeno,labelsnp,outfreqs,outsnps
c
	open(unit=25,file='selectedPopulations.dat',status='old') !input
	read(25,*) dataset
	print *
	print *, 'dataset: ', adjustl(adjustr(dataset))
c	dataset='v37.2.1240K_HumanOrigins_unzip'
	labelind=adjustl(adjustr(dataset)//'.ind')
	labelgeno=adjustl(adjustr(dataset)//'.geno')
	labelsnp=adjustl(adjustr(dataset)//'.snp')
	outfreqs=
     &   adjustl(adjustr(dataset)//'-Allele-Frequencies.dat')
	outsnps=adjustl(adjustr(dataset)//'-effective-SNP.dat')
	open(unit=10,file=labelind)	!input
	open(unit=20,file=labelgeno) 	!input
	open(unit=21,file=labelsnp) 	!input
	open(unit=30,file=outfreqs)	!output
	open(unit=31,file=outsnps) 	!output
c
	do k=1,kpopmax
	read(25,*,end=1001) pops(k)
	enddo
1001	kpop=k-1
	do k=1,indivmax
	read(10,*,end=3333) res,resmf,popula(k)
	enddo
3333	indiv=k-1
	do 456 kpopula=1,kpop
	knum=0
	do 123 kindiv=1,indiv
	if(popula(kindiv).eq.pops(kpopula)) then
	knum=knum+1
	kdicco(kpopula,knum)=kindiv !line in *.ind and column in *.geno, where the individual is located
	endif
	ndicco(kpopula)=knum	!# individuals in population 'kpopula'
123	continue  	!kindiv
456	continue	!kpopula
	print *, 'Number of populations selected:', kpop
	print *, 'Number of individuals in dataset:', indiv
	print *, 'Selected populations and number of individuals:'
	maxcol=maxval(kdicco)	! rightmost column to be read in *.geno
	print 5551, maxcol
5551	format(' Rightmost column to be read in *.geno:',I6)
	print *
	do i=1,kpop
	print 222, pops(i), ndicco(i)
	enddo
222	format(A20,I4)
	write(30,555) adjustl(pops(1:kpop))	!header in allele freqs dataset
555	format(3x,100(A14,1x))
c
	do ksnps=1,kgeno
	read(21,*,end=5552)
	enddo
5552	maxgeno=ksnps-1
	rewind(21)
	print *,'max number of SNPs to be read:',maxgeno
!	 Computing allele frequencies
	print *, 'Computing allele frequencies. Wait ...'
	kgenoeff=0
	eps=1.d-8
	top=1.d0+eps
	bot=-eps
	do 654 ksnp=1,kgeno	!SNP
	read(20,999,end=7777) ng(1:indiv)	!reads SNP values for all populations
999	format(1000000I1)
	read(21,*) snplabel,msnp,rsnp,jsnp,base1,base2
	if(mod(ksnp,10000).eq.0) print *, 'SNPs processed:',ksnp 	!do not dismay
	do 987 kpopula=1,kpop
	fq=0.d0
	npopeff=ndicco(kpopula)	!effective number of individuals in that population because N/A's(=9)
	indivpop=0
c	do 1 i=1,indiv
	do 1 i=1,maxcol
	do 2 j=1,ndicco(kpopula)
	if(i.eq.kdicco(kpopula,j)) then	    !choose column in *.geno
	if(ng(i).eq.9) then	!discard that individual because N/A
	npopeff=npopeff-1
	if(npopeff.le.0) then	!discard that SNP in the allele freqs dataset
	print *, 'SNP invalid:', ksnp , i,ng(i)
	goto 654
	endif
	goto 2
	endif
	fq=fq+(2-ng(i))/2.d0
	indivpop=indivpop+1
	endif
2	continue
1	continue
	freqs(kpopula)=fq/npopeff
987	continue
	write(30,888) freqs(1:kpop)
	kgenoeff=kgenoeff+1
	do m=1,kpop	!Checks  abs(allele freq)<1
	if(freqs(m).gt.top) print *, 'warning, freq>1',ksnp,m,freqs(m)
	if(freqs(m).lt.bot) print *, 'warning, freq<0',ksnp,m,freqs(m)
	enddo	
	write(31,1551) adjustr(snplabel),msnp,rsnp,jsnp,base1,base2 
654	continue
1551	format(A20,I8,f15.6,I20,2A2)
888	format(10000e15.6)
7777	print *, 'Effective number of SNPs:',kgenoeff
c	
	CALL CPU_TIME(TIME)
	print *,'Total CPU= ',real(TIME),'s.'
	end