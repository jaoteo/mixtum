# Mixtum
The goal of these Fortran codes is to estimate, from SNP data, the admixture proportion 
in which a 2-way or 3-way population admixture process took place. 
The *f*-formalism [1] is used following the geometric interpretation
in [2,3] to estimate the proportions and other indicators.

Schematically:

1. Two-way: AA=alpha * C1+(1-alpha) * C2, where AA: Hybrid, C1 & C2: Contributors, 
and  alpha: mixing proportion.
AA, C1 and C2 are vectors (arrays) whose components are allele frequencies.

2. Three-way: AA = alpha * C1+beta * C2+(1-alpha-beta) * C3.

Codes are in Fortran. They have been compiled with Intel Fortran (Windows) and 
GNU Fortran (Linux and Windows).
A Python/GUI versión is in progress at
https://github.com/jmcastelo

The dataset input is the *Eigenformat Triad*: XYZ.geno, XYZ,ind, XYZ.snp.
All the 3 files XYZ are ascii. The two later ones can be inspected on the display
whereas the visualization of the former one depends on the text editor capacity. 

XYZ is the label of the Triad and is read by the codes in the first line of the ascii file 
'selectedPopulations.dat'. The next 3 (or 4) lines have the labels of the populations 
AA, C1 and C2 (and eventually C3), in the 2-way (eventually 3-way) case.
All the subsequent lines are the labels of the Auxiliary populations to be used, 
one per line. This driver file is provided by the user. 

There are 4 Fortran codes:

1. **MixtuM-whatPops.f**: reads the first line in 'selectedPopulations.dat' to get the label XYZ. 
Then, it reads XYZ.ind and tells what are the populations inside, 
as well as the number of individuals in each one, and also the number of SNPs in XYZ.snp.
Outcomes are in XYZ-popula-number.dat.

2. **MixtuM-allele-frequencies.f**: reads the first line in 'selectedPopulations.dat' to get XYZ. 
Next, it reads the label of the populations, one in each line.
The outcome of the code is in the file  XYZ-allele-frequencies.dat, with the names of each 
population heading the allele frequencies in columns (one per population). Every line is a SNP. 
The number of lines may be smaller than the nominal number of SNPs in the Triad. 
This is because whenever all the individuals in a particular population have an 
entry '9' (=N/A) in XYZ.geno, then that SNP is supressed for all the
populations enlisted in selectedPopulations.dat. 
The effective copy XYZ.snp file is given in XYZ-effective-SNP.dat.
The ordering of columns in XYZ-allele-frequencies.dat is linked to its later use.
It is the same order as in the lines of selectedPopulations.dat. 
The first population (first column) is intended to be the Hybrid in the 2-way (or 3-way) 
model under escrutiny. The second and third populations (columns) correspond to 
Parental-1 and Parental-2, respectively. In a 3-way model, the fourth column is interpreted as Parental-3.
The remaining populations (columns) are allele frequencies of the Auxiliary Populations. 

3. **MixtuM-2way.f**: reads selectedPopulations.dat to ascertain the dataset label XYZ and to count 
the total number of populations involved (columns in XYZ.geno). Then it computes f-statistics for a 2-way model.
The outcomes are in the following three files:

3A. **Mixtum-f4-points.dat**: f4's values (renormalised and not-renormalised) in columns (see the header), 
f4-ratios and the corresponding Auxiliary population pair. 
This file is intended to plot the outcomes.  A summary of results is at bottom:
admixture model,SNPs, Auxiliary pairs, proportions (f4-ratio mean, pre-JL and post-JL, 
renormalised and not-renormalised), pre-JL and post-JL angle, f3 admixture test.

3B. **Mixtum-out-Brief.dat**: Summary of results. It is incremental (all runs are added).

3C. **Mixtum-f4-alphas-OK.dat**: alpha estimates via the f4-ratio that fall in the range [0,1].

4. **MixtuM-3way**: reads selectedPopulations.dat to ascertain the dataset label XYZ and to count 
the total number of populations involved (columns in XYZ.geno). Then it computes f-statistics for a 3-way model.
The outcomes are in the '*.dat' files.

CAVEAT: These are Fortran codes and therefore need to be compiled to get an executable.
Command line: 

1.     gfortran -o MixtuM-* MixtuM-*.f      (GNU Fortran, Windows & Linux. It gives MixtuM-*.exe)
2.     ifort MixtuM-*.f                     (Intel Fortran, Windows. It gives Mixtum-*.exe)
         
Source code files contain further details/comments about the codes. 

oteo@uv.es, 
JA Oteo, Departament de Fisica Teorica, Universitat de Valencia, Spain.
July, 2024

[1]: Patterson N, Moorjani P, Luo Y, Mallick S, Rohland N, Zhan Y,
Genschoreck T, Webster T, Reich D. 2012. 'Ancient admixture
in human history'. GENETICS. 192:1065–93

[2]: Oteo JA, Oteo-García G. 2024. 'The geometry of admixture in population genetics:
the blessing of dimensionality'. GENETICS (to be published)

[3]: Oteo-García G, Oteo JA. 2021. 'A Geometrical Framework for
f-Statistics'. Bulletin of Mathematical Biology. 83:14
