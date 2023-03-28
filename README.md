# PedigreeToolsF90
## Author: John Woolliams 
## Version: 

PedigreeTools.f90 is pedigree cleaning protocol made up of several programmes. Each of the programmes does some pedigree cleaning and processing and produces input for the next stage of cleaning (based on the input pedigree). This readme file can be used together with diagrams.pptx for better understanding. PedigreeTools.f90 consists of:

    pedv7.f90    – analyse and renumbers pedigree
    kosaraju.f90 – finds all circular subsets
    johnson.f90  – lists all distinct circuits

When the pedigree is cleaned two programmes can be used for: 

    wmlv5.f90   – calculating pedigree F
    colleau.f90 – calculating average relationships

### Input files:
  1.	Mating list with columns: id, sire, dam, sex (1 male, 2 female)). 
      In our toy examples it’s “testv1.csv” (and testv2.csv, testv3.csv and testv4.csv). 
      Both 0 and NA is recogniseda as missing value.

```
c	0	b	2
a	0	0	1
b	A	0	0
c	D	b	0
e	d	b	1
f	a	e	0
g	d	c	0
h	h	g	0
k	0	0	2
m	0	k	1
n	0	k	0
p	m	n	0
0	0	0	0
NA	NA	NA	0
```
  2.	pedv7_in.txt file (input list and options)
```
testv1.csv          ! (unordered) pedigree mating list, input file
0                   ! header lines to skip in input file
testv1              ! tag for output files
1                   ! missing birth records
1                   ! selfings
1                   ! mixed use parents
1                   ! multiple record incompatibilities
1                   ! self parenting
```

### Programmes

**pedv7.f90**

Takes in information from pedigree and mating list and does analysis and renumbering.
Depending on the input pedigree and what the program will find - it gives output files in your folder:

```
pedv7_log.txt  – summary of run
pedv7_err1.csv – added parents
pedv7_err2.csv – selfing alert
pedv7_err3.csv – mixed use parents
pedv7_err4.csv – duplicate errors
pedv7_err5.csv – self-parentiing
pedv7_info.csv – summary id info (mating list information)
```

If the pedigree is considered valid – it will give output files (which are input files for next stage – calculating pedigree F):
```
pedv7_ped.csv
wmlv5_in.txt
```
   
Using **wmlv5.f90** gives output files:
```
wmlv5_log.txt  – summary of run
wmlv5_pedf.csv – updated pedigree
colle_in.txt   – next input file
```

**colleau.f90** calculates average relationships of an individual with all in "active" set eg. the current breeding population.
```
colle_actv.txt – active set list
wmlv5_pedf.csv – pedigree file
colle_in.txt   – input file
```
The active set list is manually made list that one should made it based on it's own preferences. From testv1 can look like this :

```
e
f
g
```

and colleau.f90 then writes:
```
colle_log.txt  – summary of run
colle_pedc.csv – outputs 
```


## In the pedigree is NOT considered valid and has serious errors it will give next file to continue cleaning:
```
kosaraju_in.txt
```
**kosaraju.f90** finds all circular subsets and takes in input files:
```
kosaraju_in.txt – input file
pedv7_info.csv  – summary id info
```
and writes:
```
kosar_log.txt  – summary of run
kosar_subg.csv – minimal subgraph
kosar_sets.txt – discrete sets
johnson_in.txt – next input file
```
**johnson.f90** finds all circular subsets using input files:
```
johnson_in.txt
kosar_sets.csv
```
and writes output files:
```
johns_log.txt  – summary of run
johns_subg.csv – strong subgraphs
john_circ.csv  – circular paths
johns_freq.csv – frequency in paths
```

After this step one should consider manual cleaning and checking the raw pedigree. 
