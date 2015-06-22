# GenomeResearch
This repository includes the codes used for the in silico management of populations in Bosse et al (Genome Research, 2015). 

It includes the data required to run the code as well as output. 

Codes included: CebifronsNeutGeneal.f90 CebifronsNeutMol.f90 CebifronsNeutSeg.f90 CebifronsNonNeutGeneal.f90 CebifronsNonNeutMol.f90 CebifronsNonNeutSeg.f90

which were developed to analyse the filtered resequencing data from Sus Cebifrons. 

They all perform first an expansion to 10 individuals and then management using optimal contributions. 

The differences between the codes are whether OC are based on genealogies (..GeneaL.f90), on molecular IBS coancestry (...Mol.f90) or on shared segments coancestry (...Seg.f90) Neutral scenarios with the first four codes, and assuming some deleterious load in the last three.

Codes developed to analyse the genotype data from the Pietrain breed included are: PietrainNeutGeneal.f90 PietrainNeutMol.f90 PietrainNeutSeg.f90 which perform OC with genealogies, with IBS or with shared segments.
