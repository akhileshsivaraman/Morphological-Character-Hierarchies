# Morphological Character Hierarchies
Morphological characters often exist in hierarchies whereby the state of one character controls the state of a subset of characters. This violates the assumptions of phylogenetic analysis that characters are independent and applicable to all taxa. As these violations are known to have the possibility of misleading phylogenetic analyses, two programmes, Morphy and Anagallis, have been developed to address the issue of inapplicability by accounting for the violations. The phylogenies produced by Morphy and Anagallis were compared to those of PAUP*, a traditional programme of phylogenetic analysis, to investigate the merit of the new programmes.

This analysis was undertaken as my final year project at Imperial College London under the supervision of Dr. Martin Brazeau with additional inputs from Thomas Guillerme, at the University of Sheffield.

In this repository, you can find all the code and data used for the project. As I'm still very new to GitHub, this repository is a bit of mess at the moment! But, I hope that this brief description will make it easier to navigate the repository.

The files present in the "Data" folder are named after the lead author(s) whose morphological data sets I analysed with PAUP*, Morphy and Anagallis.

Within each of the study folders, you can find character matrices, character lists, files containing the trees found by each programme and the R scripts used to analyse those trees.
Files relating to PAUP* are marked "Paup", Morphy marked "Morphy" and Anagallis marked "Ana".
Character matrices can be identified by the extensions .nex and .tnt while tree lists contain "tree" in their filenames (files labelled "contree" describe the consensus tree generated by one of the programmes).

Data obtained from tree comparison analyses are pooled into .csv (and .numbers) files in the "Data" folder and the R scripts in this folder were used to carry out analyses across the data sets.
