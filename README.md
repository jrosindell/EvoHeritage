# EvoHeritage
A set of software tools for EvoHeritage calculations in R

# EvoHeritage tools
The root directory contains a set of functions for calculating EvoHeritage that can be sourced and then used in analyses
The root directory also contains a test file that can be sourced to confirm these functions are working correctly
The working directory for running the unit tests should be the root

# Example applications
Sub directories contain individual example applications of the main functions
Each example application has a directory called Data which is excluded (gitignored) from the GitHub repository
The working directory for running any of the example applications should be the directory of that application

## Living fossils application

Calculated living fossilness for mammals using various methods based on ED and EvoHeritage
Graphs the results and provides a comparison between the methods

## Primary productivity application

Compares species richness, PD and EvoHeritage as descriptors of primary productivity in the Cedar Creek grassland experiments
graphs the results comparing the methods

# Code authorship and citation
The EvoHeritage tools and living fossil applications are written by James Rosindell
The primary productivity application was written by James Rosindell and William D Pearse
These tools and the two applications are as described in "Phylogenetic Biodiversity Metrics Should Account for Both Accumulation and Attrition of Evolutionary Heritage" 2023 authored by James Rosindell, Kerry Manson, Rikki Gumbs, William D Pearse and Mike Steel
https://doi.org/10.1101/2022.07.16.499419

# Contact details
For any questions about using this code please write to j.rosindell@imperial.ac.uk