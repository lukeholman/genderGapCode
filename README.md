# genderGapCode
Code used to extract and analyse gender data from PubMed and arXiv for a forthcoming paper by Luke Holman, Devi Stuart-Fox, and Cindy Hauser.

### Description of each file

#### Author counts for journal country and position.csv.zip
Zipped .csv file containing a summary of the dataset. The spreadsheet gives the number of male, female and unknown-gender authors that were counted for each combination of year, authorship position (i.e. first/last/middle/single), country (including 'unknown', which refers either to authors with no affiliation, or those with an affiliation for which we could not identify the country), and journal (using the abbreviations favoured by PubMed). 'First' and 'Last' authors were counted from all papers with 2 or more authors. 'Middle' authors are any authors other than the first and last, on papers with three or more authors. Single authors are the authors of papers that list only one author.

##### Data making functions.R
A series of functions used to make the datasets used in the analysis (e.g. to extract useful information from PubMed's XML data, and to assign gender to author names)

#### Genderize script.R
Script to call the genderize.io API from R, and get the genders associated with each name (in a country-specific manner, where possible)

#### journal_disciplines.csv
A csv file that gives the discipline that was assigned to each journal. Use this if you want to interpret the data in 'Author counts for journal country and position.csv.zip' in terms of disciplines. 

#### Plot and analysis functions.R
Functions used to make plots or statistically analyse the data from PubMed and arXiv. Called with source() from 'Plots and analyses.R'.

#### Plots and analyses.R
Script used to produce all the figures, tables and statistical results in the paper.

#### Processing the Pubmed XML files into useable data.R
Script that processes the PubMed XML files, using functions from 'Data making functions.R'.

#### arXiv mining script.R
Script used to collect gender data from arXiv.

#### script to download a local copy of Medline.R
This R script writes a Shell script that can be used to download a local copy of the Medline database (which is >20 GB), in order to process the PubMed XML files locally.





