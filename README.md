# genderGapCode
Code used to extract and analyse gender data from PubMed and arXiv for a 2018 [PloS Biology](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2004956) paper by Luke Holman, Devi Stuart-Fox, and Cindy Hauser.

See https://lukeholman.github.io/genderGap/ for a web app that lets you explore the data set from that paper, and see https://github.com/lukeholman/genderGap for the Javascript code and .json data used by the web app. 

#### A note about the data and code
The code relies heavily on a SQLite3 database that was created from the PubMed XML files using the script **Processing the Pubmed XML files into useable data.R**, and then analysed in detail in **Plots and analyses.R**. This database is around 2.5GB and so it is not archived on Github. Intead, you can download it from the Open Science Framework at https://osf.io/bt9ya/. My R code retrieves information from the database using the dplyr and dbplyr R packages; here is a [guide](https://cran.r-project.org/web/packages/dbplyr/vignettes/dbplyr.html) to getting started accessing databases using R. 

For people who prefer to work with data in a spreadsheet, I have also made a summary of the gender data split by country, journal, year, and authorship position, called **Author counts for journal country and position.csv.zip**. This file is comprehensive enough to facilitate many of the analyses you might wish to do. However, in order to keep the file size manageable, this spreadsheet discards some of the information that you could get from the database. For example, it discards some information on author order, the precise publication date of every paper, the PMIDs, DOIs and titles of each paper, etc.

If you'd like to get a more up-to-date version of the PubMed data (which is current up to August 2016), I suggest you download a local copy of PubMed using the script **script to download a local copy of Medline.R** (note that it will be over 20GB), and then work through the R script **Processing the Pubmed XML files into useable data.R** to extract useful information from PubMed's XML files. You could also look into the R package ``RISmed``, which can be used to search PubMed from R. This latter option is definitely preferable for smaller studies focused on a particular subset of the literature, but it is not a sensible way to download the entirety of PubMed. For arXiv, I used the easy-to-use R package ``aRxiv`` to get the whole dataset.


### Description of each file

#### Author counts for journal country and position.csv.zip
Zipped .csv file containing a summary of the dataset. The spreadsheet gives the number of male, female and unknown-gender authors that were counted for each combination of year, authorship position (i.e. first/last/middle/single), country (including 'unknown', which refers either to authors with no affiliation, or those with an affiliation for which we could not identify the country), and journal (using the abbreviations favoured by PubMed). 'First' and 'Last' authors were counted from all papers with 2 or more authors. 'Middle' authors are any authors other than the first and last, on papers with three or more authors. Single authors are the authors of papers that list only one author. The unknown-gender authors are people who only gave initials, those whose names were not listed on genderize.io, or those with names that are not associated with one gender >95% of the time (e.g. Alex, Robin).

#### script to download a local copy of Medline.R
This R script writes a Shell script that can be used to download a local copy of the Medline database (which is >20 GB), in order to process the PubMed XML files locally.

##### Data making functions.R
A series of functions used to make the datasets used in the analysis (e.g. to extract useful information from PubMed's messy XML files, and to assign gender to author names)

#### Processing the Pubmed XML files into useable data.R
Script that processes the PubMed XML files, using functions from 'Data making functions.R'.

#### Genderize script.R
Script to call the genderize.io API from R, and get the genders associated with each name (in a country-specific manner, where possible). Note that genderize.io is not free to use.

#### Plot and analysis functions.R
Functions used to make plots or statistically analyse the data from PubMed and arXiv. Called with source() from 'Plots and analyses.R'.

#### Plots and analyses.R
Script used to produce the figures, tables and statistical results in the paper.

#### arXiv mining script.R
Script used to collect gender data from arXiv using R.
