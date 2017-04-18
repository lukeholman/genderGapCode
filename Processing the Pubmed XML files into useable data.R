# First set the working directory to the location of this script
setwd("~/Desktop/Gender Pubmed mining project/R scripts")

# Load libraries needed - install them first if you need to using install.packages()
library(XML)
suppressPackageStartupMessages(library(dplyr))
library(reshape2)
library(stringr)
source("Data making functions.R")


######################################################################################################
# PARSING A LOCAL COPY OF THE WHOLE OF PUBMED
######################################################################################################

# Parse all the XML files from the "Update" part of Pubmed into more readable csv files, containing just the info that we want.
# The result includes all the citatins, even those with no authors etc, which we will subsequently cull. 
# Note that you can run each of these in its own R window to greatly save time - I have 8 cores on my Mac, so I split the job into 8 chunks. Still took about a week.
parse.many.xml.files("../data/Update files - raw xml", "../data/Update files - parsed", focal.chunk = 1)
parse.many.xml.files("../data/Main files - raw xml", "../data/Main files - parsed", focal.chunk = 1)

parse.many.xml.files("../data/Update files - raw xml", "../data/Update files - parsed", focal.chunk = 2)
parse.many.xml.files("../data/Main files - raw xml", "../data/Main files - parsed", focal.chunk = 2)

parse.many.xml.files("../data/Update files - raw xml", "../data/Update files - parsed", focal.chunk = 3)
parse.many.xml.files("../data/Main files - raw xml", "../data/Main files - parsed", focal.chunk = 3)

parse.many.xml.files("../data/Update files - raw xml", "../data/Update files - parsed", focal.chunk = 4)
parse.many.xml.files("../data/Main files - raw xml", "../data/Main files - parsed", focal.chunk = 4)

parse.many.xml.files("../data/Update files - raw xml", "../data/Update files - parsed", focal.chunk = 5)
parse.many.xml.files("../data/Main files - raw xml", "../data/Main files - parsed", focal.chunk = 5)

parse.many.xml.files("../data/Update files - raw xml", "../data/Update files - parsed", focal.chunk = 6)
parse.many.xml.files("../data/Main files - raw xml", "../data/Main files - parsed", focal.chunk = 6)

parse.many.xml.files("../data/Update files - raw xml", "../data/Update files - parsed", focal.chunk = 7)
parse.many.xml.files("../data/Main files - raw xml", "../data/Main files - parsed", focal.chunk = 7)

parse.many.xml.files("../data/Update files - raw xml", "../data/Update files - parsed", focal.chunk = 8)
parse.many.xml.files("../data/Main files - raw xml", "../data/Main files - parsed", focal.chunk = 8)


# Zip all the parsed files into a single massive object, and write it to disk
pubmed <- zip.all.files.together()
# write.csv(pubmed, file = "../data/All data together.csv", row.names = F)

pubmed <- pubmed[!duplicated(pubmed$pmid), ] # Few duplicates to remove

# Throw out all journals for which we have < 100 papers fitting this description
pubmed <- pubmed[pubmed$journal %in% names(table(pubmed$journal)[table(pubmed$journal) > 100]), ]

# write.csv(pubmed, file = "../outputs/Pubmed data for analysis.csv", row.names = F)


######################################################################################################
# GETTING A RESEARCH DISCIPLINE FOR EACH OF THE JOURNALS 
######################################################################################################

# Start off by parsing the XML file provided by PubMed with all the info they have about each journal they index
journal.metadata <- parse.journal.descriptions.xml.file("../data/lsi2016.xml")
journal.metadata$n <- melt(table(pubmed$journal))$value[match(journal.metadata$short.title, melt(table(pubmed$journal))$Var1)]
#write.csv(journal.metadata, file="../data/journal.metadata.csv", row.names = F)

# Restrict the metadata to just the journals that we have in our focal dataset
journal.metadata <- read.csv("../data/journal.metadata.csv", stringsAsFactors = F)
culled.metadata <- journal.metadata[!is.na(journal.metadata$n), ]

# Many of the journals have not been assigned a "broad.journal.headings" entry by Pubmed.
# However, Pubmed have performed text mining and assigned "MeSH" terms to many journals - the MeSH terms are often the same as one of the broad.journal.headings categories.
# Therefore, let's check the MeSH column for a value that can be used to fill in the broad.journal.headings column:

# Make an approved list of topics, which is the unique values present in the "broad.journal.headings" column
topics <- sort(unique(culled.metadata$broad.journal.headings))

# If there is one term which matches the list, fill it in the broad.journal.headings column.
mesh.for.missing.ones <- strsplit(culled.metadata[is.na(culled.metadata$broad.journal.headings), ]$mesh.headings, split="_")
hits <- lapply(mesh.for.missing.ones, function(x) x[x %in% topics])
culled.metadata$broad.journal.headings[which(is.na(culled.metadata$broad.journal.headings))[sapply(hits, length) == 1]] <- unlist(hits[sapply(hits, length) == 1])

# If there is more than one matching MeSH term, pick one at random, and fill it in the broad.journal.headings column.
mesh.for.missing.ones <- strsplit(culled.metadata[is.na(culled.metadata$broad.journal.headings), ]$mesh.headings, split="_")
hits <- lapply(mesh.for.missing.ones, function(x) x[x %in% topics])
culled.metadata$broad.journal.headings[which(is.na(culled.metadata$broad.journal.headings))[sapply(hits, length) > 1]] <- sapply(hits[sapply(hits, length) > 1], function(x) x[sample(length(x))][1])

# There are still many journals with no category, and this includes many major journals (e.g. "Nature Medicine" is not classified as anything by PubMed, when it should clearly be "Medicine")
missing.ones <- with(culled.metadata, data.frame(journal = short.title, mesh = mesh.headings, row = 1:nrow(culled.metadata), n=n, stringsAsFactors = F))
missing.ones <- missing.ones[is.na(culled.metadata$broad.journal.headings), ]

# Manually fill some of the broad.journal.headings by using grep on the MeSH terms...
grep.on.Mesh <- function(grep.term, category) {
  if(!(category %in% topics)) print("check for typo - did you mean to create a new category?")
  rows <- 1:nrow(missing.ones)
  missing.ones$broad.journal.headings[rows %in% grep(grep.term, missing.ones$mesh) & is.na(missing.ones$broad.journal.headings)] <- category # Add the category, provided it is currently listed as NA (thus, the earlier replacements take precedence)
  missing.ones
}

# ...or by using grep on the title of the journal
grep.on.title <- function(grep.term, category) {
  if(!(category %in% topics)) print("check for typo - did you mean to create a new category?")
  rows <- 1:nrow(missing.ones)
  missing.ones$broad.journal.headings[rows %in% grep(grep.term, missing.ones$journal) & is.na(missing.ones$broad.journal.headings)] <- category # Add the category, provided it is currently listed as NA (thus, the earlier replacements take precedence)
  missing.ones
}

missing.ones$broad.journal.headings <- NA

# First let's do the classifications that we are very sure should take priority. First we'll do some individual journals that proved difficult to classify using other means
missing.ones <- grep.on.title("Pacing Clin Electrophysiol", "Cardiology")
missing.ones <- grep.on.title("Zhonghua Yi Xue Za Zhi", "Complementary Therapies")
missing.ones <- grep.on.title("Sensors", "Biomedical Engineering")
missing.ones <- grep.on.title("J Acoust Soc Am", "Biophysics")
missing.ones <- grep.on.title("Mater Sci Eng C Mater Biol Appl", "Biomedical Engineering")
missing.ones <- grep.on.title("Med Eng", "Biomedical Engineering")
missing.ones <- grep.on.title("Cochrane Database Syst Rev", "Medicine")
missing.ones <- grep.on.title("Bioresour Technol", "Chemistry")
missing.ones <- grep.on.title("Biomed Res Int", "Medicine")
missing.ones <- grep.on.title("J Vis Exp", "Science")
missing.ones <- grep.on.title("Trials", "Medicine")
missing.ones <- grep.on.title("Carbohydr Res", "Chemistry")
missing.ones <- grep.on.title("BMC Res Notes", "Science")
missing.ones <- grep.on.title("[Mm]acromol", "Biochemistry")
missing.ones <- grep.on.title("ScientificWorldJournal", "Science")
missing.ones <- grep.on.title("Anim Reprod Sci", "Veterinary Medicine")
missing.ones <- grep.on.title("Phys Rev E Stat Nonlin Soft Matter Phys", "Biophysics")
missing.ones <- grep.on.title("ACS Appl Mater Interfaces", "Biochemistry")
missing.ones <- grep.on.title("IEEE Trans Biomed Circuits Syst", "Biomedical Engineering")
missing.ones <- grep.on.title("Nat Commun", "Science")
missing.ones <- grep.on.title("Cell Rep", "Cell Biology")
missing.ones <- grep.on.title("Nat Prod Commun", "Chemistry")
missing.ones <- grep.on.title("Sci Signal", "Cell Biology")
missing.ones <- grep.on.title("Biophys", "Biophysics")
missing.ones <- grep.on.title("Opt Express", "Biomedical Engineering")
missing.ones <- grep.on.title("Small", "Biomedical Engineering")
missing.ones <- grep.on.title("Europace", "Cardiology")
missing.ones <- grep.on.title("J Cereb Blood Flow Metab", "Brain")
missing.ones <- grep.on.title("J R Soc Interface", "Science")
missing.ones <- grep.on.title("Shock", "Medicine")
missing.ones <- grep.on.title("Coll Antropol", "Anthropology")
missing.ones <- grep.on.title("J AAPOS", "Ophthalmology")
missing.ones <- grep.on.title("Infect Genet Evol", "Epidemiology")
missing.ones <- grep.on.title("Mol Phylogenet Evol", "Genetics")
missing.ones <- grep.on.title("Vector Borne Zoonotic Dis", "Tropical Medicine")
missing.ones <- grep.on.title("CNS Drugs", "Psychopharmacology")
missing.ones <- grep.on.title("Philos Trans A Math Phys Eng Sci", "Science")
missing.ones <- grep.on.title("J Asian Nat Prod Res", "Chemistry")
missing.ones <- grep.on.title("Comput Biol", "Computational Biology")
missing.ones <- grep.on.title("Ethn Dis", "Medicine")
missing.ones <- grep.on.title("Zhongguo Gu Shang", "Orthopedics")
missing.ones <- grep.on.title("Biol Trace Elem Res", "Biochemistry")
missing.ones <- grep.on.title("Mol Vis", "Cell Biology")

# now let's assign some categories based on specific health problems, etc
missing.ones <- grep.on.title("AIDS", "Acquired Immunodeficiency Syndrome")
missing.ones <- grep.on.title("HIV", "Acquired Immunodeficiency Syndrome")
missing.ones <- grep.on.title("Stroke", "Neurology")
missing.ones <- grep.on.title("Dermat", "Dermatology")
missing.ones <- grep.on.title("[Oo]nco", "Neoplasms")
missing.ones <- grep.on.title("Leuk", "Neoplasms")
missing.ones <- grep.on.title("Lymphoma", "Neoplasms")
missing.ones <- grep.on.title("Anaesth", "Anesthesiology")
missing.ones <- grep.on.title("Thromb Haemost", "Vascular Diseases")
missing.ones <- grep.on.title("Bone", "Osteopathic Medicine")
missing.ones <- grep.on.title("[Gg]yn[ae]", "Gynecology")
missing.ones <- grep.on.title("BMC Musculoskelet Disord", "Physical and Rehabilitation Medicine")
missing.ones <- grep.on.title("Asthma", "Pulmonary Medicine")
missing.ones <- grep.on.title("Diabetes", "Endocrinology")
missing.ones <- grep.on.title("Headache", "Neurology")
missing.ones <- grep.on.title("Chirality", "Chemistry")
missing.ones <- grep.on.title("Virus Genes", "Virology")
missing.ones <- grep.on.title("Health Serv", "Health Services Research")

# Now do some where the journal's title implies a specific specialty 
missing.ones <- grep.on.title("Gastroenterol", "Gastroenterology")
missing.ones <- grep.on.title("Pancrea", "Gastroenterology")
missing.ones <- grep.on.title("Endosc", "Gastroenterology")
missing.ones <- grep.on.title("Nephrol", "Nephrology")
missing.ones <- grep.on.title("Cerebrovasc", "Brain")
missing.ones <- grep.on.title("Urol", "Urology")
missing.ones <- grep.on.title("Periodontol", "Dentistry")
missing.ones <- grep.on.title("Otorhinolar", "Otolaryngology")
missing.ones <- grep.on.title("Otol", "Otolaryngology")
missing.ones <- grep.on.title("[Cc]ardio", "Cardiology")
missing.ones <- grep.on.title("Circ", "Cardiology")
missing.ones <- grep.on.title("Coron", "Cardiology")
missing.ones <- grep.on.title("Heart", "Cardiology")
missing.ones <- grep.on.title("Dent", "Dentistry")
missing.ones <- grep.on.title("Orthod", "Dentistry")
missing.ones <- grep.on.title("Endocrinol", "Endocrinology")
missing.ones <- grep.on.title("Public Health", "Public Health")
missing.ones <- grep.on.title("[Cc]ancer", "Neoplasms")
missing.ones <- grep.on.title("Head Neck", "Medicine")
missing.ones <- grep.on.title("Transpl", "Transplantation")
missing.ones <- grep.on.title("Vet ", "Veterinary Medicine")
missing.ones <- grep.on.title("Brain", "Brain")
missing.ones <- grep.on.title("Crystallogr", "Chemistry")
missing.ones <- grep.on.title("Chromatogr", "Chemistry")
missing.ones <- grep.on.title("Nurs", "Nursing")
missing.ones <- grep.on.title("Ophthalm", "Ophthalmology")
missing.ones <- grep.on.title("Alzheimer", "Neurology")
missing.ones <- grep.on.title("Toxicol", "Toxicology")
missing.ones <- grep.on.title("[Ii]mmun", "Allergy and Immunology")
missing.ones <- grep.on.title("Allergy", "Allergy and Immunology")
missing.ones <- grep.on.title("[Pp]harm", "Pharmacology")
missing.ones <- grep.on.title("Neurophys", "Neurology")
missing.ones <- grep.on.title("Neurochem", "Psychopharmacology")
missing.ones <- grep.on.title("Appetite", "Nutritional Sciences")
missing.ones <- grep.on.title("Obesity", "Nutritional Sciences")
missing.ones <- grep.on.title("Nutr ", "Nutritional Sciences")
missing.ones <- grep.on.title("Health Serv", "Medicine")
missing.ones <- grep.on.title("Neurobiol", "Neurology")
missing.ones <- grep.on.title("Structure", "Biochemistry")
missing.ones <- grep.on.title("Microb", "Microbiology")
missing.ones <- grep.on.title("Urban Health", "Public Health")
missing.ones <- grep.on.title("Physiol", "Physiology")
missing.ones <- grep.on.title("Parasit Vectors", "Tropical Medicine")
missing.ones <- grep.on.title("[Rr]adiol", "Radiology")
missing.ones <- grep.on.title("[Rr]adiother", "Radiotherapy")
missing.ones <- grep.on.title("Antivir", "Anti-Infective Agents")
missing.ones <- grep.on.title("Antimicrob", "Anti-Infective Agents")
missing.ones <- grep.on.title("Antiinfect", "Anti-Infective Agents")
missing.ones <- grep.on.title("Antibiot", "Anti-Infective Agents")
missing.ones <- grep.on.title("Anti Infect", "Anti-Infective Agents")
missing.ones <- grep.on.title("Anticancer", "Neoplasms")
missing.ones <- grep.on.title("Antibod", "Allergy and Immunology")
missing.ones <- grep.on.title("Androl", "Medicine")
missing.ones <- grep.on.title("Traumatol", "Medicine")
missing.ones <- grep.on.title("Psychiatry", "Psychiatry")
missing.ones <- grep.on.title("Surg", "General Surgery") # note General surgery comes last, so for example heart surgery is under cardiology

# Next let's look in the MeSH categories, which generally have very accurate info about the journal's content
missing.ones <- grep.on.title("Blood", "Hematology")
missing.ones <- grep.on.Mesh("Microbiology", "Microbiology")
missing.ones <- grep.on.Mesh("Pollution", "Toxicology")
missing.ones <- grep.on.Mesh("Stem Cells", "Cell Biology")
missing.ones <- grep.on.Mesh("Cell ", "Cell Biology")
missing.ones <- grep.on.Mesh("Cellular", "Cell Biology")
missing.ones <- grep.on.Mesh("Surgery", "General Surgery")
missing.ones <- grep.on.Mesh("Surgical", "General Surgery")
missing.ones <- grep.on.Mesh("Veterin", "Veterinary Medicine")
missing.ones <- grep.on.Mesh("Psychophysiol", "Psychophysiology")
missing.ones <- grep.on.Mesh("Sociology", "Social Sciences")
missing.ones <- grep.on.Mesh("Nursing", "Nursing")
missing.ones <- grep.on.Mesh("Nurse", "Nursing")
missing.ones <- grep.on.Mesh("Cognit", "Psychology")
missing.ones <- grep.on.Mesh("Mental", "Psychology")
missing.ones <- grep.on.Mesh("Psychology", "Psychology")
missing.ones <- grep.on.Mesh("Psychiatry", "Psychiatry")
missing.ones <- grep.on.Mesh("Health Policy", "Public Health")
missing.ones <- grep.on.Mesh("Chemistry, Pharmaceutical", "Pharmacology")
missing.ones <- grep.on.Mesh("Chemistry", "Chemistry")
missing.ones <- grep.on.Mesh("Pharmacogenetics", "Pharmacology")
missing.ones <- grep.on.Mesh("Pharmaceutical", "Pharmacology")
missing.ones <- grep.on.Mesh("Biomedical and Dental Materials", "Biomedical Engineering")
missing.ones <- grep.on.Mesh("Dentistry", "Dentistry")
missing.ones <- grep.on.Mesh("Dental", "Dentistry")
missing.ones <- grep.on.Mesh("Genom", "Computational Biology")
missing.ones <- grep.on.Mesh("Cardiovascular Diseases", "Cardiology")
missing.ones <- grep.on.Mesh("Heart", "Cardiology")
missing.ones <- grep.on.Mesh("Nervous System Diseases", "Neurology")
missing.ones <- grep.on.Mesh("Endocrine System Diseases", "Endocrinology")
missing.ones <- grep.on.Mesh("Neurosciences", "Brain")
missing.ones <- grep.on.Mesh("Neuropsychology", "Psychology")
missing.ones <- grep.on.Mesh("Midwifery", "Midwifery")      # Note there are plenty of midwifery journals and it's not a PubMed category, so I decided to create a new category
missing.ones <- grep.on.title("Midwif", "Midwifery")
missing.ones <- grep.on.Mesh("Pain", "Medicine")
missing.ones <- grep.on.Mesh("Forensic Medicine", "Medicine")
missing.ones <- grep.on.Mesh("Genetic Therapy", "Genetics, Medical")
missing.ones <- grep.on.Mesh("Orthoped", "Orthopedics")
missing.ones <- grep.on.Mesh("Neonatal", "Obstetrics")
missing.ones <- grep.on.Mesh("Nutritional Physiological Phenomena", "Nutritional Sciences")
missing.ones <- grep.on.Mesh("Epilepsy", "Neurology")
missing.ones <- grep.on.Mesh("Psychotherapy", "Psychology")
missing.ones <- grep.on.Mesh("Autis", "Psychiatry")
missing.ones <- grep.on.Mesh("Osteoporosis", "Osteopathic Medicine")
missing.ones <- grep.on.Mesh("Neoplasms", "Neoplasms")


# Now apply some broad strokes to try to mop up the remainder of unclassified journals
missing.ones <- grep.on.title("Cell", "Cell Biology")
missing.ones <- grep.on.title("Genet", "Genetics")
missing.ones <- grep.on.title("Chem", "Chemistry")
missing.ones <- grep.on.title("Drug", "Pharmacology")
missing.ones <- grep.on.title("Surg", "General Surgery")
missing.ones <- grep.on.title("Chem", "Chemistry")
missing.ones <- grep.on.title("Health", "Medicine")
missing.ones <- grep.on.title("Pathol", "Medicine")
missing.ones <- grep.on.title("Rehabil", "Medicine")
missing.ones <- grep.on.title("Neuro", "Medicine")
missing.ones <- grep.on.Mesh("Med", "Medicine")
missing.ones <- grep.on.title("Med", "Medicine")

# sort(unlist(tapply(missing.ones$n[is.na(missing.ones$broad.journal.headings)], missing.ones$journal[is.na(missing.ones$broad.journal.headings)], function(x) x))) # Some code used for bug checking the above

# Fill in the new research discipline information from 'missing.ones'
culled.metadata$broad.journal.headings[missing.ones$row] <- missing.ones$broad.journal.headings
rm(missing.ones)

# Now get rid of some of the smaller PubMed disciplines by merging with others in a more sensible grouping

# The category "technology" sounds vague, and in fact all of the journals in it could easily be placed in one of the other categories
culled.metadata$broad.journal.headings[grep("Enzyme Microb Technol", culled.metadata$short.title)] <- "Biotechnology"
culled.metadata$broad.journal.headings[grep("J Long Term Eff Med Implants", culled.metadata$short.title)] <- "Biomedical Engineering"
culled.metadata$broad.journal.headings[grep("Med Device Technol", culled.metadata$short.title)] <- "Biomedical Engineering"
culled.metadata$broad.journal.headings[grep("Soc Stud Sci", culled.metadata$short.title)] <- "Sociology"
culled.metadata$broad.journal.headings[grep("Zhongguo Yi Liao Qi Xie Za Zhi", culled.metadata$short.title)] <- "Biomedical Engineering"

# The following categories have fewer than 10 journals, so let's merge with a larger category
culled.metadata$broad.journal.headings[grep("Osteopathic Medicine", culled.metadata$broad.journal.headings)] <- "Medicine"  # only 1 journal
culled.metadata$broad.journal.headings[grep("Podiatry", culled.metadata$broad.journal.headings)] <- "Medicine" # only 2 journals
culled.metadata$broad.journal.headings[grep("Teratology", culled.metadata$broad.journal.headings)] <- "Medicine" # only 2 journals
culled.metadata$broad.journal.headings[grep("Aerospace Medicine", culled.metadata$broad.journal.headings)] <- "Medicine" # only 4 journals
culled.metadata$broad.journal.headings[grep("Military Medicine", culled.metadata$broad.journal.headings)] <- "Medicine" 
culled.metadata$broad.journal.headings[grep("Disaster Medicine", culled.metadata$broad.journal.headings)] <- "Medicine" # only 2 journals
culled.metadata$broad.journal.headings[grep("Bacteriology", culled.metadata$broad.journal.headings)] <- "Microbiology" # only 2 journals
culled.metadata$broad.journal.headings[grep("Chemistry, Clinical", culled.metadata$broad.journal.headings)] <- "Chemistry" # only 4 journals
culled.metadata$broad.journal.headings[grep("Family Planning Services", culled.metadata$broad.journal.headings)] <- "Reproductive Medicine" # only 4 journals
culled.metadata$broad.journal.headings[grep("Histocytochemistry", culled.metadata$broad.journal.headings)] <- "Cell Biology" # only 4 journals


# Here are the journals with no category assigned, in order of number of papers
sort(tapply(culled.metadata$n[is.na(culled.metadata$broad.journal.headings)], culled.metadata$short.title[is.na(culled.metadata$broad.journal.headings)], sum))

######################################################################################################
# GETTING IMPACT FACTOR FOR EACH OF THE JOURNALS 
######################################################################################################

# Load in the impact factor data from ISI
IFdata <- read.csv("../data/Impact factors 2015 from Thompson Reuters JCR.csv", stringsAsFactors = F)
IFdata <- IFdata[!duplicated(IFdata$ISO_ABBREV), ] # Remove the duplicate rows (for some reason, Thompson Reuters have included many journals multiple times - the info seems to be the same for each duplicated entry though)
IFdata$ISO_ABBREV <- gsub("[.]", "", IFdata$ISO_ABBREV) # Get rid of periods in journal names, as in my dataset

# Correct some of the non-matching names in the IF dataset. I wish they'd all use the same journal abbreviations, or that ISSN was a reliable way to do this!
IFdata$ISO_ABBREV[IFdata$ISO_ABBREV == "BMJ-British Medical Journal"] <- "BMJ"
IFdata$ISO_ABBREV[IFdata$ISO_ABBREV == "JAIDS"] <- "J Acquir Immune Defic Syndr"
IFdata$ISO_ABBREV[IFdata$ISO_ABBREV == "Chem Commun"] <- "Chem Commun (Camb)"
IFdata$ISO_ABBREV[IFdata$ISO_ABBREV == "Angew Chem-Int Edit"] <- "Angew Chem Int Ed Engl"
IFdata$ISO_ABBREV[IFdata$ISO_ABBREV == "JAMA-J Am Med Assoc"] <- "JAMA"
IFdata$ISO_ABBREV[IFdata$ISO_ABBREV == "Appl Optics"] <- "Appl Opt"
IFdata$ISO_ABBREV[IFdata$ISO_ABBREV == "ChemistryOpen"] <- "Chemistry"
IFdata$ISO_ABBREV[IFdata$ISO_ABBREV == "Can Med Assoc J"] <- "CMAJ"
IFdata$ISO_ABBREV[IFdata$ISO_ABBREV == "Nature Genet"] <- "Nat Genet"
culled.metadata$short.title[culled.metadata$short.title == "Proc Biol Sci"] <- "Proc R Soc B Biol Sci"  # oh, Proc B. Everyone has a different name for it
IFdata$ISO_ABBREV[IFdata$ISO_ABBREV == "Proc R Soc B-Biol Sci"] <- "Proc R Soc B Biol Sci"
IFdata$ISO_ABBREV[IFdata$ISO_ABBREV == "Phys Rev E"] <- "Phys Rev E Stat Nonlin Soft Matter Phys"
IFdata$ISO_ABBREV[IFdata$ISO_ABBREV == "J Bone Joint Surg-Am Vol"] <- "J Bone Joint Surg Am"

# In the PubMed data, these 'journal families' are not listed under separate names as they are in the IF data, so let's merge them and take the mean impact factor
averageIF <- function(journal.family, IFdata){
  average <- mean(IFdata$IMPACT_FACTOR[grep(journal.family, IFdata$ISO_ABBREV)])
  IFdata <- rbind(IFdata, rep(NA, ncol(IFdata))); IFdata[nrow(IFdata), 2] <- journal.family
  IFdata$IMPACT_FACTOR[IFdata$ISO_ABBREV == journal.family] <- average
  return(IFdata)
}

IFdata <- averageIF("J Phys Chem", IFdata)
IFdata <- averageIF("Biochim Biophys Acta", IFdata)
IFdata <- averageIF("J Phys Chem", IFdata)
IFdata <- averageIF("Acta Crystallogr", IFdata)
IFdata <- averageIF("J Chromatogr", IFdata)
IFdata <- averageIF("Am J Med Genet", IFdata)


# First fill in by exact match of the ISSN
culled.metadata$IF <- IFdata$IMPACT_FACTOR[match(culled.metadata$ISSN, IFdata$ISSN)] # First fill in by exact match of the title
# Where that failed, fill in based on non-case sensitive journal title
culled.metadata$IF[is.na(culled.metadata$IF)] <- IFdata$IMPACT_FACTOR[match(tolower(culled.metadata$short.title[is.na(culled.metadata$IF)]), tolower(IFdata$ISO_ABBREV))]
# Where THAT failed, fill in based on non-case sensitive journal title with the dashes removed
culled.metadata$IF[is.na(culled.metadata$IF)] <- IFdata$IMPACT_FACTOR[match(tolower(culled.metadata$short.title[is.na(culled.metadata$IF)]), gsub("-", "", tolower(IFdata$ISO_ABBREV)))]

culled.metadata[is.na(culled.metadata$IF), c(3,9)][order(culled.metadata[is.na(culled.metadata$IF), c(9)]), ]


# Write the cleaned-up journal metadata to disk
write.csv(culled.metadata, file = "../outputs/Journal metadata.csv", row.names = F)



###########################################################################
# MINE THE PUBMED XML FILES AGAIN, AND GET THE AFFILIATIONS THIS TIME
# I didn't get them the first time around, as it only occurred to me later
###########################################################################

setwd("~/Desktop/Gender Pubmed mining project/R scripts")

# devtools::install_github("ironholds/poster")
#library(poster)
library(XML)
library(stringr)
require(reshape2)
require(dplyr)
source("Data making functions.R")

# Parse the XML files to retreive the authors' addresses. This is a very big job 
# This function makes use of the 'poster' package to pre-process the addresses into fields (house, city, state, country)
# Though 'poster' very often fails to retrieive the country etc, in these cases it just dumps all of the address into the 'house' field, so no information is lost
parse.many.xml.affiliations("../data/Update files - raw xml", "../data/Update files - parsed affiliations", focal.chunk = 1, num.chunks = 3)
parse.many.xml.affiliations("../data/Main files - raw xml", "../data/Main files - parsed affiliations", focal.chunk = 1, num.chunks = 3)

parse.many.xml.affiliations("../data/Update files - raw xml", "../data/Update files - parsed affiliations", focal.chunk = 2, num.chunks = 3)
parse.many.xml.affiliations("../data/Main files - raw xml", "../data/Main files - parsed affiliations", focal.chunk = 2, num.chunks = 3)

parse.many.xml.affiliations("../data/Update files - raw xml", "../data/Update files - parsed affiliations", focal.chunk = 3, num.chunks = 3)
parse.many.xml.affiliations("../data/Main files - raw xml", "../data/Main files - parsed affiliations", focal.chunk = 3, num.chunks = 3)


# make a list of the most common countries found in the affiliation field (there are c. 118 real countries that show up semi-regularly; after that, it's all words that erroneously ended up in the "country" field)
countries <- get.unique.countries(c("../data/Update files - parsed affiliations", "../data/Main files - parsed affiliations"))
acceptable.countries <- c("japan", tail(parse.countries(countries, type = "summary"), 117)[,1])
acceptable.countries <- acceptable.countries[!(acceptable.countries %in% c("kingdom", "united"))] # check these are not in there

# Parse the 'state' field, and make a list of country-state associations
states <- get.unique.states(c("../data/Update files - parsed affiliations", "../data/Main files - parsed affiliations"))
state.countries <- states[!is.na(states$country), c(1,3)]

# Parse the 'city' field, and make a list of country-city associations
cities <- get.unique.cities(c("../data/Update files - parsed affiliations", "../data/Main files - parsed affiliations"))
city.countries <- cities[!is.na(cities$country), c(1,3)]


# Finally make the final affiliations dataset. Two columns - the pubmed id, and the author countries. This function also uses information from the 'house' address field
make.address.datafile(c("../data/Update files - parsed affiliations", "../data/Main files - parsed affiliations"), over.write = F, acceptable.countries, state.countries, city.countries)


affiliations <- do.call("rbind", lapply(list.files("../temp", full.names = T), read.csv, stringsAsFactors = F)) # zip the files together
affiliations <- affiliations[!duplicated(affiliations$pmid), ] # Few duplicates to remove
affiliations <- affiliations[!is.na(affiliations$country), ] # A typo stopped the NAs getting properly removed
affiliations$country <- str_replace_all(affiliations$country, "NA arab emirates", "united arab emirates")
write.csv(affiliations, file = "../outputs/affiliations.csv", row.names = F)

######################################################################################################
# ADDING AUTHOR GENDER INFORMATION FROM GENDERISE.IO TO THE PUBMED DATA
######################################################################################################

# Add the address data to each paper
pubmed <- read.csv("../outputs/Pubmed data for analysis.csv", stringsAsFactors = F) %>% select(pmid, forenames) # load the big data.frame
affiliations <- read.csv( "../outputs/affiliations.csv", stringsAsFactors = F) # Load the author address data
pubmed$country <- affiliations$country[match(pubmed$pmid, affiliations$pmid)]
rm(affiliations)
gc()

chunk.size <- 15000 # higher numbers cause memory problems
nChunks <- ceiling(nrow(pubmed) / chunk.size)
starts <- seq(1,1+(nChunks-1)*chunk.size, by = chunk.size)
ends <- seq(chunk.size, nChunks*chunk.size, by = chunk.size); ends[length(ends)] <- nrow(pubmed)
for(i in 1:length(starts)) write.csv(pubmed[starts[i]:ends[i], ], file = paste("../temp/names.to.gender", i, ".csv", sep = ""), row.names = F)
rm(pubmed)

all.names <- read.csv("../data/genderize.names.mastersheet.csv", stringsAsFactors = F) # load masses of names (and countries they correspond to) from genderise.io
all.names <- all.names[!is.na(all.names$prob),]
all.names$prob[all.names$gender == "M"] <- 1 - all.names$prob[all.names$gender == "M"]
all.names$country[is.na(all.names$country)] <- "unspec" # Change NAs to "country-unspecified"

files <- list.files("../temp", full.names = T)[str_detect(list.files("../temp", full.names = T), "names.to.gender")]

# Now loop over each temporary file, and add the gender information from 'all.names' in a country-sensitive manner (e.g. if we know the person's affiliation is from a particular country, try to use the gender information specific to that country where available - otherwise, use the gender info from all countries combined)
for(i in 1:length(files)){
  print(paste("Doing chunk ", i, " of ", length(files), ".", sep=""))
  df <- read.csv(files[i], stringsAsFactors = F)
  
  author.lists <- melt(strsplit(df$forenames, split = "_")) # first col ("value") is a name, second column ("L1") is the pmid of its paper
  names(author.lists) <- c("name", "pmid") # give them clearer names 
  author.lists$pmid <- df$pmid[author.lists$pmid]
  author.lists$gender <- "U"
  author.lists$gender.score <- NA
  author.lists$country <- str_replace(df$country[match(author.lists$pmid, df$pmid)], "NA", "unspec") # adds the country or countries for that name. Change NA to unspec to facilitate matching
  author.lists$country[is.na(author.lists$country)] <- "unspec"
  unique.pmids <- sort(unique(author.lists$pmid))
  country.n <- as.numeric(tapply(author.lists$country, author.lists$pmid, function(x) 1+str_count(x[1], "_"))) # count the number of countries listed for each paper
  author.n <- as.numeric(tapply(author.lists$name, author.lists$pmid, length)) # count the number of authors listed for each paper
  
  for(j in 1:length(unique.pmids)){
    # If only one country is listed for this paper, use that country for all authors. So, there is nothing to do in this for-loop and it will move to next paper
    
    if(country.n[j] != 1){   # Otherwise, if there is more than one country listed...
      if(country.n[j] == author.n[j]){ # And the number of countries listed is the same as the number of authors...
        author.lists$country[author.lists$pmid == unique.pmids[j]] <- unlist(strsplit(author.lists$country[author.lists$pmid == unique.pmids[j]][1], split = "_")) # Add in one country per author
      }
      else if(country.n[j] < author.n[j]){ # Or if the number of countries listed is LESS than the number of authors...
        author.lists$country[author.lists$pmid == unique.pmids[j]] <- NA  # Be conservative and assume we don't know where the authors are from
      }
    }
  }
  
  # Add in the gender information in a country-sepcific manner, where possible
  matches <- match(paste(tolower(author.lists$name), author.lists$country), paste(all.names$name, all.names$country))
  author.lists$gender <- all.names$gender[matches]
  author.lists$gender[is.na(author.lists$gender)] <- "U"
  author.lists$gender.score <- all.names$prob[matches]
  
  # When that fails, add in the gender information in a non-country-sepcific manner
  matches <- match(tolower(author.lists$name[author.lists$gender == "U"]), all.names$name)
  author.lists[author.lists$gender == "U", names(author.lists) %in% c("gender", "gender.score")] <- cbind(all.names$gender[matches], all.names$prob[matches])
  author.lists$gender[is.na(author.lists$gender)] <- "U"
  author.lists$gender.score <- as.numeric(author.lists$gender.score)
  
  xx <- tapply(author.lists$gender, author.lists$pmid, function(x) paste0(x, collapse="")) # Add the info from this chunk back into the focal chunk of data from the 'temp' folder
  df$gender[match(names(xx), df$pmid)] <- as.character(xx)
  xx <- tapply(round(author.lists$gender.score, 3), author.lists$pmid, function(x) paste0(x, collapse="_"))
  df$gender.score[match(names(xx), df$pmid)] <- as.character(xx)
  
  write.csv(df, file = files[i], row.names = F)
}

# Now associate the gender info with the pubmed data
gender.info <- do.call("rbind", lapply(files, read.csv, stringsAsFactors = F)) %>% select(pmid, gender, gender.score)
pubmed <- read.csv("../outputs/Pubmed data for analysis.csv", stringsAsFactors = F) 
pubmed$gender <- gender.info$gender[match(pubmed$pmid, gender.info$pmid)]
pubmed$gender.score <- gender.info$gender.score[match(pubmed$pmid, gender.info$pmid)]
write.csv(pubmed, "../outputs/Pubmed data with gender from genderize.csv", row.names = F)


#################################################################
# Writing the data to a database on disk
#################################################################
library(dplyr)
setwd("~/Desktop/Gender Pubmed mining project/R scripts")

# Add the address data to each paper
pubmed <- read.csv("../outputs/Pubmed data with gender from genderize.csv", stringsAsFactors = F) # load the big data.frame
affiliations <- read.csv( "../outputs/affiliations.csv", stringsAsFactors = F) # Load the author address data
pubmed$country <- affiliations$country[match(pubmed$pmid, affiliations$pmid)]
rm(affiliations)
gc()

# Load up the journal meta-data
journal.metadata <- read.csv("../outputs/Journal metadata.csv", stringsAsFactors = F) # load the journal data
names(journal.metadata)[names(journal.metadata) %in% "broad.journal.headings"] <- "Discipline" # give it a clearer name 
names(journal.metadata)[names(journal.metadata) %in% "country"] <- "journal.country" # give it a clearer name 

# Add some final post-hoc changes to the journal disciplines
journal.metadata$Discipline[journal.metadata$Discipline == "Science"] <- "Multidisciplinary" # Re-name a couple of categories for clarity/brevity
journal.metadata$Discipline[journal.metadata$Discipline == "Acquired Immunodeficiency Syndrome"] <- "AIDS" # Give it a shorter name, for easy plotting
journal.metadata$Discipline[journal.metadata$Discipline == "Photography"] <- "Medicine" # There is only one "photography" journal -  Journal of Visual Communication in Medicine. Let's just put this in "medicine"
journal.metadata$Discipline[journal.metadata$Discipline == "Sociology"] <- "Multidisciplinary" # There is only one "sociology" journal -  Social studies of science. Let's just put this in "multidisciplinary" (this seems apt based on their website)

# Add a discipline column to the 'papers' table, so that we can run indexed searches by discipline later
pubmed <- left_join(pubmed, journal.metadata %>% rename(journal = short.title, discipline = Discipline) %>% select(discipline, journal))

# Create and index a database using the data frames
pubmed_db <- src_sqlite("../data for analysis/pubmed_gender_data.sqlite3", create = T)
copy_to(pubmed_db, pubmed, "papers", temporary = FALSE, indexes = list("journal", "country", "discipline"))
copy_to(pubmed_db, journal.metadata, "journals", temporary = FALSE, indexes = list("Discipline"))
rm(list = c("pubmed", "journal.metadata")) # free up memory


##########################################################################################################################################################################
# Post-hoc error corrections. Somehow, there are >1000 journals in the dataset which are not listed in the journal metadata file. Must have been a bug earlier on. At any rate let's fix it now
# I think part of the problem could be that Pubmed's list of all the journals they index is not complete. Also, people mis-spell journal names, add whitespace, or use different acronyms for the same journal, making it tricky. A final issue is that PubMed duplicated the short names of a handful of journals, and gave them more than one different discipline. I fixed this by picking one discipline in these cases.
##########################################################################################################################################################################

pubmed_db <- src_sqlite("../data for analysis/pubmed_gender_data.sqlite3", create = F)
pubmed_sqlite <- tbl(pubmed_db, "papers")
journals_sqlite <- tbl(pubmed_db, "journals")

unique.papers <- (pubmed_sqlite %>% distinct(journal) %>% collect(n = Inf) %>% as.data.frame)[,1]
unique.journals <- (journals_sqlite %>% distinct(short.title) %>% collect(n = Inf) %>% as.data.frame)[,1]
unique.papers[!(unique.papers%in%unique.journals)] 

# Many of the journal names still have whitespace at the end, so remove them now
journal.column <- (pubmed_sqlite %>% select(journal) %>% collect(n = Inf) %>% as.data.frame)[,1]
last.char <- substr(journal.column, nchar(journal.column), nchar(journal.column))
journal.column[last.char == " "] <- substr(journal.column[last.char == " "], 1, nchar(journal.column[last.char == " "]) - 1)
rm(last.char)
journal.column[journal.column == "Proc Biol Sci"] <- "Proc Roy Soc B" # This is my favourite acronym for it! Seems there are loads of versions.

unique.papers <- sort(unique(journal.column))
unique.papers[!(unique.papers%in%unique.journals)] 

new.metadata <- journals_sqlite %>% collect(n = Inf) %>% as.data.frame
topics <- sort(unique(new.metadata$Discipline))
addition <- new.metadata[1:length(unique.papers[!(unique.papers%in%unique.journals)]), ] 
addition[1:nrow(addition), 1:ncol(addition)] <- NA
addition$short.title <- unique.papers[!(unique.papers%in%unique.journals)] 
new.metadata <- rbind(new.metadata,addition)
                      
# So that this works with the code that I already wrote much futher up, making an object called "missing.ones" just like I did originally
counts <- melt(table(journal.column)) %>% mutate(journal.column = as.character(journal.column))
missing.ones <- with(new.metadata, data.frame(journal = short.title, row = 1:nrow(new.metadata), stringsAsFactors = F))
missing.ones$n <- counts$value[match(missing.ones$journal, counts$journal.column)]
missing.ones <- missing.ones[is.na(new.metadata$Discipline), ]

###### Now scroll up and run all the lines of code with grep.on.title in them. Lots of the coding was already done above. #####

missing.ones$broad.journal.headings <- NA
missing.ones <- grep.on.title("Biol", "Biology")
missing.ones <- grep.on.title("Worm", "Biology")
missing.ones <- grep.on.title("Evol", "Biology")
missing.ones <- grep.on.title("Biochem", "Biochemistry")
missing.ones <- grep.on.title("Zoo", "Zoology")
missing.ones <- grep.on.title("Hepatol", "Medicine")
missing.ones <- grep.on.title("Psychol", "Psychology")
missing.ones <- grep.on.title("Cogn", "Psychology")
missing.ones <- grep.on.title("Behav", "Psychology")
missing.ones <- grep.on.title("Thromb", "Hematology")
missing.ones <- grep.on.title("Hematol", "Hematology")
missing.ones <- grep.on.title("Hepat", "Medicine")
missing.ones <- grep.on.title("Circ", "Cardiology")
missing.ones <- grep.on.title("Heart", "Cardiology")
missing.ones <- grep.on.title("Plant", "Botany")
missing.ones <- grep.on.title("Photosyn", "Botany")
missing.ones <- grep.on.title("Pediatr", "Pediatrics")
missing.ones <- grep.on.title("Virol", "Virology")
missing.ones <- grep.on.title("Thyroid", "Endocrinology")
missing.ones <- grep.on.title("Endocr", "Endocrinology")
missing.ones <- grep.on.title("Trop Parasitol", "Tropical Medicine")
missing.ones <- grep.on.title("Malar J", "Tropical Medicine")
missing.ones$broad.journal.headings[missing.ones$journal == "Trends Ecol Evol (Amst)"] <- "Biology" 
missing.ones <- grep.on.title("PLoS ONE", "Multidisciplinary")
missing.ones <- grep.on.title("Proc Natl Acad Sci USA", "Multidisciplinary")
missing.ones <- grep.on.title("PeerJ", "Multidisciplinary")
missing.ones <- grep.on.title("Nippon Rinsho", "Medicine")
missing.ones <- grep.on.title("Nippon Naika Gakkai Zasshi", "Medicine")
missing.ones <- grep.on.title("Ugeskr Laeg", "Medicine")
missing.ones <- grep.on.title("Proc Roy Soc B", "Biology")
missing.ones <- grep.on.title("EMBO", "Cell Biology")
missing.ones <- grep.on.title("Retina", "Ophthalmology")
missing.ones <- grep.on.title("FEBS", "Biochemistry")
missing.ones <- grep.on.title("FASEB J", "Cell Biology")
missing.ones <- grep.on.title("Dis", "Medicine")
missing.ones <- grep.on.title("Manag Care", "Medicine")
missing.ones <- grep.on.title("Respir", "Pulmonary Medicine")
missing.ones <- grep.on.title("Clin", "Medicine")
missing.ones <- grep.on.title("Genome", "Genetics")
missing.ones <- grep.on.title("Prz Lek", "Medicine")
missing.ones <- grep.on.title("Diagn", "Medicine")
missing.ones <- grep.on.title("Orthop", "Orthopedics")
missing.ones <- grep.on.title("Meth Enzymol", "Biochemistry")
missing.ones <- grep.on.title("Adv Mater Weinheim", "Biophysics")
missing.ones <- grep.on.title("Springerplus", "Multidisciplinary")
missing.ones <- grep.on.title("Nanoscale Res Lett", "Multidisciplinary")
missing.ones <- grep.on.title("MMWR Morb Mortal Wkly Rep", "Multidisciplinary")
missing.ones <- grep.on.title("Food Funct", "Nutritional Sciences")
missing.ones <- grep.on.title("J Food Prot", "Nutritional Sciences")
missing.ones <- grep.on.title("Stand Genomic Sci", "Genetics")
missing.ones <- grep.on.title("J Food Sci Technol", "Nutritional Sciences")
missing.ones <- grep.on.title("F1000Res", "Multidisciplinary")
missing.ones <- grep.on.title("NMR Biomed", "Biomedical Engineering")
missing.ones <- grep.on.title("J Neural Transm (Vienna)", "Neurology")
missing.ones <- grep.on.title("J Phys Ther Sci", "Physical and Rehabilitation Medicine")
missing.ones <- grep.on.title("Data Brief", "Multidisciplinary")
missing.ones <- grep.on.title("Theriogenology", "Biology")
missing.ones <- grep.on.title("Biotechnol Biofuels", "Biotechnology")
missing.ones <- grep.on.title("Haematol", "Hematology")
missing.ones <- grep.on.title("Cytogenet", "Cell Biology")
missing.ones <- grep.on.title("Mycobiology", "Microbiology")
missing.ones <- grep.on.title("Photochem", "Biology")
missing.ones <- grep.on.title("Biomed Mater", "Biomedical Engineering")
missing.ones <- grep.on.title("Mycorrhiza", "Botany")


missing.ones %>% filter(is.na(broad.journal.headings)) %>% select(journal, n) %>% arrange(n) # All the ones I have not done manually have < 3000 papers

new.metadata$Discipline[missing.ones$row] <- missing.ones$broad.journal.headings
new.metadata$n[missing.ones$row] <- missing.ones$n


# Fix these journals, which PubMed assigned more than category to (happens occasionally when there were multiple different journals that have the same shortened name)
new.metadata$Discipline[grep("Acta Crystallogr", new.metadata$short.title)] <- "Chemistry Techniques, Analytical"
new.metadata$Discipline[grep("Biophys", new.metadata$short.title)] <- "Biophysics"
new.metadata$Discipline[grep("Chin Med", new.metadata$short.title)] <- "Complementary Therapies"
new.metadata$Discipline[grep("Microbiol", new.metadata$short.title)] <- "Microbiology"
new.metadata$Discipline[grep("Med J", new.metadata$short.title)] <- "Medicine"
new.metadata$Discipline[grep("Neuro", new.metadata$short.title)] <- "Neurology"
new.metadata$Discipline[grep("Pharm", new.metadata$short.title)] <- "Pharmacology"
new.metadata$Discipline[grep("Psychophar", new.metadata$short.title)] <- "Psychopharmacology"
new.metadata$Discipline[grep("Complement Alternat", new.metadata$short.title)] <- "Complementary Therapies"
new.metadata$Discipline[grep("Vet J", new.metadata$short.title)] <- "Veterinary Medicine"
new.metadata$Discipline[grep("Physiol", new.metadata$short.title)] <- "Physiology"
new.metadata$Discipline[grep("Immunol", new.metadata$short.title)] <- "Immunology"
new.metadata$Discipline[grep("Derm", new.metadata$short.title)] <- "Dermatology"
new.metadata$Discipline[grep("Med Assoc", new.metadata$short.title)] <- "Medicine"
new.metadata$Discipline[grep("Biomed", new.metadata$short.title)] <- "Medicine"
new.metadata$Discipline[grep("[Oo]nco", new.metadata$short.title)] <- "Neoplasms"
new.metadata$Discipline[grep("Psychiatr", new.metadata$short.title)] <- "Psychiatry"
new.metadata$Discipline[grep("Surg", new.metadata$short.title)] <- "General Surgery"
new.metadata$Discipline[grep("Anaesth", new.metadata$short.title)] <- "Anesthesiology"
new.metadata$Discipline[grep("Anesthes", new.metadata$short.title)] <- "Anesthesiology"
new.metadata$Discipline[grep("Urol", new.metadata$short.title)] <- "Urology"
new.metadata$Discipline[grep("Genet", new.metadata$short.title)] <- "Genetics"
new.metadata$Discipline[grep("Paediatr", new.metadata$short.title)] <- "Pediatrics"
new.metadata$Discipline[grep("Pediatr", new.metadata$short.title)] <- "Pediatrics"
new.metadata$Discipline[grep("Ophthal", new.metadata$short.title)] <- "Ophthalmology"
new.metadata$Discipline[grep("Otolar", new.metadata$short.title)] <- "Otolaryngology"
new.metadata$Discipline[grep("Gastroent", new.metadata$short.title)] <- "Gastroenterology"
new.metadata$Discipline[grep("Org Chem", new.metadata$short.title)] <- "Chemistry"
new.metadata$Discipline[grep("Tissue Eng Part", new.metadata$short.title)] <- "Medicine"
new.metadata$Discipline[grep("Midwifery", new.metadata$short.title)] <- "Midwifery"
new.metadata$Discipline[grep("Antibiot", new.metadata$short.title)] <- "Anti-Infective Agents"
new.metadata$Discipline[grep("Antivir", new.metadata$short.title)] <- "Anti-Infective Agents"
new.metadata$Discipline[grep("Med", new.metadata$short.title)] <- "Medicine"
new.metadata$Discipline[grep("Environ Health", new.metadata$short.title)] <- "Environmental Health"
new.metadata$Discipline[grep("Nutr", new.metadata$short.title)] <- "Nutritional Sciences"
new.metadata$Discipline[grep("[Tt]rauma", new.metadata$short.title)] <- "Traumatology"
new.metadata$Discipline[grep("Pain", new.metadata$short.title)] <- "Medicine"
new.metadata$Discipline[grep("Cancer", new.metadata$short.title)] <- "Neoplasms"
new.metadata$Discipline[grep("Bioinform", new.metadata$short.title)] <- "Bioinformatics"
new.metadata$Discipline[grep("Biotech", new.metadata$short.title)] <- "Biotechnology"
new.metadata$Discipline[grep("BioTech", new.metadata$short.title)] <- "Biotechnology"
new.metadata$Discipline[grep("Cardiol", new.metadata$short.title)] <- "Cardiology"
new.metadata$Discipline[grep("Gynecol", new.metadata$short.title)] <- "Gynecology"
new.metadata$Discipline[grep("Public Health", new.metadata$short.title)] <- "Public Health"
new.metadata$Discipline[grep("Stem Cell", new.metadata$short.title)] <- "Cell Biology"
new.metadata$Discipline[grep("Speech Lang", new.metadata$short.title)] <- "Speech-Language Pathology"

new.metadata$Discipline[new.metadata$short.title == "Am Indian Alsk Native Ment Health Res"] <- "Psychology"
new.metadata$Discipline[new.metadata$short.title == "Biom"] <- "Statistics as Topic"
new.metadata$Discipline[new.metadata$short.title == "Clin Exp Hypertens"] <- "Vascular Diseases"
new.metadata$Discipline[new.metadata$short.title == "Hosp"] <- "Hospitals"
new.metadata$Discipline[new.metadata$short.title == "J Comp Psychol"] <- "Psychology"
new.metadata$Discipline[new.metadata$short.title == "J Environ Sci Health"] <- "Environmental Health"
new.metadata$Discipline[new.metadata$short.title == "J Epidemiol Community Health"] <- "Epidemiology"
new.metadata$Discipline[new.metadata$short.title == "Manag Care"] <- "Health Services"
new.metadata$Discipline[new.metadata$short.title == "Med Mycol"] <- "Microbiology"
new.metadata$Discipline[new.metadata$short.title == "Mycopathologia"] <- "Microbiology"
new.metadata$Discipline[new.metadata$short.title == "Pharm Biol"] <- "Pharmacology"
new.metadata$Discipline[new.metadata$short.title == "Q J Exp Psychol"] <- "Psychology"
new.metadata$Discipline[new.metadata$short.title == "Reprod Nutr Dev"] <- "Nutritional Sciences"
new.metadata$Discipline[new.metadata$short.title == "Soc Sci Med"] <- "Social Sciences"
new.metadata$Discipline[new.metadata$short.title == "Eukaryotic Cell"] <- "Cell Biology"
new.metadata$Discipline[new.metadata$short.title == "Mol Reprod Dev"] <- "Biology"
new.metadata$Discipline[new.metadata$short.title == "Mol Ecol Resour"] <- "Biology"
new.metadata$Discipline[new.metadata$short.title == "DNA Repair (Amst)"] <- "Biology"
new.metadata$Discipline[new.metadata$short.title == "Mitochondrial DNA"] <- "Genetics"
new.metadata$Discipline[new.metadata$short.title == "Neural Regen Res"] <- "Neurology"
new.metadata$Discipline[new.metadata$short.title == "J Bodyw Mov Ther"] <- "Physical and Rehabilitation Medicine"

new.metadata$Discipline[new.metadata$Discipline == "Pharmacy"] <- "Pharmacology" # merge these disciplines

new.metadata %>% filter(is.na(Discipline)) %>% select(short.title, n) %>% arrange(n) # All the journals I have not done manually have < 3000 papers
sort(unique(new.metadata$Discipline))

# At this point I decided to leave the remaining journals' disciplines as "unclassified" in the interests of time (small returns in sample size for continued effort at this point)
new.metadata$Discipline[is.na(new.metadata$Discipline)] <- "Unclassified"
# write.csv(new.metadata, file = "../outputs/Journal metadata.csv", row.names = F) # write the modified metadata file to disk

pubmed <- pubmed_sqlite %>% collect(n = Inf) %>% as.data.frame
pubmed$discipline <- new.metadata$Discipline[match(pubmed$journal, new.metadata$short.title)]
pubmed_db$con %>% db_drop_table(table = "papers") # remove the old papers table from the database
copy_to(pubmed_db, pubmed, "papers", temporary = FALSE, indexes = list("journal", "discipline", "country")) # Copy the revised table onto the database on disk
pubmed_db$con %>% db_drop_table(table = "journals") # remove the old journal info table from the database
copy_to(pubmed_db, new.metadata, "journals", temporary = FALSE, indexes = list("short.title", "Discipline")) # Copy the revised table onto the database on disk
rm(pubmed)

#####################################################################################
# Inspect the success rate of gendering, and make an additional gender column that
# classifies as "U" everyone for whom gender is not known with >= 95% confidence
# This column is called gender.95
#####################################################################################

# Load up the database
pubmed_db <- src_sqlite("../data for analysis/pubmed_gender_data.sqlite3", create = F)
pubmed_sqlite <- tbl(pubmed_db, "papers")

# Let's look at the proportion of authors whose names were gendered, with high and low confidence
x <- (pubmed_sqlite %>% select(gender.score) %>% collect(n = Inf) %>% as.data.frame())[,1]
split <- unlist(strsplit(x, split = "_"))
split[split == "NA"] <- NA
split <- as.numeric(split)

# We discard 100 - 58.677 - 31.823 = 9.50 % of the authors by getting rid of the authors with slightly ambiguously gendered names.
sum(split <= 0.05, na.rm=T) / sum(!is.na(split)) # 58.677% of the authors whose gender is not totally unknown have names that are >=95% male
sum(split >= 0.95, na.rm=T) / sum(!is.na(split)) # 31.823% of the authors whose gender is not totally unknown have names that are >=95% female

# Some more stats on the success rate of gendering:
sum(is.na(split)) / length(split) # There is no gender information at all for 19.48% authors (i.e. XXX)
sum(!is.na(split)) # By contrast, we have at least some gender information for XXX% authors (i.e. 39,226,243)
sum(!is.na(split) & (split <= 0.05 | split >= 0.95)) # And we have >95% confident gender information for XX% authors (i.e. 35,499,685)

# Let's make a new gender column, that classifies as "U" everyone for whose gender is not known with >= 95% confidence
split <- strsplit(x, split = "_")
paperID <- unlist(mapply(rep, 1:length(split), each = sapply(split, length)))
split <- unlist(split)
split[split == "NA"] <- NA
split <- as.numeric(split)
new.gender <- rep("U", length(split))
new.gender[split >= 0.95] <- "F"
new.gender[split <= 0.05] <- "M"
new.gender <- tapply(new.gender, paperID, function(x) paste0(x, collapse = ""))
rm(list = c("paperID", "x", "split"))

# Add the new gender column to the dataset. First load up the relevant table as a data.frame
pubmed_sqlite <- tbl(pubmed_db, "papers") %>% collect(n = Inf) %>% as.data.frame()
pubmed_sqlite$gender.95 <- new.gender # Add the gender column
pubmed_db$con %>% db_drop_table(table = "papers") # remove the old table from the database
copy_to(pubmed_db, pubmed_sqlite, "papers", temporary = FALSE, indexes = list("journal", "country")) # Copy the revised table onto the database on disk
pubmed_sqlite <- tbl(pubmed_db, "papers")
rm(new.gender)

#################################################################
# Make doubly sure there are no duplicated Pubmed IDs in the dataset
# And do some final clean-ups that slipped the net earlier
#################################################################

pubmed_db <- src_sqlite("../data for analysis/pubmed_gender_data.sqlite3", create = F)
pubmed_sqlite <- tbl(pubmed_db, "papers") %>% collect(n = Inf) %>% as.data.frame()
pubmed_sqlite$discipline[pubmed_sqlite$discipline == "Acquired Immunodeficiency Syndrome"] <- "AIDS" # rename this one
# Get rid of the erroneously-extracted countries "united" and "kingdom", which were missed earlier due to bugs (they are both ambiguous because of e.g. UK/Saudi and UK/USA). There are only a small number of each, so no big deal.
split.country <- strsplit(pubmed_sqlite$country, split = "_")   
split.country <- lapply(split.country, function(x) {
  x[x == "united"] <- "NA"
  x[x == "kingdom"] <- "NA"
  x})

# For some papers, my code found more countries than there are authors (maybe because of dual affiliation). To be simple/conservative, let's just say that all of these failed to extract the country
num.countries <- sapply(split.country, length)
num.authors <- nchar(pubmed_sqlite$gender)
split.country[num.countries > num.authors] <- "NA"
split.country <- unlist(lapply(split.country, function(x) paste0(x, collapse = "_")))
pubmed_sqlite$country <- (split.country) 
rm(list = c("split.country", "num.countries", "num.authors"))

# To match the journal meta.data file, make sure that journals whose name ends in a space have that space deleted    
num.char <- nchar(pubmed_sqlite$journal)
pubmed_sqlite$journal[substr(pubmed_sqlite$journal, num.char, num.char) == " "] <- substr(pubmed_sqlite$journal[substr(pubmed_sqlite$journal, num.char, num.char) == " "], 1,
                                                                                          nchar(pubmed_sqlite$journal[substr(pubmed_sqlite$journal, num.char, num.char) == " "])-1)
pubmed_sqlite$journal[pubmed_sqlite$journal == "Proc Biol Sci"] <- "Proc Roy Soc B"  # And use my preferred abbreviation for this one
rm(num.char)

pubmed_sqlite <- pubmed_sqlite[!duplicated(pubmed_sqlite$pmid), ] # remove the few duplicates
pubmed_db$con %>% db_drop_table(table = "papers") # remove the old table from the database
copy_to(pubmed_db, pubmed_sqlite, "papers", temporary = FALSE, indexes = list("journal", "discipline", "country")) # Copy the revised table onto the database on disk
pubmed_sqlite <- tbl(pubmed_db, "papers")


####################################################################
# Make a summary dataset to analyse, and put online. Has all the info
# for each combination of journal, discipline, country and author position
####################################################################

journals <- (pubmed_sqlite %>% select(journal) %>% distinct() %>% collect(n = Inf) %>% as.data.frame())[,1]

full.summary <- do.call("rbind", lapply(1:length(journals), function(journal.i){
  print(journals[journal.i])
  focal <- pubmed_sqlite %>% select(journal, discipline, country, date, gender.95) %>% filter(journal == journals[journal.i]) %>% collect(n=Inf) # get the data on the focal journal
  focal$country[is.na(focal$country)] <- "unknown" # note the unknown countries
  n.countries <- str_count(focal$country, "_") + 1 # count the countries (number of underscores, plus 1)
  focal$country[!(n.countries == 1 | n.countries == nchar(focal$gender.95))] <- "unknown" # To be safe, for papers with multiple countries where this number is not the same as number of authors, list all countries as 'unknown'
  focal$country <- str_replace_all(as.character(focal$country), "NA", "unknown")
  
  rm(n.countries)
  
  focal$date <- substr(focal$date, 7, 10) # convert to year
  years <- unique(focal$date)
  
  author.list <- lapply(1:length(years), function(year.i){
    focal.year <- focal %>% filter(date == years[year.i])
    
    split.country <- strsplit(focal.year$country, split = "_") # split the string holding the countries
    split.gender.95 <- strsplit(focal.year$gender.95, split = "") # split the string holding the genders
    
    if(sum(grepl("NA", unlist(split.country))) > 0) print(focal$country[str_detect(focal$country, "NA")])
    
    author.list.for.this.year <- lapply(1:length(split.gender.95), function(paper.i){
      focal.countries <- split.country[[paper.i]]
      focal.genders <- split.gender.95[[paper.i]]
      num.authors <- length(focal.genders)
      if(length(focal.countries) == 1) focal.countries <- rep(focal.countries, num.authors)
      
      if(num.authors == 1){
        author.list.for.this.paper <- paste(focal.countries, "Single", focal.genders, sep = "_")
      }
      
      if(num.authors == 2){
        author.list.for.this.paper <- c(paste(focal.countries[1], "First", focal.genders[1], sep = "_"),
                                        paste(focal.countries[num.authors], "Last", focal.genders[num.authors], sep = "_"))
        
      }
      
      if(num.authors > 2){
        
        middle.authors <- sapply(2:(num.authors - 1), function(xx) paste(focal.countries[xx], "Middle", focal.genders[xx], sep = "_"))
        
        author.list.for.this.paper <- c(paste(focal.countries[1], "First", focal.genders[1], sep = "_"),
                                        middle.authors,
                                        paste(focal.countries[num.authors], "Last", focal.genders[num.authors], sep = "_"))
      }
      
      return(melt(table(author.list.for.this.paper)))
    })
    
    author.list.for.this.year <- do.call("rbind", author.list.for.this.year) %>% group_by(author.list.for.this.paper) %>% summarise(n = sum(value)) %>% mutate(year = years[year.i])
    return(author.list.for.this.year)
  })
  
  author.list <- do.call("rbind", author.list) %>% group_by(author.list.for.this.paper, year) %>% summarise(n = sum(n)) %>% rename(author.type = author.list.for.this.paper)
  split <- do.call("rbind", strsplit(as.character(author.list$author.type), split = "_"))
  author.list <- data.frame(journal = journals[journal.i], country = split[,1], year = as.numeric(author.list$year), position = split[,2], gender = split[,3], n = author.list$n, stringsAsFactors = F)
  return(author.list)
  
})) %>% mutate(year = as.numeric(year), country = capitalise.countries(country)) %>% arrange(journal, country, year, position, gender)


write.csv(full.summary, file = "../data for analysis/Author counts for journal country and position.csv", row.names = F) # 265 MB uncompressed




########################################################################################################
# Make the file "genderizing.success.rate.by.country.csv", which we need for a couple of the plots.
# This is used to check if the gender assignment success rate is similar across countries 
########################################################################################################

pubmed_db <- src_sqlite("../data for analysis/pubmed_gender_data.sqlite3", create = F)
pubmed_sqlite <- tbl(pubmed_db, "papers")
journals_sqlite <- tbl(pubmed_db, "journals")

# Split the names, genders, and countries
columns <- pubmed_sqlite %>% select(forenames, country, gender.95) %>% collect(n = Inf) %>% as.data.frame()
split.names <- strsplit(columns$forenames, split = "_"); columns <- columns %>% select(-forenames)
split.country <- strsplit(columns$country, split = "_"); columns <- columns %>% select(-country)
split.gender <- strsplit(columns$gender.95, split = ""); rm(columns)

# Discard all papers where the country is just listed as NA
to.keep <- sapply(split.country, function(x) {
  if(length(x) == 1) {
    if(is.na(x)) return(FALSE)
   # if(x == "NA") return(FALSE)
  }
  return(TRUE)
})
to.keep[sapply(split.country, length) != sapply(split.gender, length)] <- FALSE  # Also discard all papers where the number of countries listed is not equal to the number of authors
split.names <- split.names[to.keep]
split.country <- split.country[to.keep]
split.gender <- split.gender[to.keep]
length(split.country) # There are 1,902,001 papers that fit the bill

split.names <- unlist(split.names)
split.gender <- unlist(split.gender)
split.country <- unlist(split.country)
to.I <- rep(FALSE, length(split.names))
to.I[nchar(split.names) == 1] <- TRUE # Change single-letter names to have the gender I, as these are initials
to.I[str_detect(str_replace(str_replace(split.names, "[-,/?!]", ""), "[.]", ""), "^[[:upper:]]+$")] <- TRUE # Also change the given names in ALL CAPS to I, as they are mostly initials
rm(split.names)

# Revisit the database to get some more columns
split.gender2 <- (pubmed_sqlite %>% select(gender) %>% collect(n = Inf) %>% as.data.frame())[,1]
split.gender2 <- split.gender2[to.keep]
num.authors <- nchar(split.gender2)
data.to.make.country.summary <- (pubmed_sqlite %>% select(date, discipline, journal, abstract.length, is.journal.article) %>% collect(n = Inf) %>% as.data.frame() )[to.keep, ]
data.to.make.country.summary$date <- substr(data.to.make.country.summary$date, 7, 10)
data.to.make.country.summary <- data.frame(gender = unlist(strsplit(split.gender2, split = "")),
                                           gender.95 = split.gender, 
                                           country = split.country,
                                           year = as.numeric(as.character(unlist(mapply(rep, data.to.make.country.summary$date, each = num.authors)))),
                                           discipline = unlist(mapply(rep, data.to.make.country.summary$discipline, each = num.authors)),
                                           journal = unlist(mapply(rep, data.to.make.country.summary$journal, each = num.authors)),
                                           include.in.rigorous = unlist(mapply(rep, data.to.make.country.summary$abstract.length > 0 & data.to.make.country.summary$is.journal.article == 1, each = num.authors)),
                                           stringsAsFactors = F)
for(i in c(1:3,5:6))data.to.make.country.summary[, i] <- as.character(data.to.make.country.summary[, i])
data.to.make.country.summary$gender[to.I] <- "I"
data.to.make.country.summary$gender.95[to.I] <- "I"
rm(list = c("to.keep", "to.I", "split.gender", "split.gender2", "split.country", "num.authors"))

length(data.to.make.country.summary$gender.95) # For gender.95, there are 8,252,264 authors where we know their country (applying the strict criterion to discard all papers where the number of addresses was different to number of authors) 
length(data.to.make.country.summary$gender[data.to.make.country.summary$gender != "U" & data.to.make.country.summary$gender != "I"]) # Of these, 6,735,080 have known gender with some confidence
length(data.to.make.country.summary$gender.95[data.to.make.country.summary$gender.95 != "U" & data.to.make.country.summary$gender.95 != "I"]) # Of these, 6,098,400 have known gender to >95% confidence

# Same stats for the rigorous data:
length(data.to.make.country.summary$gender.95[data.to.make.country.summary$include.in.rigorous]) # For gender.95, there are 4991946 authors where we know their country (applying the strict criterion to discard all papers where the number of addresses was different to number of authors)
length(data.to.make.country.summary$gender[data.to.make.country.summary$include.in.rigorous & data.to.make.country.summary$gender != "U" & data.to.make.country.summary$gender != "I"]) # Of these, 4,130,675 have known gender with some confidence
length(data.to.make.country.summary$gender.95[data.to.make.country.summary$include.in.rigorous & data.to.make.country.summary$gender.95 != "U" & data.to.make.country.summary$gender.95 != "I"]) # Of these, 3,745,092 have known gender to >95% confidence

# Make a master summary that counts author gender (all 4 kinds: M, F, unknown due name not recognised "U", and unknown due to initial, "I")
country.summary.by.year.and.journal <- data.to.make.country.summary %>% group_by(country, gender.95, year, journal, discipline) %>% summarise(count = n()) %>% as.data.frame()

# Same thing for the rigorous data - commented out as I didn't plot this in the paper. Very similar-looking
# country.summary.by.year.and.journal.rigorous <- data.to.make.country.summary %>% filter(include.in.rigorous) %>% group_by(country, gender.95, year, journal, discipline) %>% summarise(count = n()) %>% as.data.frame()


# Analyse these results

# Only 2.94% of authors use initials instead of a forename only
100 * sum(country.summary.by.year.and.journal$count[country.summary.by.year.and.journal$gender.95 == "I"]) / sum(country.summary.by.year.and.journal$count) 

# If I include all authors (including those using initials as forenames), then the name genderisation worked with >=95% success for 73.9% authors, across all countries
100 * sum(country.summary.by.year.and.journal$count[country.summary.by.year.and.journal$gender %in% c("M", "F")]) / sum(country.summary.by.year.and.journal$count) 

# If I discount all the authors whose forename is initials, then the name genderisation worked for 76.1% of authors, across all countries
100 * sum(country.summary.by.year.and.journal$count[country.summary.by.year.and.journal$gender %in% c("M", "F")]) / sum(country.summary.by.year.and.journal$count[country.summary.by.year.and.journal$gender %in% c("M", "F", "U")]) 

# Here is the genderisation success rate split up by country
success.rate <- country.summary.by.year.and.journal %>% group_by(country) %>% summarise(success.rate = sum(count[gender.95 %in% c("M", "F")]) / sum(count[gender.95 %in% c("M", "F", "U")]), n = sum(count[gender.95 %in% c("M", "F", "U")])) %>% as.data.frame() %>% arrange(success.rate) %>% mutate(country = as.character(country))
success.rate$country <- as.character(sapply(success.rate$country, function(x) {s <- strsplit(x, " ")[[1]]; paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")})) # Fix the capitalisation
success.rate$country <- str_replace(success.rate$country, " And ", " and ")
success.rate$country <- str_replace(success.rate$country, "Usa", "USA")
success.rate <- success.rate[!is.na(success.rate$country) & success.rate$country != "NA", ]
write.csv(success.rate, "../outputs/genderizing.success.rate.by.country.csv", row.names = F)