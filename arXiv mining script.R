setwd("~/Desktop/Gender Pubmed mining project/R scripts")

# library(devtools)
# install_github("ropensci/aRxiv")
library(aRxiv)
library(reshape2)
library(stringr)
library(lme4)
source("Plot and analysis functions.R")

#######################################################################
# FIRST DEFINE A COUPLE OF FUNCTIONS WE WILL NEED LATER
#######################################################################

# Function to search ArXiv when provided with a dataframe of search terms and the number of hits for each of them
# limit argument controls how many searche hits to ask for in each request. Keeps memory usage down, 1000 works fine
search.arxiv <- function(searches, limit = 1000){
  
  # Start of internal function that is used to parse arxiv items into the data we need
  process.arxiv.item <- function(arxiv.item){
    arxiv.item <- arxiv.item[, c("id", "submitted", "authors", "primary_category")] # Discard all but these columns
    arxiv.item$authors <- gsub("[.]", "", arxiv.item$authors) # remove periods from names
    arxiv.item$parsed.authors <- NA
    multi.author <- 1:nrow(arxiv.item) %in% grep("_", arxiv.item$authors) # identify multi author papers
    
    parse.split.names <- function(split.names){
      if(length(split.names)==0) return(NULL)
      chars <- lapply(split.names, function(x) nchar(gsub("[-!?;]", "", x)))
      use.first.one <- sapply(chars, function(y) y[1] > 1) # Use the first listed name if it's >1 character long (excluding punctuation - no "J-" etc)
      use.second.one <- sapply(chars, function(y) y[1] == 1 & y[2] > 1 & length(y) > 2) # else use the second, as in C. Montgomery Burns
      parsed.names <- rep(NA, length(split.names))
      if(length(use.first.one) > 0)  parsed.names[use.first.one] <- unlist(sapply(split.names[use.first.one], function(x) x[1]))
      if(length(use.second.one) > 0) parsed.names[use.second.one] <- unlist(sapply(split.names[use.second.one], function(x) x[2]))
      parsed.names
    }
    
    if(sum(!multi.author) > 0) arxiv.item$parsed.authors[!multi.author] <- parse.split.names(strsplit(arxiv.item$authors[!multi.author], split = " "))
    if(sum(multi.author) > 0) {
      split.author.lists <- melt(strsplit(arxiv.item$authors[multi.author], split = "_"))
      split.author.lists$parsed <- parse.split.names(strsplit(as.character(split.author.lists[,1]), split = " "))
      arxiv.item$parsed.authors[multi.author] <- as.character(with(split.author.lists, tapply(parsed, L1, paste0, collapse = "_")))
    }
    
    arxiv.item <- with(arxiv.item, data.frame(id=id, date=submitted, forenames = parsed.authors, category=primary_category, stringsAsFactors = F))
    arxiv.item <- arxiv.item[!is.na(arxiv.item$forenames), ] # Remove NAs
    arxiv.item <- arxiv.item[str_count(arxiv.item$forenames, "NA")*2 + str_count(arxiv.item$forenames, "_") != nchar(arxiv.item$forenames), ]
    arxiv.item$date <- sapply(strsplit(arxiv.item$date, split = " "), function(x) x[1])
    arxiv.item$date <- gsub("-", "_", arxiv.item$date)
    arxiv.item
  } # end of internal function
  
  searches$output.file <- paste("../data/Arxiv data/search", 1:nrow(searches), ".csv", sep = "")
  done.files <- list.files("../data/Arxiv data", full.names=T)
  searches <- searches[!(searches$output.file %in% done.files), ]
  num.to.do <- nrow(searches)
  
  do.one.search <- function(start){
    search.results <- arxiv_search(searches$search.term[i], batchsize = 1000, limit = limit, sep = "_", start = start)
    if(nrow(search.results) == 0) return(NULL)
    process.arxiv.item(search.results)
  } 
  
  for(i in 1:num.to.do){
    print(paste("Doing file ", i, " out of ", num.to.do, ".", sep=""))
    print(searches$output.file[i]) # print which file we're doing
    starts <- 0; if(searches$n[i] > limit) starts <- seq(0, searches$n[i], limit) # set the start points, allowing retreival of all the data (search limit is 50,000)
    write.csv(do.call("rbind", lapply(starts, do.one.search)), file = searches$output.file[i], row.names = F)
  }
}

# Plot that shows the % female authors in every field on ArXiv
make.summary.plot.arxiv <- function(subset.for.analysis, type, xlim, cutoff.n){
  
  summary.df <- with(subset.for.analysis,
                     as.data.frame(do.call("rbind", tapply(gender, paste(category, sub.category), find.prop.female, return.numbers = T, type = type))))
  names(summary.df) <- c("prop.female", "lowerCI", "upperCI", "nFemales", "nMales")
  cats <- as.data.frame(do.call("rbind", strsplit(sort(unique(apply(subset.for.analysis[,grep("category", names(subset.for.analysis))], 1, paste0, collapse = "_"))), split = "_")))
  names(cats) <- c("category", "sub.category")
  for(i in 1:ncol(cats)) cats[,i] <- as.character(cats[,i])
  summary.df <- cbind(cats, summary.df)
  rownames(summary.df) <- NULL
  rm(cats)
  
  # Merge statistics into mathematics for easier display:
  summary.df$sub.category[summary.df$category == "Statistics"] <- paste("Statistics -", summary.df$sub.category[summary.df$category == "Statistics"])
  summary.df$category[summary.df$category == "Statistics"] <- "Mathematics and Statistics"
  summary.df$category[summary.df$category == "Mathematics"] <- "Mathematics and Statistics"
  
  plot.arxiv.categories <- function(categories, cutoff.n){
    dataset <- summary.df[with(summary.df, nMales + nFemales > cutoff.n & category %in% categories), ]
    dataset <- dataset[order(dataset$prop.female), ] # order by category then prop women authors
    dataset <- dataset[order(dataset$category, dataset$prop.female), ] 
    dataset$sub.category <- as.character(dataset$sub.category)
    dataset$sub.category <- factor(dataset$sub.category, levels = dataset$sub.category)
    ggplot(dataset, aes(x = 100 * prop.female, y = sub.category)) + geom_errorbarh(aes(xmin = 100*lowerCI, xmax = 100*upperCI), height = 0) + geom_point() + facet_wrap(~category, scales = "free_y") + ylab(NULL) + xlab(NULL) + xlim(xlim)
  }
  
  p1 <- plot.arxiv.categories("Computer Science", cutoff.n)
  p2 <- plot.arxiv.categories("Mathematics and Statistics", cutoff.n)
  p3 <- plot.arxiv.categories("Physics", cutoff.n)
  p4 <- plot.arxiv.categories("Astrophysics", cutoff.n)
  p5 <- plot.arxiv.categories("Nonlinear Sciences", cutoff.n)
  p6 <- plot.arxiv.categories("Quantitative Biology", cutoff.n)
  p7 <- plot.arxiv.categories("Quantitative Finance", cutoff.n)
  grid.arrange(arrangeGrob(p1,p2,p3,ncol=3), cbind(rbind(ggplotGrob(p4), ggplotGrob(p5)), rbind(ggplotGrob(p6), ggplotGrob(p7))), bottom = "% female authors \u00B1 95% CIs")
}

# Plot that looks at the change in female authorship over time - for the overall authship gender ratio only
arxiv.year.plot <- function(data, cutoff.n = 50, years = 1995:2016){
  data <- data[with(data, nMales+nFemales > 50 & year %in% years), ]
  ggplot(data, aes(x = year, y = pF.overall*100)) + geom_errorbar(aes(ymin = 100*lCI.overall, ymax = 100*uCI.overall), width = 0) + geom_point() + facet_wrap(~sub.category) + xlab("Year") + ylab("% female authors \u00B1 95% CIs") 
}

# Function to calculate the number of years in the past an ArXiv paper was published (if "before2016")
convert.dates.arxiv <- function(dates, type = "before2016") {
  if(type == "since2000") dates <- as.numeric(as.Date(dates, "%Y_%m_%d")) - as.numeric(as.Date("2000_01_01","%Y_%m_%d"))
  else if(type == "before2016") dates <- as.numeric(as.Date(dates, "%Y_%m_%d")) - as.numeric(as.Date("2016_12_16","%Y_%m_%d"))
  return(dates / 365)  # Convert from days to years
}


#######################################################################
# NOW LET'S CALL UP ARXIV WITH THE 'aRxiv' PACKAGE
#######################################################################

# First let's make a list of the search terms. Big disciplines on Arxiv need to be split up into year ranges,
# because the server won't return a search if there are >50,000 hits
# The ArXiv categories are listed in 'arxiv_cats' from the aRxiv package
big.ones <- c("astro-ph", "gr-qc", "hep-ph", "hep-th", "quant-ph","cat:cond-mat.mes-hall","cat:math.MP","cat:math-ph","cat:cond-mat.stat-mech","cat:cond-mat.mtrl-sci","cat:cond-mat.str-el"," cat:nucl-th","cat:hep-ex","cat:math.AG", "cat:math.AG", "cat:math.CO")
searches <- c(paste("cat:", big.ones, " AND submittedDate:[", 1900," TO ", 2002, "]", sep=""),
                 paste("cat:", big.ones, " AND submittedDate:[", 2003," TO ", 2007, "]", sep=""),
                 paste("cat:", big.ones, " AND submittedDate:[", 2008," TO ", 2012, "]", sep=""),
                 paste("cat:", big.ones, " AND submittedDate:[", 2013," TO ", 2016, "]", sep=""))
searches <- c(searches, paste("cat:", arxiv_cats[,1][!(arxiv_cats[,1] %in% big.ones)], sep=""))
searches <- data.frame(search.term = searches, stringsAsFactors = F)
searches$n <- as.numeric(sapply(searches$search.term, arxiv_count))
searches <- searches[searches$n > 0, ]
searches <- searches[order(searches$n), ]
write.csv(searches, file = "../data/Arxiv data/searches.csv", row.names = F)

#######################################################################
# NOW RUN THE SEARCH AND WRITE RESULTS TO DISK - TAKES A WHILE
#######################################################################
searches <- read.csv("../data/Arxiv data/searches.csv", stringsAsFactors = F)
search.arxiv(searches)
arxiv.data <- do.call("rbind", lapply(list.files("../data/Arxiv data", full.names = T)[grep("search[[:digit:]]", list.files("../data/Arxiv data", full.names = T))], read.csv, stringsAsFactors = F))
arxiv.data <- arxiv.data[!duplicated(arxiv.data$id), ] # remove duplicates. There are 579999 unique papers

#######################################################################
# NOW ADD GENDER INFORMATION TO THE DATA FROM ARXIV
#######################################################################

genderize.data <- read.csv("../data/genderize.names.mastersheet.csv", stringsAsFactors = F)
genderize.data <- genderize.data[genderize.data$gender != "U", ] %>% arrange(name) 
genderize.data <- genderize.data[unlist(tapply(-genderize.data$count, genderize.data$name, order)) == 1, ] # restrict to the instance of each name that has the most data (we won't use the country info)

all.names <- strsplit(arxiv.data$forenames, split = "_") # Remove the underscores to get a character vector of all the author's forenames
paperID <- unlist(mapply(rep, 1:length(all.names), each = sapply(all.names, length))) # paper mappings for all the authors
all.names <- unlist(all.names)
all.names[all.names %in% c("van", "Van")] <- "dkufghsidhgoisdfhgiu" # There are lots of people where my code thinks their first name is "van" or "Van", but this is more likely part of the last name. Score the gender of these people as "unknown" to be safe.
genders <- genderize.data$gender[match(tolower(all.names), genderize.data$name)] 
genders[is.na(genders)] <- "U"
genders[str_detect(all.names, "^[[:upper:]]+$")] <- "U" # to guard against initials that look like names (e.g. AMY Smith), remove genders associated with names in ALL CAPS
arxiv.data$gender <- tapply(genders, paperID, function(x) paste0(x, collapse = ""))
scores <- genderize.data$prob[match(tolower(all.names), genderize.data$name)]
scores[genders == "M"] <- 1 - scores[genders == "M"]
arxiv.data$gender.score <- tapply(scores, paperID, function(x) paste0(x, collapse = "_"))

# Success of gender assignment for ArXiv:

genders[str_detect(all.names, "^[[:upper:]]+$")] <- "I" # Change the names in ALL CAPS to I, as they are mostly initials
genders[nchar(all.names) < 4 & str_detect(all.names, "[~-]")] <- "I" # Names with fewer than 3 real letters are generally not real names, e.g. A-Y Jones

sum(str_detect(genders, "[MF]")) / sum(str_detect(genders, "[MFU]")) # 85.0% gendering success rate, when we exclude the putative initials
sum(str_detect(genders, "[MF]")) / sum(str_detect(genders, "[MFUI]")) # 57.7% gendering success rate, when we include the putative initials

#######################################################################
# NOW CLEAN THE 'CATEGORY' FIELD - IT'S NOT NEATLY-CURATED 
#######################################################################

lukes.cats <- arxiv_cats
lukes.cats$real.abbrev <- overall$category[match(lukes.cats$abbreviation, overall$category)]
lukes.cats <- rbind(lukes.cats[,2:3], data.frame(description = NA, real.abbrev = overall$category[!( overall$category %in% arxiv_cats$real.abbrev)]))
abbrevs <- lukes.cats$real.abbrev[is.na(lukes.cats$description)]
insertions <- c("Physics - Plasma Physics", "Statistics - Theory", "Nonlinear Sciences - Adaptation and Self-Organizing Systems", "Physics - Accelerator Physics", 
                "Physics - Atomic Physics", "Physics - Chemical Physics", "Physics - Superconductivity", "Quantitative Finance - Economics", "Nonlinear Sciences - Cellular Automata and Lattice Gases",
                "Statistics - Other", "Nonlinear Sciences - Exactly Solvable and Integrable Systems", "Quantitative Finance - Mathematical Finance", "Physics - Materials Science", 
                "Quantitative Finance - Computational Finance", "Quantitative Finance - Trading and Market Microstructure", "Quantitative Finance - Risk Management", 
                "Mathematics - Functional Analysis", "Quantitative Finance - Portfolio Management", "Nonlinear Sciences - Adaptation and Self-Organizing Systems", 
                "Quantitative Finance - Pricing of Securities", "Computer Science - Emerging Technologies", "Astrophysics - Astrophysics of Galaxies", "Quantitative Finance - General Finance", 
                "Quantitative Finance - Statistical Finance", "Mathematics - Differential Geometry", "Nonlinear Sciences - Pattern Formation and Solitons", "Mathematics - Quantum Algebra", 
                "Computer Science - Formal Languages and Automata Theory", "Condensed Matter - Other", "Mathematics - Algebraic Geometry", "Astrophysics - Earth and Planetary Astrophysics", 
                "Nonlinear Sciences - Chaotic Dynamics", "Computer Science - Computation and Language", "Astrophysics - Solar and Stellar Astrophysics",
                "Computer Science - Systems and Control", "Astrophysics - High Energy Astrophysical Phenomena", "Astrophysics - Instrumentation and Methods for Astrophysics", 
                "Astrophysics - Cosmology and Nongalactic Astrophysics", "Condensed Matter - Quantum Gases", "Computer Science - Social and Information Networks")

lukes.cats$description[lukes.cats$real.abbrev %in% abbrevs] <- insertions
lukes.cats <- lukes.cats[!(is.na(lukes.cats$real.abbrev)), ] 
arxiv.data$long.cat <- lukes.cats$description[match(arxiv.data$category, lukes.cats$real.abbrev)]

# Sort out the category naming system
names(arxiv.data)[names(arxiv.data) == "category"] <- "cat.abbreviation"
arxiv.data$category <- sapply(strsplit(as.character(arxiv.data$long.cat), split = " - "), head, 1)
arxiv.data$sub.category <- sapply(strsplit(as.character(arxiv.data$long.cat), split = " - "), function(x) x[2])
arxiv.data$category[arxiv.data$category == "General Relativity and Quantum Cosmology"] <- "Physics"
arxiv.data$sub.category[arxiv.data$category == "General Relativity and Quantum Cosmology"] <- NA
arxiv.data$category[arxiv.data$category == "High Energy Physics"] <- "Physics"
arxiv.data$sub.category[arxiv.data$category == "High Energy Physics"] <- NA
arxiv.data$category[arxiv.data$category == "Quantum Physics"] <- "Physics"
arxiv.data$sub.category[arxiv.data$category == "Quantum Physics"] <- NA
arxiv.data$category[arxiv.data$category == "Mathematical Physics"] <- "Physics"
arxiv.data$sub.category[arxiv.data$category == "Mathematical Physics"] <- NA
arxiv.data$category[arxiv.data$category == "Nuclear Theory"] <- "Physics"
arxiv.data$sub.category[arxiv.data$category == "Nuclear Theory"] <- NA
arxiv.data$category[arxiv.data$category == "Nuclear Experiment"] <- "Physics"
arxiv.data$sub.category[arxiv.data$category == "Nuclear Experiment"] <- NA
arxiv.data$category[arxiv.data$sub.category == "Quantum Gases"] <- "Physics"      # merged the condensed matter list into physics
arxiv.data$category[arxiv.data$category == "Condensed Matter"] <- "Physics"
arxiv.data$sub.category[is.na(arxiv.data$sub.category)] <- "Unspecified"
arxiv.data$sub.category[arxiv.data$sub.category == "Statistics"] <- "Unspecified"
arxiv.data <- arxiv.data[, !(names(arxiv.data) %in% "long.cat")]


# Lastly, add a gender.95 column, which classifies as "U" all those people whose gender is not known with >= 95% certainty (uses the scores from genderize.io)
split.gender <- strsplit(arxiv.data$gender, split = "")
split.gender <- data.frame(gender = unlist(split.gender), paper = unlist(mapply(rep, 1:length(split.gender), sapply(split.gender, length))))
split.gender$score <- unlist(strsplit(arxiv.data$gender.score, split = "_"))
split.gender$score[split.gender$score == "NA"] <- NA
split.gender$score <- as.numeric(split.gender$score)
split.gender$gender.95 <- split.gender$gender
split.gender$gender.95[split.gender$score > 0.05 & split.gender$score < 0.95] <- "U"
split.gender$score[split.gender$gender == "U"] <- NA # Some of the unknown genders are erroneously showing a score of e.g. 1. Fix that while we are here
arxiv.data$gender.95 <- (split.gender %>% group_by(paper) %>% summarise(gender.95 = paste0(gender.95, collapse = "")) %>% as.data.frame())[,2] # Add the gender.95 and fixed score columns back to the dataframe
arxiv.data$gender.score <- (split.gender %>% group_by(paper) %>% summarise(gender.score = paste0(gender.score, collapse = "_")) %>% as.data.frame())[,2]
arxiv.data <- arxiv.data %>% select(id, date, forenames, cat.abbreviation, gender, gender.95, gender.score, category, sub.category) # reorder columns for neatness
rm(split.gender)

# Write the final dataset to disk
# write.csv(arxiv.data, file = "../data for analysis/arxiv data.csv", row.names = F)


#######################################################################
#  MAKE A USEFUL SUMMARY OF THE ARXIV DATA 
# (% women per discipline, per year, for all 4 authorship types - with confidence intervals)
#######################################################################

arxiv.data <- read.csv("../data for analysis/arxiv data.csv", stringsAsFactors = F)

arxiv.data$year <- as.numeric(substr(arxiv.data$date, 1, 4)) # add a year column

pasted.categories <- paste(arxiv.data$category, arxiv.data$sub.category, sep = "_")
year.summary.subcategories <- melt(lapply(tapply(arxiv.data$year, pasted.categories, unique), melt))[,2:3] %>% 
  mutate(nFemales.first = NA, nMales.first = NA, pF.first = NA, lCI.first = NA, uCI.first = NA, 
         nFemales.last = NA, nMales.last = NA, pF.last = NA, lCI.last = NA, uCI.last = NA, 
         nFemales.single = NA, nMales.single = NA, pF.single = NA, lCI.single = NA, uCI.single = NA, 
         nFemales.overall = NA, nMales.overall = NA, pF.overall = NA, lCI.overall = NA, uCI.overall = NA) %>% rename(year = value, category = L1)
for(i in 1:nrow(year.summary.subcategories)) {
  focal <- arxiv.data$gender.95[arxiv.data$year == year.summary.subcategories$year[i] & pasted.categories == year.summary.subcategories$category[i]] # for the focal combo of year and sub-category...
  year.summary.subcategories[i, grep(".first", names(year.summary.subcategories))] <- find.prop.female(focal, type = "first_but_not_single", return.numbers = T)[c(4:5,1:3)] # get the number of M+F authors
  year.summary.subcategories[i, grep(".last", names(year.summary.subcategories))] <- find.prop.female(focal, type = "last_but_not_single", return.numbers = T)[c(4:5,1:3)] # and the 95% CIs on % females
  year.summary.subcategories[i, grep(".single", names(year.summary.subcategories))] <- find.prop.female(focal, type = "single", return.numbers = T)[c(4:5,1:3)]
  year.summary.subcategories[i, grep(".overall", names(year.summary.subcategories))] <- find.prop.female(focal, type = "overall", return.numbers = T)[c(4:5,1:3)]
}
year.summary.subcategories <- year.summary.subcategories[with(year.summary.subcategories, !is.na(pF.overall) | !is.na(nMales.overall) | !is.na(nFemales.overall)), ] # get rid of combinations with no data at all
for(i in 1:4) {  # Replace the NAs with zeros
  focal <- year.summary.subcategories[, grep("nFemales", names(year.summary.subcategories))[i]] 
  focal[is.na(focal)] <- 0
  year.summary.subcategories[, grep("nFemales", names(year.summary.subcategories))[i]]  <- focal
  focal <- year.summary.subcategories[, grep("nMales", names(year.summary.subcategories))[i]] 
  focal[is.na(focal)] <- 0
  year.summary.subcategories[, grep("nMales", names(year.summary.subcategories))[i]]  <- focal
}

year.summary.subcategories$pF.first[with(year.summary.subcategories, nFemales.first + nMales.first == 0)] <- NA   # Make sure the CIs are nice and sensible
year.summary.subcategories$lCI.first[with(year.summary.subcategories, nFemales.first + nMales.first == 0)] <- NA
year.summary.subcategories$uCI.first[with(year.summary.subcategories, nFemales.first + nMales.first == 0)] <- NA
year.summary.subcategories$pF.last[with(year.summary.subcategories, nFemales.last + nMales.last == 0)] <- NA
year.summary.subcategories$lCI.last[with(year.summary.subcategories, nFemales.last + nMales.last == 0)] <- NA
year.summary.subcategories$uCI.last[with(year.summary.subcategories, nFemales.last + nMales.last == 0)] <- NA
year.summary.subcategories$pF.single[with(year.summary.subcategories, nFemales.single + nMales.single == 0)] <- NA
year.summary.subcategories$lCI.single[with(year.summary.subcategories, nFemales.single + nMales.single == 0)] <- NA
year.summary.subcategories$uCI.single[with(year.summary.subcategories, nFemales.single + nMales.single == 0)] <- NA
year.summary.subcategories <- data.frame(do.call("rbind", strsplit(year.summary.subcategories$category, split = "_")), year.summary.subcategories[,-2]) %>% rename(category = X1, sub.category = X2)
year.summary.subcategories <- year.summary.subcategories %>% arrange(category, sub.category, year)

rm(pasted.categories)

year.summary.categories <- year.summary.subcategories %>% group_by(category, year) %>% summarise(nFemales.first = sum(nFemales.first), nMales.first = sum(nMales.first),  pF.first = NA, lCI.first = NA, uCI.first = NA, nFemales.last = sum(nFemales.last), nMales.last = sum(nMales.last), pF.last=NA, lCI.last=NA, uCI.last=NA, nFemales.single = sum(nFemales.single), nMales.single = sum(nMales.single), pF.single=NA, lCI.single=NA, uCI.single=NA, nFemales.overall = sum(nFemales.overall), nMales.overall = sum(nMales.overall), pF.overall=NA, lCI.overall=NA, uCI.overall=NA) %>% as.data.frame()
year.summary.categories$pF.first[with(year.summary.categories, nFemales.first + nMales.first > 0)] <- with(year.summary.categories[year.summary.categories$nFemales.first + year.summary.categories$nMales.first > 0, ], nFemales.first/(nFemales.first + nMales.first))
year.summary.categories$pF.last[with(year.summary.categories, nFemales.last+ nMales.last > 0)] <- with(year.summary.categories[year.summary.categories$nFemales.last + year.summary.categories$nMales.last > 0, ], nFemales.last/(nFemales.last + nMales.last))
year.summary.categories$pF.single[with(year.summary.categories, nFemales.single + nMales.single > 0)] <- with(year.summary.categories[year.summary.categories$nFemales.single + year.summary.categories$nMales.single > 0, ], nFemales.single/(nFemales.single + nMales.single))
year.summary.categories$pF.overall[with(year.summary.categories, nFemales.overall + nMales.overall > 0)] <- with(year.summary.categories[year.summary.categories$nFemales.overall + year.summary.categories$nMales.overall > 0, ], nFemales.overall/(nFemales.overall + nMales.overall))

for(i in 1:nrow(year.summary.categories)){
  try(year.summary.categories[i, grep("CI.first", names(year.summary.categories))] <- with(year.summary.categories[i, ], binom.test(nFemales.first, nFemales.first+nMales.first)$conf.int), silent = T)
  try(year.summary.categories[i, grep("CI.last", names(year.summary.categories))] <- with(year.summary.categories[i, ], binom.test(nFemales.last, nFemales.last+nMales.last)$conf.int), silent = T)
  try(year.summary.categories[i, grep("CI.single", names(year.summary.categories))] <- with(year.summary.categories[i, ], binom.test(nFemales.single,  nFemales.single + nMales.single)$conf.int), silent = T)
  try(year.summary.categories[i, grep("CI.overall", names(year.summary.categories))] <- with(year.summary.categories[i, ], binom.test(nFemales.overall, nFemales.overall + nMales.overall)$conf.int), silent = T)
}

# Make some summarises for major categories in the years 2012-2016 inclusive, for plotting alongside the PubMed data
year.summary.categories.since2012 <- year.summary.categories %>% filter(year >= 2012) %>% group_by(category) %>% summarise(nFemales.first = sum(nFemales.first), nMales.first = sum(nMales.first),  pF.first = NA, lCI.first = NA, uCI.first = NA, nFemales.last = sum(nFemales.last), nMales.last = sum(nMales.last), pF.last=NA, lCI.last=NA, uCI.last=NA, nFemales.single = sum(nFemales.single), nMales.single = sum(nMales.single), pF.single=NA, lCI.single=NA, uCI.single=NA, nFemales.overall = sum(nFemales.overall), nMales.overall = sum(nMales.overall), pF.overall=NA, lCI.overall=NA, uCI.overall=NA) %>% as.data.frame()

year.summary.categories.since2012$pF.first <- with(year.summary.categories.since2012, nFemales.first / (nFemales.first + nMales.first))
year.summary.categories.since2012$pF.last <- with(year.summary.categories.since2012, nFemales.last / (nFemales.last + nMales.last))
year.summary.categories.since2012$pF.single <- with(year.summary.categories.since2012, nFemales.single / (nFemales.single + nMales.single))
year.summary.categories.since2012$pF.overall <- with(year.summary.categories.since2012, nFemales.overall / (nFemales.overall + nMales.overall))

for(i in 1:nrow(year.summary.categories.since2012)){
  try(year.summary.categories.since2012[i, grep("CI.first", names(year.summary.categories.since2012))] <- with(year.summary.categories.since2012[i, ], binom.test(nFemales.first, nFemales.first+nMales.first)$conf.int), silent = T)
  try(year.summary.categories.since2012[i, grep("CI.last", names(year.summary.categories.since2012))] <- with(year.summary.categories.since2012[i, ], binom.test(nFemales.last, nFemales.last+nMales.last)$conf.int), silent = T)
  try(year.summary.categories.since2012[i, grep("CI.single", names(year.summary.categories.since2012))] <- with(year.summary.categories.since2012[i, ], binom.test(nFemales.single,  nFemales.single + nMales.single)$conf.int), silent = T)
  try(year.summary.categories.since2012[i, grep("CI.overall", names(year.summary.categories.since2012))] <- with(year.summary.categories.since2012[i, ], binom.test(nFemales.overall, nFemales.overall + nMales.overall)$conf.int), silent = T)
}
# Reorder the columns same as PubMed data
year.summary.categories.since2012 <- year.summary.categories.since2012 %>% select(category, nFemales.overall, nMales.overall, nFemales.first, nMales.first, nFemales.last, nMales.last, nFemales.single, nMales.single, lCI.overall, uCI.overall, lCI.first, uCI.first,  lCI.last,  uCI.last, lCI.single, uCI.single)



# Make some summarised for sub-categories in the years 2012-2016 inclusive, for plotting in the supplementary material
year.summary.subcategories.since2012 <- year.summary.subcategories %>% filter(year >= 2012) %>% mutate(sub.category = paste(category, sub.category, sep = "_")) %>% group_by(sub.category) %>% summarise(nFemales.first = sum(nFemales.first), nMales.first = sum(nMales.first),  pF.first = NA, lCI.first = NA, uCI.first = NA, nFemales.last = sum(nFemales.last), nMales.last = sum(nMales.last), pF.last=NA, lCI.last=NA, uCI.last=NA, nFemales.single = sum(nFemales.single), nMales.single = sum(nMales.single), pF.single=NA, lCI.single=NA, uCI.single=NA, nFemales.overall = sum(nFemales.overall), nMales.overall = sum(nMales.overall), pF.overall=NA, lCI.overall=NA, uCI.overall=NA) %>% as.data.frame()

year.summary.subcategories.since2012$pF.first <- with(year.summary.subcategories.since2012, nFemales.first / (nFemales.first + nMales.first))
year.summary.subcategories.since2012$pF.last <- with(year.summary.subcategories.since2012, nFemales.last / (nFemales.last + nMales.last))
year.summary.subcategories.since2012$pF.single <- with(year.summary.subcategories.since2012, nFemales.single / (nFemales.single + nMales.single))
year.summary.subcategories.since2012$pF.overall <- with(year.summary.subcategories.since2012, nFemales.overall / (nFemales.overall + nMales.overall))

for(i in 1:nrow(year.summary.subcategories.since2012)){
  try(year.summary.subcategories.since2012[i, grep("CI.first", names(year.summary.subcategories.since2012))] <- with(year.summary.subcategories.since2012[i, ], binom.test(nFemales.first, nFemales.first+nMales.first)$conf.int), silent = T)
  try(year.summary.subcategories.since2012[i, grep("CI.last", names(year.summary.subcategories.since2012))] <- with(year.summary.subcategories.since2012[i, ], binom.test(nFemales.last, nFemales.last+nMales.last)$conf.int), silent = T)
  try(year.summary.subcategories.since2012[i, grep("CI.single", names(year.summary.subcategories.since2012))] <- with(year.summary.subcategories.since2012[i, ], binom.test(nFemales.single,  nFemales.single + nMales.single)$conf.int), silent = T)
  try(year.summary.subcategories.since2012[i, grep("CI.overall", names(year.summary.subcategories.since2012))] <- with(year.summary.subcategories.since2012[i, ], binom.test(nFemales.overall, nFemales.overall + nMales.overall)$conf.int), silent = T)
}
# Reorder the columns same as PubMed data
year.summary.subcategories.since2012 <- year.summary.subcategories.since2012 %>% select(sub.category, nFemales.overall, nMales.overall, nFemales.first, nMales.first, nFemales.last, nMales.last, nFemales.single, nMales.single, lCI.overall, uCI.overall, lCI.first, uCI.first,  lCI.last,  uCI.last, lCI.single, uCI.single)


# save summaries to disk
write.csv(year.summary.subcategories, "../data for analysis/arxiv data summarised by year and subcategory.csv", row.names = F)
write.csv(year.summary.categories, "../data for analysis/arxiv data summarised by year and category.csv", row.names = F)
write.csv(year.summary.categories.since2012, "../data for analysis/arxiv data since 2012 summarised by category.csv", row.names = F)
write.csv(year.summary.subcategories.since2012, "../data for analysis/arxiv data since 2012 summarised by sub category.csv", row.names = F)


#######################################################################
# DESCRIPTIVE STATISTICS OF THE ARXIV DATA
#######################################################################

# Load up the ArXiv data
arxiv.data <- read.csv("../data for analysis/arxiv data.csv", stringsAsFactors = F)
nrow(arxiv.data) # 579999 preprint in total
arxiv.data <- arxiv.data[grep("[MF]", arxiv.data$gender), ] # Discard all papers where author gender is unknown for all authors
nrow(arxiv.data) # 538688 preprints where at least one author's gender is known
sum(str_count(arxiv.data$gender, "[MF]")) # 1183805 authors with at least partially known gender 
sum(str_count(arxiv.data$gender.95, "[MF]")) # 1043330 authors with gender known to >= 95% certainty
sum(nchar(arxiv.data$gender)) # 1938039 authors in total

# Load up the summary datasets
year.summary.subcategories <- read.csv("../data for analysis/arxiv data summarised by year and subcategory.csv", stringsAsFactors = F)
year.summary.categories <- read.csv("../data for analysis/arxiv data summarised by year and category.csv", stringsAsFactors = F)

# Sample size per sub-category histogram
data.frame(table(paste(arxiv.data$category, arxiv.data$sub.category))) %>% ggplot(aes(x = Freq)) + geom_histogram(fill = "skyblue", colour = "black", bins = 30) + xlab("Number of articles in sub-category") + ylab("Frequency") + scale_x_log10(breaks = c(10,100,1000,10000)) 




#######################################################################
# PLOTTING THE ARXIV DATA
#######################################################################

# Let's plot the % women in every sub-category in the years 2012 - 2016, for different author positions
subset.for.analysis <- arxiv.data[as.numeric(substr(arxiv.data$date, 1, 4)) >= 2012, ]
make.summary.plot.arxiv(subset.for.analysis, "overall", xlim = c(0,49), cutoff.n = 10)
make.summary.plot.arxiv(subset.for.analysis, "first_but_not_single", xlim = c(0,46), cutoff.n = 10)
make.summary.plot.arxiv(subset.for.analysis, "last_but_not_single", xlim = c(0,40), cutoff.n = 10)
make.summary.plot.arxiv(subset.for.analysis, "single", xlim = c(0,40), cutoff.n = 10)

# Plot the change in % women with time in each major category except Quantitative Finance (not much data there)
years.to.plot <- seq(1995, 2030,by=5)
year.summary.categories %>% filter(nMales+nFemales > 200 ) %>% filter(category != "Quantitative Finance") %>%
  ggplot(aes(x = year, y = 100*pF.overall, weight = nMales+nFemales)) + stat_smooth(method = "lm", fullrange = T) + geom_errorbar(aes(ymin = 100*lCI.overall, ymax = 100*uCI.overall), width = 0) + geom_point() + facet_wrap(~category) + scale_x_continuous(breaks = years.to.plot, limits = range(years.to.plot)) + ylab("% female authors \u00B1 95% CIs")



# Plot the change in % women with time in each sub-category
arxiv.year.plot(year.summary.subcategories[year.summary.subcategories$category %in% c("Physics", "Astrophysics"), ])
arxiv.year.plot(year.summary.subcategories[year.summary.subcategories$category == "Computer Science", ])
arxiv.year.plot(year.summary.subcategories[year.summary.subcategories$category %in% c("Mathematics") , ])
arxiv.year.plot(year.summary.subcategories[year.summary.subcategories$category == "Quantitative Biology", ])
arxiv.year.plot(year.summary.subcategories[year.summary.subcategories$category == "Quantitative Finance", ])
arxiv.year.plot(year.summary.subcategories[year.summary.subcategories$category == "Nonlinear Sciences", ])


# Try to plot all 4 authorship positions - not a pretty plot, too few data
arxiv.year.plot.multi <- function(year.summary.categories, cutoff.n = 100, years = 1995:2016){
  cols <- c(brewer.pal(3, "Set1"), "black")
  data <- year.summary.categories[with(year.summary.categories, nMales+nFemales > cutoff.n & year %in% years), ]
  ggplot(data, aes(x = year, y = pF.overall*100)) + 
    geom_line(colour = cols[4]) + 
    geom_line(aes(y = pF.first*100), colour = cols[1]) + 
    geom_line(aes(y = pF.last*100), colour = cols[2]) + 
    geom_line(aes(y = pF.single*100), colour = cols[3]) + 
    facet_wrap(~category) + 
    ylab("% female authors \u00B1 95% CIs") 
}
arxiv.year.plot.multi(year.summary.categories)




#######################################################################
# PREDICTIVE MODELS FOR THE ARXIV DATA
#######################################################################



# Fit a binomial GLMM with the first-listed author's gender as 0/1 (includes single and multi-author papers)
# Random effects: sub-category and the year:sub.category interaction
# Fixed effect: year (note 'year' is time before the present, to the nearest day, measured in years) 
# Fit one model per major category (e.g. physics, maths, CompSci...)
temp.arxiv.data <-  arxiv.data %>%
  mutate(category = replace(category, category %in% c("Mathematics", "Statistics"), "Mathematics and Statistics")) %>% # First merge a couple of category pairs
  mutate(category = replace(category, category %in% c("Physics", "Astrophysics"), "Physics and Astrophysics"))
# Calculate the number of years before 2016_12_16 each paper was published on ArXiv (that's that date that I text-mined the data)
temp.arxiv.data$years.ago <- convert.dates.arxiv(temp.arxiv.data$date)
temp.arxiv.data$first.author.female <- NA
temp.arxiv.data$first.author.female[substr(temp.arxiv.data$gender, 1, 1) == "M"] <- 0
temp.arxiv.data$first.author.female[substr(temp.arxiv.data$gender, 1, 1) == "F"] <- 1
temp.arxiv.data <- temp.arxiv.data[!is.na(temp.arxiv.data$first.author.female), ] # Discard papers where first listed author's gender is not known

# Do a model within each major category
models <- lapply(sort(unique(temp.arxiv.data$category)), function(x) glmer(first.author.female ~ years.ago + (years.ago|sub.category), data = temp.arxiv.data[temp.arxiv.data$category == x, ], family = "binomial", control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e5))))
summaries <- lapply(models, summary)
names(summaries) <- sort(unique(arxiv.data$category))
summaries

# Model all the data, treating major categories as the random effect
summary(glmer(first.author.female ~ years.ago + (years.ago|category), family = "binomial", data = temp.arxiv.data, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun=2e6))))

