# This function is used to call Pubmed from R. In this project, I just use it to look at the Abstracts as a quality control step.

# If clean.up is FALSE, it just returns the entire PubMed record in a big, confusing list of lists. Useful for bug-checks etc
# If TRUE, it extracts only the data that we want and returns it as a nice data frame
# If save.abstracts=TRUE, we take the entire abstract text (or the first-listed chunk of text, for abstracts split over multiple fields)
# Else, we jsut count the number of characters in the abstract
pubmed_fetch <- function(ids, clean.up = T, save.abstracts = F, save.titles = F){
  
  require(RCurl)
  require(XML)
  
  # Internal function that gets all the data that we want out of a pubmed record
  retrieve.data <- function(record, save.abstracts, save.titles){
    
    forenames <- get.first.names(record)   # Get the first names using the get.first.names function
    record <- unlist(record)  # Now unlist the record
    field <- names(record) # Store the names of each Pubmed field
    record <- as.character(record) # Ensure the record is a character not a factor (needed for string matching)
    
    # Prepare the output data frame. It will contain NAs unless we can find the data to fill each field in the PubMed record
    output <- data.frame(id=NA, doi=NA, date=NA, journal = NA, type = NA, abstract = NA)
    
    id <- record[grep("PMID.text", field)] # Try to get the PubMed ID
    if(length(id) >= 1) output$id <- id[1] # If there is one, add it to the output
    
    # Get the date, and store it in this format: 15_01_2012
    date <- paste0(c(record[field == "MedlineCitation.DateCreated.Day"], record[field == "MedlineCitation.DateCreated.Month"], record[field == "MedlineCitation.DateCreated.Year"]), collapse = "_")
    if(length(date) >= 1) output$date <- date[1]
    
    journal <- record[field == "MedlineCitation.Article.Journal.ISOAbbreviation"] # Get the journal
    if(length(journal) == 1) output$journal <- journal
    
    # Get the publication type list
    type <- record[grep("PublicationTypeList.PublicationType", field)]
    if(length(type) == 1) output$type <- type
    else if(length(type) > 1) output$type <- paste0(type, collapse = "_") # If there are multiples, collapse to a single character string with underscores
    
    abstract <- record[grep("Abstract.AbstractText", field)] # Get the entire abstract(s)
    if(length(abstract) >= 1 & !(save.abstracts)) output$abstract <- max(nchar(abstract, keepNA = T)) # Save the length of the abstract (or the longest abstract chunk, for papers with multi-part abstracts split over multiple fields). Discarding the abstracts makes the files much smaller (MB rather than GB)
    else if (length(abstract) >= 1 & save.abstracts) output$abstract <- paste0(abstract, collapse = "___")
    
    if(save.titles) output$title <- record[grep("Article.ArticleTitle", field)] 
    
    # Searching for the DOI is slightly hard, because annoyingly there is no obvious, consistent 'doi' field in PubMed!
    # Here, I define a DOI is defined as a string that contains "10." followed by at least one digit, then a slash, then some letters and/or numbers. 
    to.search.for.doi <- record[c(grep("ArticleIdList", field), grep("ELocationID", field))] # Search only the "ArticleIdList" and "ELocationID" fields
    doi <- unique(to.search.for.doi[grep("10[.][[:digit:]]+/[[:alnum:]]", to.search.for.doi)]) 
    if(length(doi) >= 1) output$doi <- doi[1] # If there are multiple possible DOIs, just take the first candidate DOI in the "ArticleIdList" and "ELocationID" fields
    
    return(data.frame(output[,!(names(output) %in% "abstract")], forenames, abstract = output$abstract, stringsAsFactors = F))
  }
  
  
  
  args <- c(id = paste(ids, collapse = ","), db = "pubmed", rettype = "xml", email = entrez_email, tool = "rpubmed")
  
  url_args <- paste(paste(names(args), args, sep="="), collapse = "&")
  base_url <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?retmode=full"
  url_string <- paste(base_url, url_args, sep = "&")
  records <- getURL(url_string)
  Sys.sleep(0.33)     #NCBI limits requests to three per second
  records <- xmlToList(xmlTreeParse(records, useInternalNodes = TRUE))
  if(!clean.up) return(records)  
  
  output <- do.call("rbind", lapply(records, retrieve.data, save.abstracts=save.abstracts, save.titles=save.titles))
  rm(records) # Try to minimise pressure on memory by deleting this as early as we can
  rownames(output) <- NULL
  return(output)
}




# If type = "since2000", we represent the date at which the paper appeared on Pubmed as number of days since the 1st January 2000, divided by 365
# Else if type = "beforePresent", it's the number of days before 20/08/2016, divided by 365 - this is the date that I downloaded a local copy of the Medline database.
convert.dates <- function(dates, type = "before2016") {
  if(type == "since2000") dates <- as.numeric(as.Date(dates, "%d_%m_%Y")) - as.numeric(as.Date("01_01_2000","%d_%m_%Y"))
  else if(type == "before2016") dates <- as.numeric(as.Date(dates, "%d_%m_%Y")) - as.numeric(as.Date("20_08_2016","%d_%m_%Y"))
  return(dates / 365.25)  # Convert from days to years
}


# Function to convert a log odds ratio to a probability
# toProb <- function(x) exp(as.numeric(x)) / (1 + exp(as.numeric(x))) 



# Properly capitalise country names
capitalise.countries <- function(countries){
  countries <- as.character(sapply(as.character(countries), function(x) {s <- strsplit(x, " ")[[1]]; paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")})) 
  countries <- str_replace(countries, " And ", " and ")
  str_replace(countries, "Usa", "USA")
}

# Used to predict the gender ratio (as proportion), given a date t, and a value for the two shape parameters r and c
pfunc <- function(t, r, c){exp(0.5*r*t)/(2*exp(0.5*r*t)+c)} # function for a sigmoid curve that converges to 0 or 0.5 at +ve infinity, and 0.5 or 0 at -ve infinity

# Get the estimated gender ratio, its rate of change, and the estimated number of years until gender parity (defined as 45-55% female)
gender.stats <- function(item, filter.type, position, country = "all", data.source, chunk.size = 10, nChunks = 100, plot = T, print.each.one = F, run.sim = F, export.JSON = F, verbose = T){
  
  # Print which item we are doing
  if(verbose) print(paste("Doing ", item, " (", position, ")", sep = ""))
  
  ############ DEFINE INTERNAL FUNCTIONS 
  pfunc.deriv <- function(p, r) r*p*(0.5-p) # the first deriviative of the pfunc function (needed to find the rate of change in gender ratio on a specified date)
  
  # log likelihood of data given the parameters, assuming binomial distribution
  # returns the negative log-likelihood to be minimised with optim()
  find.ll<-function(data, par) suppressWarnings(with(data, -1 * sum(dbinom(x = nFemales, size = n, prob = pfunc(t = date, r = par[1], c = par[2]), log = TRUE))))
  
  # Internal function to perform 1D optimisation and find the parameters r and c that optimises the log likelihood for the focal set of data. Arbitrarily chosen starting values
  run.optimiser <- function(data) optim(par = c(0.1, 1), find.ll, data = data)
  
  # Used to find the 95% CIs on a proportion (uses default method of the binom.test function). Works for a pair of columns nFemales and nMales
  get.CIs <- function(nFemales.column, nMales.column){
    CIs <- data.frame(lower = numeric(length(nFemales.column)), upper = NA)
    for(i in 1:length(nFemales.column)){
      CIs.foc <- as.numeric(binom.test(nFemales.column[i], nFemales.column[i] + nMales.column[i])$conf.int)
      CIs$lowerCI[i] <- CIs.foc[1]
      CIs$upperCI[i] <- CIs.foc[2]
    }
    return(100 * CIs)
  }
  
  # Function to export the results to disk in JSON format
  run.JSON.export <- function(data.for.web.app=NULL){  # , return.file.name.only = F
    if(position == "first") position <- "First"
    if(position == "last") position <- "Last"
    if(position == "single") position <- "Single"
    if(position == "overall") position <- "Overall"
    country <- capitalise.countries(country)
    
    if(filter.type == "disc" & country == "All"){
      data.for.web.app <- list(Discipline = item, Journal = "allJournals", Country = "allCountries", Position = position, Points = data.for.web.app[[1]], Curve = data.for.web.app[[2]])
      file.name <- paste("allJournals", item, position, "allCountries", sep = "_")
    }
    
    if(filter.type == "disc" & country != "All"){
      data.for.web.app <- list(Discipline = item, Journal = "allJournals", Country = country, Position = position, Points = data.for.web.app[[1]], Curve = data.for.web.app[[2]])
      file.name <- paste("allJournals", item, position, country, sep = "_")
    }
    
    if(filter.type == "journal" & country == "All"){
      data.for.web.app <- list(Discipline = discipline, Journal = item, Country = "allCountries", Position = position, Points = data.for.web.app[[1]], Curve = data.for.web.app[[2]])
      file.name <- paste(item, discipline, position, "allCountries", sep = "_")
    }
    
    if(filter.type == "journal" & country != "All"){
      data.for.web.app <- list(Discipline = discipline, Journal = item, Country = country, Position = position, Points = data.for.web.app[[1]], Curve = data.for.web.app[[2]])
      file.name <- paste(item, discipline, position, country, sep = "_")
    }
    
    if(filter.type == "everything.for.one.country"){
      data.for.web.app <- list(Discipline = "allDisciplines", Journal = "allJournals", Country = country, Position = position, Points = data.for.web.app[[1]], Curve = data.for.web.app[[2]])
      file.name <- paste("allJournals", "allDisciplines", position, country, sep = "_")
    }
    
    if(filter.type == "everything.for.all.countries"){
      data.for.web.app <- list(Discipline = "allDisciplines", Journal = "allJournals", Country = "allCountries", Position = position, Points = data.for.web.app[[1]], Curve = data.for.web.app[[2]])
      file.name <- paste("allJournals", "allDisciplines", position, "allCountries", sep = "_")
    }
    
    # if(return.file.name.only) return(file.name)
    
    file.name <- paste("../outputs/json files for web app/", file.name, ".json", sep = "")
    write(toJSON(data.for.web.app, pretty = T), file.name)
  }
  ############ END OF INTERNAL FUNCTIONS
  
  
  ###### GET THE FOCAL SET OF DATA
  if(!(filter.type %in% c("arxiv.cat", "arxiv.subcat"))){ # If we're doing PubMed data...
    # Get the data out of the data.frame or database connection using dplyr. Also restrict the data to 2002 onwards, since the data are very sporadic before this 
    if(filter.type == "disc") real.data <- data.source %>% filter(discipline == item) %>% select(date, country, gender.95) %>% rename(gender = gender.95)
    if(filter.type == "journal") real.data <- data.source %>% filter(journal == item) %>% select(date, country, gender.95) %>% rename(gender = gender.95)
    if(filter.type %in% c("disc", "journal") & country == "all") real.data <- real.data %>% select(-country)
    if(filter.type == "everything.for.all.countries") real.data <- data.source # for when we want to know the change in gender ratio for all of PubMed, irrespective of journal/disc/country!
    if(filter.type != "everything.for.one.country") real.data <- real.data %>% collect(n = Inf) %>% as.data.frame %>% mutate(date = convert.dates(date, type = "since2000")) %>% filter(date >= 2)
    
    if(filter.type == "everything.for.one.country") {
      if(country == "all") return("check the arguments") # safeguard against accidental massive operation
      real.data <- data.source %>% filter(country != "NA" & !(country %in% unique.countries[unique.countries != country])) # discard all papers with no country, and those where country == wrong country
      real.data <- real.data %>% select(date, country, gender.95) %>% rename(gender = gender.95) %>% collect(n = Inf) %>% as.data.frame %>% mutate(date = convert.dates(date, type = "since2000")) %>% filter(date >= 2)
      real.data <- real.data[grep(country, real.data$country),] # throw out papers where the focal country is not listed
    }
  }
  
  if(filter.type == "arxiv.cat"){ # for ArXiv data, there is more data from the 1990s than on PubMed. Restrict to 1998 onwards instead of 2002 (1998 has 8119 preprints - drops to 5492 in 1997)
    real.data <- data.source %>% filter(category == item) %>% select(date, gender.95) %>% mutate(date = convert.dates(date, type = "since2000")) %>% filter(as.numeric(substr(date, 7, 10)) >= 1998) %>% rename(gender = gender.95) 
  }
  
  if(filter.type == "arxiv.subcat"){
    real.data <- data.source %>% filter(paste(category, sub.category, sep = ": ") == item) %>% select(date, gender.95) %>% mutate(date = convert.dates(date, type = "since2000")) %>% filter(as.numeric(substr(date, 7, 10)) >= 1998) %>% rename(gender = gender.95) 
    category <- strsplit(item, split = ": ")[[1]][1]
    sub.category <- strsplit(item, split = ": ")[[1]][2]
  }
  
  # Restrict the data to the specified authorship position
  if(position == "single") real.data <- real.data[nchar(real.data$gender) == 1, ] # if doing single, first or last, throw out some of the data
  else if(position != "overall"){
    real.data <- real.data[nchar(real.data$gender) > 1, ] # restrict to multi-author papers for "first" and "last" author analyses
    if(position == "first") real.data$gender <- substr(real.data$gender, 1, 1)
    if(position == "last") {
      length.author.list <- nchar(real.data$gender)
      real.data$gender <- substr(real.data$gender, length.author.list, length.author.list)
    }
  }
  
  # If specified, trim the data down to a single country as well
  if(country != "all"){
    real.data <- real.data[str_count("U", real.data$gender) != nchar(real.data$gender), ] # Throw out papers where all authors' genders are unknown
    real.data <- real.data[!is.na(real.data$country), ] # Throw out papers where the country is recorded as NA
    real.data <- real.data[str_detect(real.data$country, country), ] # Throw out papers where the focal country was not listed anywhere
    if(nrow(real.data) == 0) {
      # print("There are no data for this country, giving up")
      return(NULL) #  If there are no data, just end the whole function and return nothing
    }
    split.country <- strsplit(real.data$country, split = "_") # split the string holding the countries
    n.countries <- sapply(split.country, length)    
    to.keep <- n.countries == 1 | n.countries == nchar(real.data$gender) # To be safe, throw out papers with multiple countries where this number is not the same as number of authors
    real.data <- real.data[to.keep, ]
    split.country <- split.country[to.keep]; rm(to.keep)
    n.countries <- sapply(split.country, length)
    if(position != "overall") { # for first/last/single
      if(position == "last") real.data$country <- sapply(split.country, function(x) x[length(x)]) # get last listed country
      else real.data$country <- sapply(split.country, function(x) x[1]) # get first/single country
      real.data <- real.data[str_detect(real.data$country, country), ] # throw out papers where the focal country did not belong to the author we wanted (e.g. first or last)
    }
    else{ # for "overall"...
      multi.country.real.data <- split.country[n.countries > 1]
      if(length(multi.country.real.data) > 0){
        # Trim off all but the authors that are from the correct country
        # We will make the assumption that all the authors from a multi-author paper for which one country is listed are from that country
        paperID <- unlist(mapply(rep, 1:length(multi.country.real.data), each = sapply(multi.country.real.data, length))) 
        multi.country.real.data <- unlist(multi.country.real.data)
        keep.author.or.not <- multi.country.real.data == country
        multi.country.split.gender <- unlist(strsplit(real.data$gender[n.countries > 1], split = ""))     
        multi.country.split.gender[unlist(keep.author.or.not)] <- "U" # Set author genders from the wrong country to "U"
        real.data$gender[n.countries > 1] <- as.character(tapply(multi.country.split.gender, paperID, function(x) paste0(x, collapse = "")))
      }
    }
    real.data <- real.data %>% select(-country) # we are done with the country column now
  }
  
  # Throw out papers with all unknown authors, and count the number of papers remaining
  real.data <- real.data[str_count("U", real.data$gender) != nchar(real.data$gender) & !is.na(real.data$gender), ] 
  n.papers <- nrow(real.data)
  
  
  ###################  If there are no data, just end the function and return nothing
  if(n.papers == 0) {
    #  print("There are no data for this item, giving up")
    return(NULL)
  }
  
  ################## Count the number of unique papers in each year - needed for the web app mouse-over information
  n.papers.per.year <- real.data %>% mutate(year = floor(date)) %>% group_by(year) %>% summarise(n = n()) %>% as.data.frame()
  n.papers.per.year$year <- 2000 + n.papers.per.year$year
  
  ###################   Count the number of male and female authors observed on each unique date - we use this data to generate resamples efficiently
  real.data <- real.data %>% group_by(date) %>% summarise(F = sum(str_count(gender, "F")), M =  sum(str_count(gender, "M"))) %>% as.data.frame %>% gather(gender, count, -date) 
  probabilities <- with(real.data, count / sum(count))
  n.authors <- sum(real.data$count) # The number of M and F authors in the dataset - each resample will consist of this many authors, drawn with replacement 
  
  ################### If there are too few data, just end the function and return nothing
  if(n.papers < 100 | n.authors < 250 | max(real.data$date) - min(real.data$date) < 5) { # <- note the criteria for attempting to fit a curve
    #  print("There are too few papers, authors, or years to fit a curve, giving up")
    return(NULL)   
  }
  
  ###################  If doing a journal, we need to also get the discipline for the focal journal from the journals database
  if(filter.type == "journal") discipline <- (journals_sqlite %>% filter(short.title == item) %>% select(Discipline) %>% collect() %>% as.data.frame())[,1][1] # note the [1] ensures we get only 1 disc in cases where there are two journals with the same "short.title" listed in journals_sqlite. For example, PubMed lists "Journal of acquired immune deficiency syndromes" and "Journal of acquired immune deficiency syndromes (1999)" as having the same abbreviation, namely "J Acquir Immune Defic Syndr"
  
  ################### NOW RESAMPLE THE DATA TO DETERMINE THE 95% CIs on parameters like current gender ratio, rate of change, and years to parity (takes a while, hence option to skip it)
  if(run.sim){ # If run.sim = T, then we will be doing lots of time-consuming resampling and exporting the data to make Figures 1 and 2 in the paper
    
    # Function to make "chunk.size" resampled datasets containing n.authors each - this will be done 'n.chunks' times (too big for memory if done in one go)
    make.resampled.data <- function(real.data, n.authors, chunk.size){
      resampled.data <- sample(1:nrow(real.data), size = chunk.size * n.authors, replace = T, prob = probabilities)
      resampled.data <- melt(table(factor(resampled.data, levels = 1:nrow(real.data)), rep(1:chunk.size, each = n.authors))) %>% rename(row = Var1, boot.replicate = Var2, count = value) %>% mutate(row = as.numeric(as.character(row)))
      resampled.data <- data.frame(date = real.data$date[resampled.data$row], gender = real.data$gender[resampled.data$row], boot.replicate = resampled.data$boot.replicate, count = resampled.data$count)
      resampled.data <- resampled.data %>% group_by(boot.replicate, date) %>% summarise(nFemales = sum(count[gender == "F"]), nMales = sum(count[gender == "M"]), n = nFemales+nMales) %>% as.data.frame %>% filter(n > 0)
      resampled.data
    }
    
    # This function runs the optimiser on a focal bit of resampled data, and gets the current gender ratio, rate of change, 
    # and the year at which the ratio got within 5% of gender parity using the pfunc model with optimised parameters
    find.response.variables <- function(data){
      optim.results <- run.optimiser(data)
      r <- optim.results$par[1]
      c <- optim.results$par[2]
      gender.ratio.at.present <- pfunc(16.63518, r, c)    # Compute the predicted gender ratio at the "present", i.e. the day I collected the data, 20/8/16
      current.rate.of.change <- pfunc.deriv(p = gender.ratio.at.present, r = r) # Compute the rate of change in the gender ratio at the "present"
      year.range <- 0:9999
      predicted.gender.ratio <- pfunc(year.range, r, c)   # Compute the predicted gender ratio from year 2000 way into the future
      ifelse(predicted.gender.ratio[2] > predicted.gender.ratio[1], rising <- TRUE, rising <- FALSE)  # Check if it's rising or falling towards 50:50
      if(rising) {
        if(gender.ratio.at.present < 0.5) parity.year <- year.range[which(predicted.gender.ratio > 0.45)][1]                 # record first year it got within 5% of parity
        else(parity.year <- "Female-biased and becoming more so")
      }
      else {
        if(gender.ratio.at.present > 0.5) parity.year <- year.range[which(predicted.gender.ratio < 0.55)][1]
        else(parity.year <- "Male-biased and becoming more so")
      }
      return(data.frame(gender.ratio.at.present=gender.ratio.at.present, current.rate.of.change=current.rate.of.change, parity.year=parity.year, stringsAsFactors = F))
    }
    
    # Now run "find.response.variables" on each set of resampled data, and get the estimated 95% CIs using quantile
    resampled.data <- make.resampled.data(real.data, n.authors, chunk.size)
    boot.estimates <- do.call("rbind", lapply(1:chunk.size, function(i) find.response.variables(resampled.data %>% filter(boot.replicate == i))))
    rm(resampled.data) # help keep memory use low
    if(nChunks > 1){
      for(ii in 2:nChunks){
        resampled.data <- make.resampled.data(real.data, n.authors, chunk.size)
        boot.estimates <- rbind(boot.estimates, do.call("rbind", lapply(1:chunk.size, function(i) find.response.variables(resampled.data %>% filter(boot.replicate == i)))))
        rm(resampled.data) # help keep memory use low
      }
    }
    
    CIs.1 <- as.numeric(quantile(boot.estimates$gender.ratio.at.present, probs = c(0.025, 0.975)))
    CIs.2 <- as.numeric(quantile(boot.estimates$current.rate.of.change, probs = c(0.025, 0.975)))
    if(is.numeric(boot.estimates$parity.year)) {
      parity.year <- median(boot.estimates$parity.year + 2000) - 2016.63518
      CIs.3 <- as.numeric(quantile(boot.estimates$parity.year + 2000, probs = c(0.025, 0.975))) - 2016.63518
      parity.year[parity.year < 0] <- 0
      CIs.3[CIs.3 < 0] <- 0
    }
    else {
      Mode <- function(x) {ux <- unique(x); ux[which.max(tabulate(match(x, ux)))]} # if not going to hit parity sometimes, pick the modal (i.e. most likely) non-parity outcome
      parity.year <- Mode(boot.estimates$parity.year)
      CIs.3 <- c(NA, NA)
    }
    
  } ######## end of if(run.sim)
  
  ############################# GET R AND C VALUES FROM THE REAL DATA, AS WELL AS THE DATA POINTS NEEDED FOR A GRAPH (here, or in the web app) 
  # This is the main thing needed if export.JSON = T
  real.data <- real.data %>% group_by(date) %>% summarise(nFemales = sum(count[gender == "F"]), nMales = sum(count[gender == "M"]), n = nFemales+nMales, proportion.F = nFemales/n) %>% as.data.frame %>% filter(n > 0)
  result <- run.optimiser(real.data)
  real.r <- result$par[1]
  real.c <- result$par[2]
  
  
  # Find the gender ratio per year from 2002 to 2016 - note that only years with more than 50 authors measured will appear
  overlay.points <- real.data %>% mutate(date = floor(date)) %>% group_by(date) %>% 
    summarise(percent.F = 100 * sum(nFemales) / (sum(nFemales) + sum(nMales)), # % females in focal year
              nFemales = sum(nFemales), nMales = sum(nMales)) %>% # number of male and female authors in the focal year
    filter(nFemales + nMales > 50) %>% # Only keep years with at least 50 authors measured
    mutate(Year = date + 2000) %>% # add a nicer year column, and reorder the columns neatly
    select(Year, percent.F, nFemales, nMales)
  
  overlay.points$n <- n.papers.per.year$n[match(overlay.points$Year, n.papers.per.year$year)] # Add the number of papers sampled in each year (saved earlier)
  overlay.points$n[is.na(overlay.points$n)] <- 0
  
  if(nrow(overlay.points) < 5) {
    # print("There are too few years, giving up")
    return(NULL) # return nothing if we do not have at least 50 gendered authors from at least 5 different years (in addition to having at least 250 authors and 100 papers overall - see above)
  }
  
  
  overlay.points <- cbind(overlay.points, get.CIs(overlay.points$nFemales, overlay.points$nMales))
  
  if(plot){ # Make a nice plot if requested 
    model.fit <- data.frame(time = seq(0, 50, by = 0.1)) %>% mutate(percent.F = 100 * pfunc(time, real.r, real.c), Year = time + 2000)
    print( ggplot() + geom_errorbar(data = overlay.points, aes(x = Year, ymin = lowerCI, ymax = upperCI), colour = "darkgrey", width = 0) + 
             geom_point(data = overlay.points, aes(x = Year, y = percent.F)) + geom_line(data = model.fit, aes(x = Year, y = percent.F)) + ylab("% women authors")  )
  }
  
  if(export.JSON){    # Note: function(t, r, c){100 * exp(0.5*r*(t-2000))/(2*exp(0.5*r*(t-2000))+c)} can be used to predict the % females, where t is the date such that "2004.25" is quarter of the way through 2004. Use this for the web app Javascript
    data.for.web.app <- list(
      overlay.points %>% select(Year, nFemales, nMales, n, percent.F, lowerCI, upperCI) %>%
        rename(Y = Year, F = nFemales, M = nMales, GR = percent.F, lc = lowerCI, uc = upperCI) %>%  
        mutate(GR = round(GR,1), lc = round(lc,1),  uc = round(uc,1)),    # note rounding to save hard disk space and download size, important for the web
      data.frame(r = signif(real.r,2), c = signif(real.c,2)))
    run.JSON.export(data.for.web.app)
  }
  
  if(run.sim){
    notes <- ""
    if(is.numeric(parity.year) & is.na(CIs.3[1])) notes <- "Varying or unmeasurable trends in gender ratio"
    if(is.numeric(parity.year) & is.na(CIs.3[2])) notes <- "Varying or unmeasurable trends in gender ratio"
    
    out <- data.frame(journal = item, discipline = item, position = position, n.authors = n.authors, n.papers = n.papers,
                      gender.ratio.at.present = median(boot.estimates[,1]), low.CI.1 = CIs.1[1], up.CI.1 = CIs.1[2],
                      current.rate.of.change = median(boot.estimates[,2]), low.CI.2 = CIs.2[1], up.CI.2 = CIs.2[2],
                      years.to.parity = parity.year, low.CI.3 = CIs.3[1], up.CI.3 = CIs.3[2], r=real.r, c=real.c, notes = notes, stringsAsFactors = F)
    out[6:11] <- out[6:11] * 100 # express as percentages
    if(filter.type == "disc") out$journal <- "All journals"
    if(filter.type == "journal") out$discipline <- discipline  
    if(filter.type == "arxiv.cat") {
      out <- out %>% rename(category = discipline, sub.category = journal)
      out$sub.category <- "All subcategories" 
    }
    if(filter.type == "arxiv.subcat") {
      out <- out %>% rename(category = discipline, sub.category = journal)
      out$category <- category   # The cat and sub-cat were saved way up above
      out$sub.category <- sub.category
    }
    
    if(print.each.one) print(out)
    return(out)
  }
  
} # End of big gender.stats function

# wrapper for the gender.stats function used to make the web app data - feeds in a dataframe of arguments, and does one row from it. For use with lapply
make.web.app.data <- function(row, arguments){
  if(row %% 10 == 0) print(row)
  output <- "not done yet"
  attempts.left <- 5   # Very occasionally fails when the curve is almost totally flat, due to negative numbers in the likelihood or something. Just re-run when this happens
  while(attempts.left > 0) {
    try(output <- with(arguments, gender.stats(item = item[row], filter.type = filter.type[row], position = position[row], country = country[row], 
                                               data.source = pubmed_sqlite, run.sim = F, print.each.one = F, export.JSON = T, plot = F, verbose = F)))
    if(!is.character(output)) attempts.left <- 0 # If the function worked, "output" will have been overwritten with NULL, so we can quit already
    attempts.left <- attempts.left - 1
  }
}

# wrapper for the gender.stats function used to make the web app data - feeds in vector of disciplines or journals, or a custom set (dataframe with columns "item" and "position")
make.figure.one.dataset <- function(items = NULL, filter.type, author.positions = c("overall", "first", "last", "single"), custom.set = NULL, data.source = pubmed_sqlite, chunk.size = 10, nChunks = 100){
  if(is.null(custom.set) & is.null(items)) return("Enter something to do")
  if(!is.null(custom.set) & !is.null(items)) return("Pick one - 'custom set' or 'items'")
  if(is.null(custom.set) & !is.null(items)) to.do <- expand.grid(item = items, position = author.positions) %>% mutate(item = as.character(item), position = as.character(position))
  else to.do <- custom.set # you can enter a dataframe with columns "item" and "position"
  output <- do.call("rbind", lapply(1:nrow(to.do), function(xx) {
    successful <- FALSE; attempts <- 0; answer <- 0
    while(!successful & attempts < 10){ # Because they can fail stochastically, try each one a few times.
      try(answer <- gender.stats(item = to.do$item[xx], filter.type=filter.type, position = to.do$position[xx], data.source = data.source, chunk.size=chunk.size, nChunks=nChunks, 
                                 plot = F, print.each.one = T, export.JSON = F, run.sim = T))
      if(!is.numeric(answer)) successful <- TRUE
      attempts <- attempts + 1
    }
    if(attempts == 10) {
      print("Tried 10 times, giving up")
      return(NULL)
    }
    else return(answer)
  }
  ))
  output <- output %>% mutate(position = factor(position, levels = c("first", "last", "single", "overall")))
  if("discipline" %in% names(output)) output <- output %>% arrange(discipline, position)
  if("category" %in% names(output)) output <- output %>% arrange(category, sub.category, position)
  output
}


# Get all possible combinations of journal/disc/country/author position, which we will need to make the data for the web app
make.list.of.combos.for.web.app <- function(unique.discs, unique.journals, unique.countries){
  arguments <- data.frame(item = c(unique.discs, unique.journals), 
                          filter.type = as.character(unlist(mapply(rep, c("disc", "journal"), each = c(length(unique.discs), length(unique.journals))))), 
                          stringsAsFactors = F)
  nr <- nrow(arguments)
  arguments <- arguments[rep(1:nr, 4), ]
  arguments$position <- rep(c("first", "last", "single", "overall"), each = nr)
  nr <- nrow(arguments)
  arguments <- arguments[rep(1:nr, 1 + length(unique.countries)), ]
  arguments$country <- rep(c("all", unique.countries), each = nr)
  suppl.args <- data.frame(item = "all", filter.type = "everything.for.one.country", position = NA, country = unique.countries)
  nr <- nrow(suppl.args)
  suppl.args <- suppl.args[rep(1:nr, 4), ]
  suppl.args$position <- rep(c("first", "last", "single", "overall"), each = nr)
  arguments <- rbind(arguments, suppl.args, data.frame(item = "all", filter.type = "everything.for.all.countries", position = c("first", "last", "single", "overall"), country = "all"))
  rm(suppl.args)
  to.keep <- rep(TRUE, nrow(arguments))
  to.keep[arguments$filter.type == "journal" & arguments$country != "all"] <- FALSE
  for(i in 1:length(unique.journals)){
    print(i)
    country.column.for.this.journal <- (pubmed_sqlite %>% filter(journal == unique.journals[i]) %>% select(country) %>% collect(n = Inf) %>% as.data.frame())[,1]
    country.column.for.this.journal <- paste0(unique(country.column.for.this.journal), collapse = "") # paste all unique affiliation lists it into a single item
    countries.present <- unique.countries[str_detect(country.column.for.this.journal, unique.countries)] # the names of all the countries represented in this journal
    to.keep[arguments$filter.type == "journal" & arguments$item == unique.journals[i] & arguments$country %in% countries.present] <- TRUE
  }
  
  arguments <- arguments[to.keep, ]
  arguments <- arguments[sample(nrow(arguments)), ] # shuffle the rows so that work is distributed evenly across cores
  row.names(arguments) <- NULL
  arguments
}

# Function used to tidy the data prior to export as supplementary files
neaten.the.data <- function(dd){
  dd[, names(dd) %in% c("gender.ratio.at.present", "low.CI.1", "up.CI.1", "current.rate.of.change", "low.CI.2", "up.CI.2", "low.CI.3", "up.CI.3", "r", "c")] <- round(dd[, names(dd) %in% c("gender.ratio.at.present", "low.CI.1", "up.CI.1", "current.rate.of.change", "low.CI.2", "up.CI.2", "low.CI.3", "up.CI.3")], 3)  # some rounding
  dd <- dd %>% select(-notes) # empty column removed
  dd$qualitative.result <- ""
  dd$qualitative.result[dd$gender.ratio.at.present < 50 & dd$current.rate.of.change > 0] <- "Male-biased and unchanging"   # Add these easier-to-read results
  dd$qualitative.result[dd$gender.ratio.at.present > 50 & dd$current.rate.of.change > 0] <- "Female-biased and unchanging"
  dd$qualitative.result[dd$gender.ratio.at.present < 50 & dd$current.rate.of.change > 0] <- "Male-biased but improving"   
  dd$qualitative.result[dd$gender.ratio.at.present > 50 & dd$current.rate.of.change < 0] <- "Female-biased but improving"
  dd$qualitative.result[dd$gender.ratio.at.present < 50 & dd$current.rate.of.change < 0] <- "Male-biased and becoming more so"
  dd$qualitative.result[dd$gender.ratio.at.present > 50 & dd$current.rate.of.change > 0] <- "Female-biased and becoming more so"
  dd$qualitative.result[dd$gender.ratio.at.present > 45 & dd$gender.ratio.at.present < 55] <- "Essentially at parity"
  dd$years.to.parity <- round(suppressWarnings(as.numeric(dd$years.to.parity)),3)  # Tidy this, which originally had characters mixed in
  dd
}


clean.up.data.for.figures <- function(df, filter.type){
  df$qualitative.result <- NA
  df$qualitative.result[!is.na(suppressWarnings(as.numeric(df$years.to.parity)))] <- NA
  df$qualitative.result[df$gender.ratio.at.present < 50 & df$current.rate.of.change < 0] <- "Male-biased and becoming more so"
  df$qualitative.result[df$gender.ratio.at.present > 50 & df$current.rate.of.change > 0] <- "Female-biased and becoming more so"
  df$qualitative.result[df$gender.ratio.at.present > 45 & df$gender.ratio.at.present < 55] <- "Essentially at parity"
  df$low.CI.3[df$gender.ratio.at.present > 45 & df$gender.ratio.at.present < 55] <- NA
  df$up.CI.3[df$gender.ratio.at.present > 45 & df$gender.ratio.at.present < 55] <- NA
  
  # Remove the dots for items that might be either going up or going down, and we're not sure (i.e. the rate of change is indistinguishable from zero)
  df$years.to.parity[is.na(df$low.CI.3)] <- NA
  df$years.to.parity[is.na(df$up.CI.3)] <- NA
  
  df$years.to.parity <- suppressWarnings(as.numeric(df$years.to.parity))
  df$years.to.parity[df$qualitative.result == "Essentially at parity"] <- 0
  df$position[df$position == "overall"] <- "Overall"
  df$position[df$position == "first"] <- "First"
  df$position[df$position == "last"] <- "Last"
  df$position[df$position == "single"] <- "Single"
  df$position <- factor(df$position, levels = c("Overall", "First", "Last", "Single"))
  if(filter.type == "disc") df <- df %>% arrange(position, gender.ratio.at.present) %>% mutate(discipline = factor(discipline, levels = unique(discipline)))
  if(filter.type == "journal") df <- df %>% arrange(position, discipline, gender.ratio.at.present) %>% mutate(journal = factor(journal, levels = unique(journal)))
  df$position <- factor(df$position, levels = c("First", "Last", "Single", "Overall"))
  return(df)
}


# Big function to make figures lik Figs 1 and 2
main.figure <- function(dataset, xlim1, xlim2 = c(-2,2), add.density = T, add.density3 = T, custom.remove.first = F, single.present = T){
  cols <- c(brewer.pal(3, "Set1"), "black")
  if(!single.present) cols <- c(brewer.pal(3, "Set1")[1:2], "black")
  longest.name <- dataset$discipline[nchar(as.character(dataset$discipline)) == max(nchar(as.character(dataset$discipline)))][1]
  
  d1 <- dataset %>% ggplot(aes(x = gender.ratio.at.present, fill = position, colour = position)) + theme_classic() + 
    theme(legend.position = "none", title = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank(), axis.text.y = element_text(colour = "white"), axis.line.y = element_blank()) + 
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0), breaks = 0, labels = longest.name) + coord_cartesian(xlim = xlim1) + 
    geom_density(alpha=0.1) + xlab(NULL) + ylab("Frequency") + scale_fill_manual(values = cols) + scale_colour_manual(values = cols) 
  
  d2 <- dataset %>% ggplot(aes(x = current.rate.of.change,  fill = position, colour = position)) + theme_classic() + 
    theme(legend.position = "none", text = element_blank(), title = element_blank(), axis.ticks = element_blank(), axis.line.y = element_blank()) +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + coord_cartesian(xlim = xlim2) +  
    geom_density(alpha=0.1) + xlab(NULL) + ylab(" ") + scale_fill_manual(values = cols) + scale_colour_manual(values = cols) + theme(legend.position = "none")
  
  d3.data <- dataset
  if(custom.remove.first){
    d3.data <- dataset %>% filter(position != "First") %>% mutate(position = factor(position, levels = c("Last", "Single", "Overall")))
    new.cols <- c(brewer.pal(3, "Set1")[2:3], "black")
  }
  
  d3 <- d3.data %>% ggplot(aes(x = years.to.parity, fill = position, colour = position)) + theme_classic() + 
    theme(legend.position = "none", text = element_blank(), title = element_blank(), axis.ticks = element_blank(), axis.line.y = element_blank()) +
    scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0)) + coord_cartesian(xlim = c(0,100)) + xlab(NULL) + ylab(" ") + theme(legend.position = "none")
  
  
  if(add.density3) {
    if(!custom.remove.first) d3 <- d3 + geom_density(alpha=0.1)  + scale_fill_manual(values = cols) + scale_colour_manual(values = cols)
    if(custom.remove.first) d3 <- d3 + geom_density(alpha=0.1)  + scale_fill_manual(values = new.cols) + scale_colour_manual(values = new.cols)
  }
  else d3 <- d3 + theme(axis.line.x = element_blank())
  
  p1 <- dataset %>% ggplot(aes(x = gender.ratio.at.present, y = discipline, colour = position)) + geom_vline(xintercept = 50, linetype=2) +
    geom_errorbarh(alpha = 0.6, height = 0, aes(xmin = low.CI.1, xmax = up.CI.1)) + geom_point(alpha=0.8) + xlab("% women authors\nin 2016") + ylab(NULL) + scale_colour_manual(values = cols) + labs(colour = "Author\nposition") + coord_cartesian(xlim = xlim1)
  
  p2 <- dataset %>% ggplot(aes(x = current.rate.of.change, y = discipline, colour = position)) + geom_vline(xintercept = 0, linetype=2) +
    geom_errorbarh(alpha = 0.6, height = 0, aes(xmin = low.CI.2, xmax = up.CI.2)) + geom_point(alpha=0.8) + coord_cartesian(xlim = xlim2) +
    xlab("Change in % women\nauthors per year") + ylab(NULL)  + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none") + scale_colour_manual(values = cols)
  
  p3 <- dataset  %>% ggplot(aes(x = years.to.parity, y = discipline, colour = position))  +
    geom_errorbarh(alpha = 0.6, height = 0, aes(xmin = low.CI.3, xmax = up.CI.3)) + geom_point(alpha=0.8) +
    xlab("Years until parity") + ylab(NULL) + theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())  + scale_colour_manual(values = cols) + coord_cartesian(xlim = c(0, 100))
  
  g <- ggplotGrob(p1 + theme(legend.position = "right"))$grobs # Get the legend from p1 
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lwidth <- sum(legend$width) # save the width of the legend
  
  
  gl.d <- lapply(list(d1,d2,d3), function(x) ggplotGrob(x + theme(legend.position="none"))) # convert the plots to ggplotGrobs to allow cbind
  gl <- lapply(list(p1,p2,p3), function(x) ggplotGrob(x + theme(legend.position="none"))) # convert the plots to ggplotGrobs to allow cbind
  
  if(add.density == T){
    plots.combined <- arrangeGrob(cbind(gl.d[1][[1]], gl.d[2][[1]], gl.d[3][[1]]), 
                                  cbind(gl[1][[1]], gl[2][[1]], gl[3][[1]]),
                                  nrow = 2,
                                  heights = c(0.1,0.9))
  }
  else plots.combined <- arrangeGrob(cbind(gl[1][[1]], gl[2][[1]], gl[3][[1]]))
  
  combined <- arrangeGrob(plots.combined, legend, ncol = 2, widths = unit.c(unit(1, "npc") - lwidth, lwidth)) # add legend
  grid.newpage()
  grid.draw(combined) # draw it
}