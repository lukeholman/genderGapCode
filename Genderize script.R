library(stringr)
library(dplyr)
library(jsonlite)
api.key <- "4da1b848571732e83d106e2e5c44ccdb"

# The Genderize API often gives a time-out error, so when we enter a search and it fails, try, try again!
# In my experience, it never fails more than once, and the second try works fine. But I added 50 tries to be safe
call.api <- function(search){   
  for(i in 1:50){
    try(return(fromJSON(search)), silent = T)
    print("The API failed, trying again")
    Sys.sleep(2)
  }
  print("The API failed 50 times, giving up")
}

# Function that dials up the genderise API, and asks for information on some names (used inside 'do.many.names')
# Optionally, one can specify the country these names are from - this accounts for e.g. Kim being a male name in Denmark, female elsewhere
# names.already.done is an optional vector of names for which we have already searched the API for gender info
do.chunk.of.names <- function(names, country.code = NULL, names.already.done = NULL){
  base.url <- paste("https://api.genderize.io/?apikey=", api.key, sep="") # The address and API key
  
  if(is.null(country.code)) search.result <- call.api(paste0(c(base.url, names), collapse = "&name[]="))
  else search.result <- call.api(paste(paste0(c(base.url, names), collapse = "&name[]="), "&country_id=", country.code, sep = ""))
  
  # standardise the search output into a dataframe with these columns
  output <- data.frame(name = search.result$name, gender=NA, prob = NA, count = NA, country= NA, stringsAsFactors = F)
  if(!is.null(country.code)) output$country <- country.code
  if("gender" %in% names(search.result)) {
    if(!is.null(search.result$gender)) output$gender <- search.result$gender} # NOT THIS
  if("probability" %in% names(search.result)) { # These columns are not always in the output, so only add them if present
    output$prob <- search.result$probability
    output$count <- search.result$count
  }
  
  if(!is.null(country.code)){ # If a country was specified, do another search to try to get the missing names with an unlocalised search
    if(sum(is.na(output$gender)) > 0){
      
      names <- output$name[is.na(output$gender)] # These are the names that were not successfully gendered yet
      
      # First, let's check if the names are in 'names.already.done'. If so, don't search for them again
      names <- names[!(names %in% names.already.done)]
      
      if(length(names) > 0){ # If there are still some names to search
        
        # Append the new names to 'names.already.done', and do the search
        names.already.done <- c(names.already.done, names)
        new.search <- call.api(paste0(c(base.url, names), collapse = "&name[]="))
        
        # Format the new search properly
        extra.output <- data.frame(name = new.search$name, gender=NA, prob = NA, count = NA, country= "unspec", stringsAsFactors = F)
        if("gender" %in% names(new.search)) { 
          if(!is.null(new.search$gender)) extra.output$gender <- new.search$gender}
        if("probability" %in% names(new.search)) { # These columns are not always in the output, so only add them if present
          extra.output$prob <- new.search$probability
          extra.output$count <- new.search$count}
        extra.output <- extra.output[!is.na(extra.output$gender), ]
        # Add in the new search hits, if any
        if(nrow(extra.output) > 0) output[match(extra.output$name, output$name), ] <- extra.output
      }
    }
  }
  
  if(is.null(names.already.done)) return(output)
  else(return(list(output, names.already.done))) # return the name hits and the up-to-date version of names.already.done
}

# Takes a big list of names, and calls up the API for answers in chunks of 10 names or fewer (10 is the max allowed)
do.many.names <- function(names, country.code = NULL, names.already.done = NULL){
  if(length(names) > 10){
    chunk <- rep(1:(floor(length(names) / 10)), each = 10) # vector of same length as 'names', saying what chunk the name belongs to
    chunk <- c(chunk, rep(max(chunk)+1, length(names) %% 10))
    # Pre-allocate memory for the answer
    output <- data.frame(name = names, gender=NA, prob = NA, count = NA, country= NA, stringsAsFactors = F)
    
    if(is.null(country.code)) { # If we didn't specify a country...
      for(i in 1:max(chunk)) { # Simply save the output
        output[chunk == i, ] <- do.chunk.of.names(names[chunk == i], country.code=country.code, names.already.done=names.already.done)
      }
    }
    else{ # if we DID specify a country...
      for(i in 1:max(chunk)) {
        focal.answer <- do.chunk.of.names(names[chunk == i], country.code=country.code, names.already.done=names.already.done)
        output[chunk == i, ] <- focal.answer[[1]] # store the answer
        names.already.done <- focal.answer[[2]] # update names.already.done
      }
    }
  }
  
  # Or if we entered 10 or fewer names...
  else output <- do.chunk.of.names(names=names, country.code=country.code)
  
  output$gender[is.na(output$gender)] <- "U"         # Relabel the genders to M, F or U
  output$gender[output$gender == "female"] <- "F"
  output$gender[output$gender == "male"] <- "M"
  
  if(is.null(names.already.done)) return(output)
  else(return(list(output, names.already.done))) # return the name hits and the up-to-date version of names.already.done
}


do.many.names.in.blocks <- function(names, output.dir, country.code = "unspec", names.already.done = NULL, chunk.size = 1000){
  
  if(length(names) > chunk.size){
    # Define start and end points of each chunk that we will write to disk
    chunk.starts <- seq(from = 1, to = 1 + chunk.size*floor(length(names)/chunk.size), by = chunk.size) 
    chunk.ends <- seq(from = chunk.size, to = chunk.size*ceiling(length(names)/chunk.size), by = chunk.size)
    chunk.ends[length(chunk.ends)] <- length(names)
  }
  else{
    chunk.starts <- 1
    chunk.ends <- length(names)
  }
  
  # Define filepaths for the output files
  output.files <- paste(output.dir, "/", country.code, 1:length(chunk.starts), ".csv", sep = "")
  done.files <- list.files(output.dir, full.names = T) # list all the files that are already completed
  remaining.to.do <- (1:length(chunk.starts))[!(output.files %in% done.files)] # numbers of files that are not yet done
  
  if(length(remaining.to.do) == 0){
    print(paste("Country", country.code, "has been completed."))
    return("Completed already") # returns something if NO NEW DATA were produced
  }
  
  for(i in 1:length(remaining.to.do)){
    focal.chunk <- remaining.to.do[i]
    print(paste("Doing country '", country.code, "', chunk ", focal.chunk, ".", sep = ""))
    focal.names <- names[chunk.starts[focal.chunk]:chunk.ends[focal.chunk]]
    if(country.code != "unspec") {
      answer <- do.many.names(focal.names, country.code = country.code, names.already.done = names.already.done)
      chunk.o.names <- answer[[1]]
      names.already.done <- answer[[2]]
    }
    else chunk.o.names <- do.many.names(focal.names)
    write.csv(chunk.o.names, file = output.files[focal.chunk], row.names = F)
  }
  return(NULL) # returns null if new data were produced
}

# Code used to generate a list of unique names, split up by country, from the big dataset
# Load up the database
pubmed_db <- src_sqlite("../data for analysis/Holman_pubmed_db.sqlite3", create = F)
pubmed_sqlite <- tbl(pubmed_db, "papers")
countries <- get.unique.countries(c("../data/Update files - parsed affiliations", "../data/Main files - parsed affiliations"))
acceptable.countries <- c("japan", tail(parse.countries(countries, type = "summary"), 117)[,1])
acceptable.countries <- acceptable.countries[!(acceptable.countries %in% c("kingdom", "united"))] # check these are not in there

unique.names.by.country <- lapply(acceptable.countries, function(focal.country){
  print(focal.country)
  focal <- (pubmed_sqlite %>% filter(country == focal.country) %>% select(forenames) %>% collect(n=Inf) %>% as.data.frame())[,1] %>% strsplit(split="_") %>% unlist() %>% unique()
  focal <- focal %>% tolower() %>% unique()
  focal <- gsub("[',./]", "", focal)
  focal <- focal[nchar(focal) > 2]
  hyphenated <- which(str_detect(focal, "-"))
  for(i in hyphenated){
    split <- strsplit(focal[i], split = "-")[[1]]
    focal[i] <- paste0(split[nchar(split) > 1], collapse = "-")
  }
  rm(hyphenated)
  unique(focal[nchar(focal) > 2])
})
names(unique.names.by.country) <- acceptable.countries
save(unique.names.by.country, file = "../outputs/unique.names.by.country.Rdata")


# Load the names to gender
load("../outputs/unique.names.by.country.Rdata")
# Split the names up
unique.names.by.country <- lapply(unique.names.by.country, function(x) unique(unlist(strsplit(x, split = "-"))))

# Here is the list of the country codes supported by the genderize API. Where possible, let's do a country-specific search first. If it gives no hits, move on to non-specific search.
code.list <- c("AF","AL","AM","AR","AZ","BA","BE","BG","BR","BY","CA","CL","CN","CO","CZ","DE","DK","EE","EO","ES","FI","FO","FR","GB","GE","GR","HK","HR","HU","ID","IE","IL","IN","IR","IS","IT","JP","KE","KH","KR","KS","KZ","LA","LK","LT","LV","MK","MM","MN","MT","MX","MY","NL","NO","NP","PE","PH","PI","PK","PL","PT","RO","RS","RU","SE","SI","SK","SO","TH","TR","TW","UA","UD","US","UZ","VA","VE","VN","ZA")
x <-read.csv("../data/country codes.csv", stringsAsFactors = F) %>% filter(code %in% code.list) %>% mutate(country = tolower(country))
x$country[x$country == "bosnia"] <- "bosnia and herzegovina"
x$country[x$country == "czechia"] <- "czech republic"
names(unique.names.by.country)[names(unique.names.by.country) %in% x$country] <- x$code[match(names(unique.names.by.country), x$country)][!is.na(x$code[match(names(unique.names.by.country), x$country)])]
countryless.names <- unique(unlist(unique.names.by.country[nchar(names(unique.names.by.country)) > 2]))
unique.names.by.country <- unique.names.by.country[nchar(names(unique.names.by.country)) == 2]


# First do the unique names from all the countries that are not indexed in Genderize. This just looks for the name across all the countries
do.many.names.in.blocks(countryless.names, output.dir = "../temp")
countryless.names <- do.call("rbind", lapply(list.files("../temp", full.names = T), read.csv, stringsAsFactors = F))
countryless.names$country <- "unspec"
write.csv(countryless.names, file = "../outputs/names/countryless.names.csv")

# Now loop over every country on the API
countries.to.do <- names(unique.names.by.country)
names.already.done <- do.call("rbind", lapply(list.files("../outputs/names", full.names = T), read.csv, stringsAsFactors = F))
names.already.done <- names.already.done %>% filter(country == "unspec") %>% select(name)
# unique.names.by.country <- unique.names.by.country[!names(unique.names.by.country)%in%"TW"]
for(i in 1:length(countries.to.do)){
  # Do the focal country
  focal.country <- countries.to.do[i]
  tester <- do.many.names.in.blocks(unique.names.by.country[[which(names(unique.names.by.country) == focal.country)]], output.dir = "../temp", country.code = focal.country, names.already.done = names.already.done)
  
  if(is.null(tester)){ # If the do.many.names.in.blocks function made some new data...
    # Zip all the little temporary files together and save with a genders.csv ending in the 'outputs' directory
    zipped.names <- do.call("rbind", 
                            lapply(list.files("../temp", full.names = T)[grep(paste("/", focal.country, sep=""), list.files("../temp", full.names = T))],
                                   read.csv, stringsAsFactors = F))
    write.csv(zipped.names, file = paste("../outputs/names/", focal.country, "_genders.csv", sep = ""), row.names = F)
    
    # Now open up all the files that have been done so far and get all the unique names that have been searched already to make an up-to-date version of names.already.done, before we loop back to the start
    names.already.done <- do.call("rbind", lapply(list.files("../outputs/names", full.names = T), read.csv, stringsAsFactors = F))
    names.already.done <- names.already.done %>% filter(country == "unspec") %>% select(name)
  }
}


# Make a big table of all the names and countries they are from
all.country.names <- lapply(lapply(list.files("../outputs/names/", full.names = T), read.csv, stringsAsFactors = F), function(x) x[!is.na(x$gender), ])
names(all.country.names) <- x$country[match(gsub("_genders.csv", "", list.files("../outputs/names/")), x$code)]
names(all.country.names)[is.na(names(all.country.names))] <- "countryless.names"
all.country.names <- lapply(all.country.names, function(x) x[,!(names(x) %in% "X")])
all.country.names <- do.call("rbind", all.country.names)
all.country.names <- all.country.names[!duplicated(paste(all.country.names$name, all.country.names$country)), ]
all.country.names <- all.country.names %>% arrange(country, name)
all.country.names$country <- x$country[match(all.country.names$country, x$code)]

#write.csv(all.country.names, file = "../data/genderize.names.mastersheet.csv", row.names = F)
