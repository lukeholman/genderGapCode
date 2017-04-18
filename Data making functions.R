########################################################################################################
# PARSING THE PUBMED XML FILES INTO USEFUL DATA
########################################################################################################

# Function that pull out all the potentially useful information from the XML file entry for a single Medline citation. We need stuff like the author's forenames, the DOI, the PMID, the publication date and the journal. We also check if there is an abstract (or rather, how long it is), and if the paper is described as a 'journal article' (0=no, 1=yes)
get.items <- function(x){
    output <- data.frame(pmid=NA, doi=NA, date=NA, journal=NA, ISSN=NA, forenames=NA, abstract.length=NA, is.journal.article=NA)
    
    # Throughout, I use this tentative style of coding so the function doesn't fail when one or more fields are missing from the XML file (which they often are).
    
    pmid.local <- try(x$PMID$text) 
    if(exists("pmid.local")) {
      if(length(pmid.local) == 1){
      output$pmid <- paste("f1", pmid.local, sep ="___")
    }}
    
    # Searching for the DOI is slightly hard, because annoyingly there is no obvious, consistent 'doi' field in PubMed! Here, I define a DOI as a string that contains "10." followed by at least one digit, then a slash, then some letters and/or numbers. I search two different pubmed fields, which can contain the DOI and potentially also some other crap.
    record <- unlist(x) # unlist the record, to make it easier to search for DOI
    field <- names(record) # Store the names of each Pubmed field
    to.search.for.doi <- record[c(grep("ArticleIdList", field), grep("ELocationID", field))] # Search only the "ArticleIdList" and "ELocationID" fields
    doi.local <- try(unique(to.search.for.doi[grep("10[.][[:digit:]]+/[[:alnum:]]", to.search.for.doi)])[1]) # if there are multiple DOIs, take first one listed
    if(exists("doi.local")) {
      if(length(doi.local) == 1){
      output$doi <- paste("f2", doi.local, sep ="___")
    }}
    
    DMY.local <- try(paste0(rev(unlist(x$DateCreated)), collapse = "_")) # date in day_month_year format. R adds extra underscores for some puzzling reason!
    if(exists("DMY.local")) {
      if(length(DMY.local) == 1){
        output$DMY <- paste("f3", DMY.local, sep ="___")
    }}
        
    journal.local <- try(gsub("[.,]", "", x$Article$Journal$ISOAbbreviation))
    if(exists("journal.local")) {
      if(length(journal.local) == 1){
        output$journal <-paste("f4", journal.local, sep ="___")
    }}
    
    ISSN.local <- try(x$Article$Journal$ISSN$text)
    if(exists("ISSN.local")){
      if(length(ISSN.local) == 1){
      output$ISSN <- paste("f5", ISSN.local, sep ="___")
    }}
    
    forenames.local <- try(unlist(x$Article$AuthorList)) # Get all the names
    forenames.local <- try(forenames.local[grep("[Ff]ore[Nn]ame", names(forenames.local))]) # Restrict to just the forenames
    if(class(forenames.local) == "character") forenames.local <- try(paste0(as.character(sapply(strsplit(forenames.local, split = " "), head, 1)), collapse = "_"))
    if(exists("forenames.local")) {
      if(length(forenames.local) == 1){
      output$forenames <- paste("f6", forenames.local, sep ="___")
    }}

    abstract.parts <- try(unlist(x$Article$Abstract))
    abstract.parts <- try(abstract.parts[-c(grep("attrs", names(abstract.parts)), grep("[Cc]opyright", names(abstract.parts)))])
    abstract.length.local <- try(sum(nchar(abstract.parts))) # Record number of characters in the abstract.
    if(exists("abstract.length.local")){ # If it exists...
      if(length(abstract.length.local) == 1){ # Check it has proper length
      output$abstract.length <- paste("f7", abstract.length.local, sep ="___")
    }}
    
    is.journal.article.local <- try(1 * (length(grep("[Jj]ournal [Aa]rticle", unlist(x$Article$PublicationTypeList))) > 0))
    if(exists("is.journal.article.local")) { # If it exists...
      if(length(is.journal.article.local) == 1){ # Check it has proper length
        output$is.journal.article <- paste("f8", is.journal.article.local, sep ="___")
    }}
    return(output)
}

# This function takes an XML file describing many citations and parses into a nice R dataframe with the info that we need.
# "filepath" argument - the filepath to the XML file to be parsed.
parse.xml.file <- function(filepath){
  
  # This function is called on each individual paper in the XML file. It gets out the data we want, in an ugly text string. I separate fields with two tildas for easy cutting later.
  branch.fun <- function(x) print(paste0(get.items(xmlToList(x)),collapse = "~~"))
  
  # Run the branch function on every paper in the XML file to produce a rough and ready synopsis of each paper. Use capture.output() + print() as this seems the easiest way to deal with the XML package (handlers?? WTF?)
  rough.data <- capture.output(xmlEventParse(file = filepath, 
                               handlers = NULL, 
                               branches = list(MedlineCitation = branch.fun)))
  
  # Internal function that cleans up the ugly data into a sensible dataframe, ready to be returned by the main function.
  clean.rough.data <- function(rough.data){
    cleaned <- substr(rough.data, 1, nchar(rough.data)-1)
    cleaned <- substr(cleaned, 6, nchar(cleaned))
    cleaned <- strsplit(cleaned, split = "~~")
    cleaned <- do.call("rbind", lapply(cleaned, function(x){
      output <- data.frame(pmid=NA, doi=NA, date=NA, journal=NA, ISSN=NA, forenames=NA, abstract.length=NA, is.journal.article=NA)
      
      # I chose to add these text tags (e.g. f1___ means field one)
      # so that there was no chance of mixing up different fields
      # The tidy function looks for the right tag, then trims it off
      tidy <- function(item, tag){ 
        focal <- item[grep(paste(tag, "___", sep=""),item)]
        if(length(focal) == 1) return(gsub(paste(tag, "___", sep=""), "", focal))
        else return(NA)
      }
      
      output$pmid <- tidy(x, "f1")
      output$doi <- tidy(x, "f2")
      output$date <- tidy(x, "f3")
      output$journal <- tidy(x, "f4")
      output$ISSN <- tidy(x, "f5")
      output$forenames <- tidy(x, "f6")
      output$abstract.length <- tidy(x, "f7")
      output$is.journal.article <- tidy(x, "f8")
      return(output)
    }))
    cleaned$date <- gsub("__", "_", cleaned$date) # For some reason extra underscores get added to the date field. Trim them off them now.
    cleaned$date <- substr(cleaned$date, 2, nchar(cleaned$date)-1)
    
    # Delete duplicate PMIDs. If there are multiples, keep the PMID nearest the bottom (allowing later records to replace earlier ones).
    cleaned <- cleaned[rev(!duplicated(rev(cleaned$pmid))), ]
    
    return(cleaned)
  }
  return(clean.rough.data(rough.data))
}

# Function that runs 'parse.xml.file' on many xml files in a smart way. It is designed such that you can stop and start it without losing work already done. Should take only a few days to parse the entirity of Medline/Pubmed on a normal desktop computer.
# Arguments:
# - input.dir: The full filepath of a directory filled with XML files to be parsed. Note that the function can deal with .xml and also .xml.gz files (R unzips them for you, though this does take time). Should be something like "~/Desktop/Pubmed xml files"
# - output.dir: The full filepath of a directory to hold the parsed output files. Should be something like "~Desktop/Parsed Pubmed data"
# focal.chunk: If left as the default, this function will work through all the xml files. If set to a number less than or equal to num.chunks, it will do the specified chunk. For example you can open up 6 R processes, and run the fuction with a number from 1 to 6, to speed things up 6-fold.
# num.chunks = 6: The number of chunks to split the job into. Change this to a higher number than 6 if you have loads of cores or whatever.
# over.write: If set to TRUE, will overwrite existing files. If FALSE, the function checks if a file has already been parsed and placed in the output directory, and does not re-do it if it's there already.
parse.many.xml.files <- function(input.dir, output.dir, focal.chunk=NULL, num.chunks = 8, over.write = F){
  xml.files <- list.files(input.dir, full.names = T) # list all the files in the input
  xml.files <- xml.files[grep("medline",xml.files)] # Make sure they are all Medline files
  num.files <- length(xml.files) # count them
  if(!is.null(focal.chunk)){ # Specify the focal chunk of files
    chunk.size <- floor(num.files / num.chunks)
    if(focal.chunk > num.chunks) return("Check your inputs, they're wrong!")
     
     # Make a sequence of files. Spread evenly over the available ones
     sequence <- seq(from = focal.chunk, to = (num.chunks*chunk.size)+focal.chunk-1, by = num.chunks)
     #if final chunk, add on the bonus files
    if(focal.chunk == num.chunks) sequence <- c(sequence, (1:num.files)[(1:num.files) > max(sequence)])

    xml.files <- xml.files[sequence]
    num.files <- length(xml.files) # re-count the number of files in the subset
  }

  if(!over.write) completed.files <- list.files(output.dir, full.names = T) # list of file paths that were processed already
  # Loop over the focal files and do them in order, following overwrite rules. Write parsed files to disk in the output directory
  for(i in 1:num.files){ 
    print(paste("Doing file ", i, " out of ", num.files, ".", sep=""))
    output.file.path <- paste(output.dir, "/", gsub("[.]xml[.]gz","",tail(strsplit(xml.files[i], split = "/")[[1]],1)), ".csv", sep="")
    if(!(output.file.path %in% completed.files) | over.write) write.csv(parse.xml.file(xml.files[i]), file = output.file.path, row.names = F)
  }
}




# Function for reaeding the file downloaded from ftp://ftp.nlm.nih.gov/online/journals/, which has info about every journal archived on Pubmed (we mostly want the "MeSH headings", which give info on the scientific discipline of the journal)
parse.journal.descriptions.xml.file <- function(filepath){
  
  # Internal function to get the necessary info about each journal from the XML file of journal info provided by Pubmed, and available at their FTP: ftp://ftp.nlm.nih.gov/online/journals/
  get.journal.info <- function(xml.item){
    output <- data.frame(NLMid=NA, title=NA, short.title=NA, ISSN=NA, mesh.headings=NA, broad.journal.headings=NA, language=NA, country=NA)
    
    try(NLMid.local <- xml.item$NlmUniqueID[1])  # Note that the [1]s ensure that I only grab the first item if for some reason multiples are present
    if(exists("NLMid.local")) {
      if(length(NLMid.local) == 1){
        output$NLMid <- paste("f1", as.character(NLMid.local), sep ="___")
      }}
    
    try(title.local <- xml.item$Title[1])
    if(exists("title.local")) {
      if(length(title.local) == 1){
        output$title <- paste("f2", as.character(title.local), sep ="___")
      }}
    
    try(short.title.local <- xml.item$MedlineTA[1])
    if(exists("short.title.local")) {
      if(length(short.title.local) == 1){
        output$short.title <- paste("f3", as.character(short.title.local), sep ="___")
      }}
    
    try(ISSN.local <- xml.item$ISSN$text[1])
    if(exists("ISSN.local")) {
      if(length(ISSN.local) == 1){
        output$ISSN <- paste("f4", as.character(ISSN.local), sep ="___")
      }}
    
    # Hard bit to get the "MeSH headings", which are PubMed descriptors of the journal. Sometimes these are descriptors of its major topic, sometimes they are other stuff that we don't want
    # E.G. a medical history journal might say its major topic is history, but there will also be a heading saying that it is about English medical history - we just want the history part.
    try(mesh.headings.list <- xml.item$MeshHeadingList)
    if(exists("mesh.headings.list")) { # if the list exists and has something in it...
      if(length(mesh.headings.list) > 0){
        
        mesh.headings <- rep(NA, length(mesh.headings.list)) # declare empty holder for the items
        
        for(i in 1:length(mesh.headings.list)) { 
          unlisted <- unlist(mesh.headings.list[[i]]) # unlist the focal entry
          if(unlisted[grep("MajorTopicYN", names(unlisted))[1]] == "Y") mesh.headings[i] <- as.character(unlisted[grep("DescriptorName.text", names(unlisted))]) # add to the holder if "MajorTopicYN" == "Y"
        }
        mesh.headings <- mesh.headings[!is.na(mesh.headings)] # get rid of the NAs, if present
        if(length(mesh.headings) == 0) mesh.headings <- NA # Fill in NA if that got rid of everything
        else if (length(mesh.headings) > 1) mesh.headings <- paste0(mesh.headings, collapse = "_") # if multiples, paste them together with underscores
        
        output$mesh.headings <- paste("f5", as.character(mesh.headings), sep ="___")
      }}
    
    try(broad.journal.headings <- paste0(as.character(unlist(xml.item$BroadJournalHeading), collapse = "_")))
    if(exists("broad.journal.headings")) {
      if(length(broad.journal.headings) == 1){
        output$broad.journal.headings <- paste("f6", broad.journal.headings, sep ="___")
      }}
    
    try(lang.local <- xml.item$Language$text[1])
    if(exists("lang.local")) {
      if(length(lang.local) == 1){
        output$language <- paste("f7", as.character(lang.local), sep ="___")
      }}
    
    try(country.local <- xml.item$PublicationInfo$Country[1])
    if(exists("country.local")) {
      if(length(country.local) == 1){
        output$country <- paste("f8", as.character(country.local), sep ="___")
      }}
    
    # If the last character of the journal name is a space (often is), trim it off. 
    set.with.spaces <- which(substr(output$short.title, nchar.name, nchar.name) == " ")
    nchar.name <- nchar(output$short.title[set.with.spaces])
    output$short.title[set.with.spaces] <- substr(output$short.title[set.with.spaces], 1, nchar.name - 1) 
    
    return(output)
  }
  
  # This function is called on each individual journal in the XML file. It gets out the data we want, in an ugly text string. I separate fields with two tildas for easy cutting later.
  branch.fun.journals <- function(x) print(paste0(get.journal.info(xmlToList(x)),collapse = "~~"))
  
  # Run the branch function on every paper in the XML file to produce a rough and ready synopsis of each paper. Use capture.output() + print() as this seems the easiest way to deal with the XML package
  rough.data <- capture.output(xmlEventParse(file = filepath, 
                                             handlers = NULL, 
                                             branches = list(Serial = branch.fun.journals)))
  
  # Internal function that cleans up the ugly data into a sensible dataframe, ready to be returned by the main function.
  clean.rough.data <- function(rough.data){
    cleaned <- substr(rough.data, 1, nchar(rough.data)-1)
    cleaned <- substr(cleaned, 6, nchar(cleaned))
    cleaned <- strsplit(cleaned, split = "~~")
    cleaned <- do.call("rbind", lapply(cleaned, function(x){
      output <- data.frame(NLMid=NA, title=NA, short.title=NA, ISSN=NA, mesh.headings=NA, broad.journal.headings=NA, language=NA, country=NA)
      
      # I chose to add these text tags (e.g. f1___ means field one)
      # so that there was no chance of mixing up different fields
      # The tidy function looks for the right tag, then trims it off
      tidy <- function(item, tag){ 
        focal <- item[grep(paste(tag, "___", sep=""),item)]
        if(length(focal) == 1) return(gsub(paste(tag, "___", sep=""), "", focal))
        else return(NA)
      }
      
      output$NLMid <- tidy(x, "f1")
      output$title <- tidy(x, "f2")
      output$short.title <- tidy(x, "f3")
      output$ISSN <- tidy(x, "f4")
      output$mesh.headings <- tidy(x, "f5")
      output$broad.journal.headings <- tidy(x, "f6")
      output$language <- tidy(x, "f7")
      output$country <- tidy(x, "f8")
      return(output)
    }))
    
    return(cleaned)
  }
  return(clean.rough.data(rough.data))
}



# Functioned used to get all the files produced by pase.many.xml.files() into a single massive spreadsheet
zip.all.files.together <- function(file.list = c(list.files("../data/Main files - parsed", full.names = T), list.files("../data/Update files - parsed", full.names = T))){
  
  # make sure the files are opened in the right order - newest update files, then oldest update files, then newest main files, and finally oldest main files.
  file.list <- c(rev(sort(file.list[grep("Update", file.list)])), rev(sort(file.list[grep("Main", file.list)])))
  
  # Define function that opens a file and returns only papers for which we have at least one author whose forename is greater than 1 character long
  discard.all.initial.only.papers <- function(file){
    df <- read.csv(file, stringsAsFactors = F) # open the file
    df <- df[!is.na(df$forenames), ] # get rid of papers with no authors
    split.names <- strsplit(df$forenames, split="_") # split the forenames
    to.discard <- sapply(split.names, function(x){ # Make a logical vector of the files to discard
      if(length(x) == 0) return(TRUE)
      if(max(nchar(x)) == 1) return(TRUE)
      return(FALSE)
      })
    return(df[!to.discard, ]) 
  }
  
  dummy <- read.csv(file.list[1], stringsAsFactors = F, nrows = 1)
  mega.file <- as.data.frame(matrix(NA, ncol = ncol(dummy), nrow = 30000000)) # pr-allocate a massive spreadsheet to store the data. Assumes we have less than 30m papers :)
  names(mega.file) <- names(dummy)
  for(i in 1:ncol(mega.file)) mega.file[,i] <- as.character(mega.file[,i])
  start <- 1
  print(paste("Doing file", 1, "out of", length(file.list)))
  
  for(i in 1:length(file.list)){
    if(i %% 10 == 0) { # print a report every 10 files
      print(paste("Doing file", i, "out of", length(file.list)))
      print(paste("Processed", sum(!is.na(mega.file[,1])), "papers so far"))
    }
    focal <- discard.all.initial.only.papers(file.list[i])
    nrow.focal <- nrow(focal)
    if(nrow.focal > 0){
      mega.file[start:(start + nrow.focal - 1), ] <- focal
      start <-  start + nrow.focal
    }
  }
  mega.file <- mega.file[!is.na(mega.file$pmid), ] # trim off the empty space at the end of the file
  mega.file <- mega.file[with(mega.file, order(journal, pmid)), ] # re-order by journal and PMID
  return(mega.file)
}




 
# First, we need to make a list of all the journals, plus list the files that they are in. Since this takes a while, we write a file with this information for use in the next function.
# make.journal.key.file <- function(file.list = c(list.files("../data/Main files - parsed", full.names = T), list.files("../data/Update files - parsed", full.names = T))){
#   journal.key <- data.frame(journal = rep(NA, 20000000), file.plus.rows = NA)
#   start <- 1
#   tally <- 0
#   num.files <- length(file.list)
#   print(paste("Doing file", 1, "out of", num.files))
#   for(i in 1:num.files){
#     if(i %% 10 == 0) print(paste("Doing file", i, "out of", num.files))
#     journal.column <- read.csv(file.list[i], stringsAsFactors = F)$journal # open the focal file and get the 'journal' column
#     journal.column <- journal.column[!is.na(journal.column)] # Get rid of the NAs
#     unique.journals <- unique(journal.column)
#     num.rows <- length(unique.journals)
#     cells <- start:(start+num.rows-1)
#     journal.key$journal[cells] <- unique.journals
#     journal.key$file.plus.rows[cells] <- paste(file.list[i], lapply(unique.journals, function(x) paste0(range(which(journal.column == x)), collapse = "-")), sep = "_")
#     start <- start + num.rows
#     tally <- tally + num.rows
#     if(tally > 20000000) print("uh-oh, you need to declare a bigger empty in the make.journal.key.file() function")
#   }
#   
#   journal.key <- journal.key[!is.na(journal.key$journal), ] # remove the excess of NAs at the end of the file
#   journal.key <- journal.key[order(journal.key$journal), ] # sort by journal
#   return(journal.key)
# }
# 
# 
# # This function takes all the parsed files, groups them into a dataframe by journal, removes duplicated PMIDs (using the updates one where available), and finally writes a new csv file to disk (one per journal) in the specified directory
# group.by.journal <- function(output.dir, journal.key, focal.chunk = NULL, num.chunks = 8){
#   
#   unique.journals <- sort(unique(journal.key$journal)) # Here are all the unique journal found earlier
#   num.journals <- length(unique.journals)
#   
#   if(!is.null(focal.chunk)){ # Specify the focal chunk of files
#     chunk.size <- floor(num.journals / num.chunks)
#     if(focal.chunk > num.chunks) return("Check your inputs, they're wrong!")
#     
#     # Make a sequence of journals. Spread evenly over the available ones
#     sequence <- seq(from = focal.chunk, to = (num.chunks*chunk.size)+focal.chunk-1, by = num.chunks)
#     #if final chunk, add on the bonus files
#     if(focal.chunk == num.chunks) sequence <- c(sequence, (1:num.journals)[(1:num.journals) > max(sequence)])
# 
#     unique.journals <- unique.journals[sequence] # Restrict to just the focal chunk of journals
#     num.journals <- length(unique.journals) # re-count the number of journals in the subset
#   }
#   
#   
#   column.names <- names(read.csv(list.files("../data/Main files - parsed/", full.names = T)[1],nrows =1))
#   
#   file.opener <- function(filepath.plus.row) { # Internal function that can open the combination of filepath and row names that is passed to it in the journal.key file 
#     split <- strsplit(filepath.plus.row, split = "_")[[1]]
#     split.rows <- as.numeric(strsplit(split[2], split = "-")[[1]])
#     focal <- read.csv(file = split[1], 
#                     skip = split.rows[1] - 1, 
#                     nrows = 1 + split.rows[2] - split.rows[1],
#                     stringsAsFactors = F)
#     names(focal) <- column.names
#     return(focal)
#   }
#   
#   print(paste("Doing journal", 1, "of", num.journals))
#   for(i in 1:num.journals){
#     if(i %% 50 == 0) print(paste("Doing journal", i, "of", num.journals))
#     files.to.open <- journal.key$file.plus.rows[journal.key$journal == unique.journals[i]]
#     # make sure they're in right order - newest update files, then oldest update files, then newest main files, and finally oldest main files.
#     files.to.open <- c(rev(sort(files.to.open[grep("Update", files.to.open)])), rev(sort(files.to.open[grep("Main", files.to.open)])))
#     data.for.focal.journal <- do.call("rbind", lapply(files.to.open, file.opener)) # Open each file and get the rows that have the focal journal in them
#     data.for.focal.journal <- data.for.focal.journal[!(duplicated(data.for.focal.journal$pmid)), ] # Remove any duplicate papers (will use the newest one, if there is a Pubmed update)
#     
#     focal.journal <- unique.journals[i] # Now write an output file, named after the journal
#     nchar.foc <- nchar(focal.journal)
#     if(substr(focal.journal, nchar.foc, nchar.foc) == " ") focal.journal <- substr(focal.journal, 1, nchar.foc - 1) # If the last character of the journal name is a space (often is), trim it off. 
#     write.csv(data.for.focal.journal, file = paste(output.dir, "/", focal.journal, ".csv", sep=""), row.names = F)
#   }
# }
# 
# 
# 
# # Discard journals that have fewer than 'sample.size.cutoff' papers that meet the following requirements:
# # - Have at least one author with a forename listed
# # - Are a 'journal article'
# ####### - Have an Abstract (i.e. character count of abstract is > 0)
# cull.some.papers.and.journals <- function(input.dir, output.dir, sample.size.cutoff){
#   files <- list.files(input.dir, full.names = T)
#   journals <- gsub(".csv","",sapply(strsplit(files, split="/"),tail,1))
#   sample.size.list <- data.frame(journal = journals, num.papers=NA, num.with.authors=NA)
#   for(i in 1:length(files)){
#     
#     if(exists("focal")) rm(focal)
#     try(focal <- read.csv(files[i], stringsAsFactors = F), silent=T)
#     if(!exists("focal"))  sample.size.list[i, 2:3] <- 0 # If there are no papers in the file, just record that and move on to the next file
#     
#     else{ # Otherwise, apply selection criteria, save the file, and record the sample size for that journal
#       # Restrict to journal articles with abstracts, and discard those columns
#       # focal <- focal[focal$abstract.length > 0 & focal$is.journal.article == 1, !(names(focal) %in% c("abstract.length", "is.journal.article"))]
#       focal$forenames <- as.character(focal$forenames) # R will interpret the forename F or T as logical, so strsplit gets confused :)
#       focal <- focal[focal$is.journal.article == 1, !(names(focal) %in% c("is.journal.article"))] # new line of code that just checks JA status
#       focal <- focal[!is.na(focal$forenames), ] # Remove papers with no authors at all
#       focal <- focal[!is.na(focal$journal), ] # Remove papers with no known journal attached (should have been done already but no harm double checking)
#       
#       sample.size.list$num.papers[i] <- nrow(focal) # Number of papers, including those with initials only
#       
#       if(sample.size.list$num.papers[i] == 0) sample.size.list$num.with.authors[i] <- 0
#       else{
#         papers.with.a.name <- sapply(strsplit(focal$forenames, split="_"), function(x) max(nchar(x))) > 1
#         sample.size.list$num.with.authors[i] <- sum(papers.with.a.name)
#         # If there are enough papers with names, trim the initial-only papers off and write a file to disk containing the culled dataset for this particular journal
#         if(sample.size.list$num.with.authors[i] > sample.size.cutoff) write.csv(focal[papers.with.a.name, ], row.names = F, file = paste(output.dir, "/", journals[i], ".csv", sep=""))
#       }
#     }
#   }
#   # Save a list of the number of papers per journal at this stage of culling
#   write.csv(sample.size.list, file = "../outputs/Sample size at culling stage 1.csv",row.names = F) 
# }
# 

# 
# # Function to make one combined spreadsheet with only forename papers from correct journals
# make.combined.datasheet <- function(input.dir, output.file){
#   combo <- do.call("rbind", lapply(input.dir, read.csv, stringsAsFactors = F))
#   print(paste("Number of duplicate pmids:", sum(duplicated(combo$pmid)))) # should be zero
#   combo <- focal[order(focal$journal, convert.dates(focal$date), focal$pmid), ] # Reorder by journal, date, and PMID
#   write.csv(combo, file = output.file, row.names = F)
# }


########################################################################################################
# TRY TO ASSIGN A GENDER TO EVERY AUTHOR IN THE DATASET
########################################################################################################

# Functions that assign a gender based on first names, using the 'gender' and 'genderdata' packages. See here for info: https://cran.rstudio.com/web/packages/gender/vignettes/predicting-gender.html

# This first function attempts to assign gender to every unique name that was recovered, and makes a big list of names and genders (this list can then be fed into the function "genderize.paper.data")
# Arg: takes every author list that I have retrieved (i.e. names separated by underscores), which are assumed to be housed in a data frame with a column called 'forenames'
# My function uses the "ssa" method for the gender function; this method uses name-gender associations from USA census data from 1930 to 2012.
make.gender.table <- function(df, male.cutoff = 0.95, female.cutoff = 0.05, return.cutoff.only = T){
  all.names <- unlist(strsplit(df$forenames, split = "_")) # Remove the underscores
  authors <- all.names
  number.of.authors <- length(authors)
  all.names <- all.names[nchar(all.names) > 1] # Exclude all the single-letter names, which are initials and cannot be genderized
  number.of.authors.with.forenames <- length(all.names)
  all.names <- unique(all.names) # Get the unique names
  number.of.unique.names <- length(all.names)
  gender.table <- gender(all.names, method = "ssa")
  output <- data.frame(name = all.names, gender = "U", prop.male = 0, stringsAsFactors = F)
  output <- data.frame(name = gender.table$name, gender = "U", prop.male = gender.table$proportion_male, stringsAsFactors = F)
  output$gender[gender.table$proportion_male > male.cutoff] <- "M"
  output$gender[gender.table$proportion_male < female.cutoff] <- "F"
  
  number.of.names.confidently.gendered <- sum(gender.table$gender != "U")
  number.of.authors.confidently.gendered <- sum(authors %in% output$name[output$gender != "U"])
  
  print(paste("Number of authors:", number.of.authors))
  print(paste("Number of authors with forenames:", number.of.authors.with.forenames))
  print(paste("Number of unique names:", number.of.unique.names))
  print(paste("Number of names confidently gendered:", number.of.names.confidently.gendered))
  print(paste("Number of authors confidently gendered:", number.of.authors.confidently.gendered))
  print(paste("% authors confidently gendered:", round(100 * number.of.authors.confidently.gendered / number.of.authors, 2)))
  
  if(return.cutoff.only) return(output[output$gender != "U", 1:2])
  
  allnames <- union(output$name, all.names)
  excess <- length(allnames) - nrow(output)
  output <- data.frame(name = allnames, gender = c(output$gender, rep("U",excess)), prop.male =  c(output$prop.male, rep(NA,excess)), stringsAsFactors = F)
  return(output)
}


# Args: the huge dataframe, and a table of all the known genders for names in the dataset (made with the previous function)
genderize.df <- function(df, gender.table){
  df$gender <- "U"
  df$gender.score <- NA
  if(nrow(df) > 15000){
    chunk.size <- 15000 # higher numbers cause memory problems
    nPapers <- nrow(df)
    nChunks <- ceiling(nPapers / chunk.size)
    starts <- seq(1,1+(nChunks-1)*chunk.size, by = chunk.size)
    ends <- seq(chunk.size, nChunks*chunk.size, by = chunk.size); ends[length(ends)] <- nPapers
    for(i in 1:length(starts)){
      print(paste("Doing chunk ", i, " of ", length(starts), ".", sep=""))
      split <- strsplit(df$forenames[starts[i]:ends[i]], split = "_")
      names(split) <- df$pmid[starts[i]:ends[i]]
      author.lists <- melt(split) # first col ("value") is a name, second ("L1") is the rownumber of its paper
      author.lists$gender <- "U"
      author.lists$gender.score <- "U"
      
      author.lists$gender <- gender.table$gender[match(author.lists$value, gender.table$name)]
      author.lists$gender[is.na(author.lists$gender)] <- "U"
      
      author.lists$gender.score <- gender.table$prop.male[match(author.lists$value, gender.table$name)]
      
      xx <- tapply(author.lists$gender, author.lists$L1, function(x) paste0(x, collapse=""))
      df$gender[match(names(xx), df$pmid)] <- as.character(xx)
      xx <- tapply(round(author.lists$gender.score, 3), author.lists$L1, function(x) paste0(x, collapse="_"))
      df$gender.score[match(names(xx), df$pmid)] <- as.character(xx)
    }
  }
  
  else author.lists$gender<- gender.table$gender[match(author.lists$value, gender.table$name)]
  
  return(df)
}


########################################################################################################
# FUNCTIONS TO GET AUTHOR AFFILIATIONS FROM THE PUBMED XML FILES (I added these after getting all the other data -
# ideally you'd get the addresses at the same time as the author list, date, etc. But it took me a while to work
# out how to parse addresses effectively)
########################################################################################################

# This function parses one XML file and attempts to retrieve the "house", city, state, and country for each listed address on each paper using parse_addr() function from the poster library
# When the parse_addr function cannot do it, it tends to dump all of the address in the "house" field, so we should save that too.
# "filepath" argument - the filepath to the XML file to be parsed.
parse.xml.affiliations <- function(filepath){
  
  address.parse <- function(addresses){
    parsed <- parse_addr(addresses)
    with(parsed, data.frame(house = house, city = city, state = state, country = country))
  }
  
  get.affiliation <- function(x){
    output <- data.frame(pmid=NA, house=NA, city=NA, state=NA, country=NA)
    
    pmid.local <- try(x$PMID$text) 
    if(exists("pmid.local")) {
      if(length(pmid.local) == 1){
        output$pmid <- paste("f1", pmid.local, sep ="___")
      }}
    
    x <- unlist(x)
    field.names <- names(x)   # We look for any field that has "Author" AND "Affiliation" in its name
    affils.local <- try(as.character(x[intersect(grep("Author", field.names), grep("Affiliation", field.names))]))
    affils.local <- try(affils.local[nchar(affils.local) > 0]) # Often the address field is present but has nothing in it
    
    if(exists("affils.local")) {
      if(length(affils.local) >= 1){
        parsed.addresses <- address.parse(affils.local)
        output$house <- paste("f2", paste0(parsed.addresses$house, collapse = "__"), sep ="___")
        output$city <- paste("f3", paste0(parsed.addresses$city, collapse = "__"), sep ="___")
        output$state <- paste("f4", paste0(parsed.addresses$state, collapse = "__"), sep ="___")
        output$country <- paste("f5", paste0(parsed.addresses$country, collapse = "__"), sep ="___")
      }}
    return(output)
  }
  
  # This function is called on each individual paper in the XML file. It gets out the data we want, in an ugly text string. I separate fields with two tildas for easy cutting later.
  branch.fun <- function(x) print(paste0(get.affiliation(xmlToList(x)),collapse = "~~"))
  
  # Run the branch function on every paper in the XML file to produce a rough and ready synopsis of each paper. Use capture.output() + print() as this seems the easiest way to deal with the XML package (handlers?? WTF?)
  rough.data <- capture.output(xmlEventParse(file = filepath, 
                                             handlers = NULL, 
                                             branches = list(MedlineCitation = branch.fun)))
  # Internal function that cleans up the ugly data into a sensible dataframe, ready to be returned by the main function.
  clean.rough.data <- function(rough.data){
    cleaned <- substr(rough.data, 1, nchar(rough.data)-1)
    cleaned <- substr(cleaned, 6, nchar(cleaned))
    cleaned <- strsplit(cleaned, split = "~~")
    cleaned <- do.call("rbind", lapply(cleaned, function(x){
      
      # I chose to add these text tags (e.g. f1___ means field one)
      # so that there was no chance of mixing up different fields
      # The tidy function looks for the right tag, then trims it off
      tidy <- function(item, tag){ 
        focal <- item[grep(paste(tag, "___", sep=""),item)]
        if(length(focal) == 1) return(gsub(paste(tag, "___", sep=""), "", focal))
        else return(NA)
      }
      data.frame(pmid=tidy(x, "f1"), house = tidy(x, "f2"), city = tidy(x, "f3"), state = tidy(x, "f4"), country = tidy(x, "f5"))
    }))
    cleaned
  }
  output <- clean.rough.data(rough.data)
  for(i in 2:5) output[, i] <- gsub("[.]", "", output[, i]) # remove full stops from addresses
  
  output <- output[!is.na(output$pmid), ] # Now remove the rows that we cannot use (saves space on disk)
  collapsed <- apply(output[,2:5], 1, paste0, collapse="") 
  test.NA <-  nchar(collapsed) - (str_count(collapsed, "NA")*2 + str_count(collapsed, "_"))
  output[test.NA > 0, ]
}


parse.many.xml.affiliations <- function(input.dir, output.dir, focal.chunk=NULL, num.chunks = 8, over.write = F){
  xml.files <- list.files(input.dir, full.names = T) # list all the files in the input
  xml.files <- xml.files[grep("medline",xml.files)] # Make sure they are all Medline files
  num.files <- length(xml.files) # count them
  if(!is.null(focal.chunk)){ # Specify the focal chunk of files
    chunk.size <- floor(num.files / num.chunks)
    if(focal.chunk > num.chunks) return("Check your inputs, they're wrong!")
    
    # Make a sequence of files. Spread evenly over the available ones
    sequence <- seq(from = focal.chunk, to = (num.chunks*chunk.size)+focal.chunk-1, by = num.chunks)
    #if final chunk, add on the bonus files
    if(focal.chunk == num.chunks) sequence <- c(sequence, (1:num.files)[(1:num.files) > max(sequence)])
    
    xml.files <- xml.files[sequence]
    num.files <- length(xml.files) # re-count the number of files in the subset
  }
  
  if(!over.write) completed.files <- list.files(output.dir, full.names = T) # list of file paths that were processed already
  # Loop over the focal files and do them in order, following overwrite rules. Write parsed files to disk in the output directory
  for(i in 1:num.files){ 
    print(paste("Doing file ", i, " out of ", num.files, ".", sep=""))
    output.file.path <- paste(output.dir, "/", gsub("[.]xml[.]gz","",tail(strsplit(xml.files[i], split = "/")[[1]],1)), ".csv", sep="")
    if(!(output.file.path %in% completed.files) | over.write) write.csv(parse.xml.affiliations(xml.files[i]), file = output.file.path, row.names = F)
  }
}
 


# Function used to see which country names appear frequently, in order to find the right way to parse the data
get.unique.countries <- function(dirs){
  require(reshape2)
  files <- do.call("c", lapply(dirs, list.files, full.names = T))
  countries <- vector(mode = "list", length = length(files))
  for(i in 1:length(files)){
    if(i %% 10 == 0) print(paste("Doing ", i, " out of ", length(files), ".", sep = ""))
    focal <- read.csv(files[i], stringsAsFactors = F)
    if(nrow(focal) > 0 & is.character(focal$country)){
      focal <-  unlist(strsplit(focal$country, split = "__"))
      countries[[i]] <- melt(table(focal[!is.na(focal)]))
    }
  }
  countries <- do.call("rbind", countries);  names(countries) <- c("country", "count")
  countries <- melt(with(countries, tapply(count, country, sum)));   names(countries) <- c("country", "count")
  countries <- countries[!is.na(countries$country), ]
  countries$country <- as.character(countries$country)
  countries
}


# flexible function that can either take a vector of un-parsed colony names (and returns them parsed), or takes the output of get.unique.countries (in which case it returns a neat summary of how common each country is in the dataset, which we need to make 'acceptable.countries')
parse.countries <- function(countries, type = "vector"){
  countries <- country.parsing(countries, type = type)
  if(type == "vector") {
    countries[!(countries %in% acceptable.countries)] <- NA
    return(countries)
  }
  
  if(type == "summary"){
    to.remove <- c("and","the","frg","are","che","ind","box","hi-tech","cantos","republica","republic","republic of","electronic","universitaet","africa","sylvester","south")
    countries <- countries[nchar(countries$country) > 2, ] # remove abbreviations
    countries <- countries[!(countries$country %in% to.remove), ]
    countries <- melt(with(countries, tapply(count, country, sum)))
    names(countries) <- c("country", "count")
    countries <-countries[order(countries$count), ]
    countries$country <- as.character(countries$country)
    return(countries)
  }
}



# Tedious function that cleans up country names in standardised format. Used inside 'parse.countries'
country.parsing <- function(countries, type){
  if(type == "vector"){
    replace <- function(original, new) { # replaces exact matches with new term
      countries[countries == original] <- new
      countries}
    grep.term <- function(original, new){ # does a grep search for a character string, replaces with new term
      countries[grep(original, countries)] <- new
      countries}
    replace.foreign <- function(original, new, number.stars = 1){ # replace function that works on accented characters, entered as a wildcard
      countries[intersect(grep(glob2rx(original), countries), which(nchar(countries) == (nchar(original)+number.stars)))] <- new
      countries}
    remove.item <- function(to.remove){       # Simply removes the specified item
      countries[countries == to.remove] <- NA 
      countries}
  }
  else if(type == "summary"){
    replace <- function(original, new) { # replaces exact matches with new term
      countries$country[countries$country == original] <- new
      countries}
    grep.term <- function(original, new){ # does a grep search for a character string, replaces with new term
      countries$country[grep(original, countries$country)] <- new
      countries}
    replace.foreign <- function(original, new, number.stars = 1){ # replace function that works on accented characters, entered as a wildcard
      countries$country[intersect(grep(glob2rx(original), countries$country), which(nchar(countries$country) == (nchar(original)+number.stars)))] <- new 
      countries}
    remove.item <- function(to.remove)  {
      countries$country[countries$country == to.remove] <- NA # Simply removes the specified item
      countries}
  }
  
  countries <- replace("united states of america", "usa")
  countries <- replace("united states", "usa")
  countries <- replace("united states air", "usa")
  countries <- replace("united states electronic", "usa")
  countries <- replace("states", "usa")
  countries <- replace("estados", "usa")
  countries <- replace("estados unidos", "usa")
  countries <- replace("federal", "usa")
  countries <- replace("usa usa", "usa")
  countries <- replace("maine", "usa")
  countries <- replace("uk", "united kingdom")
  countries <- replace("reino unido", "united kingdom")
  countries <- replace("england", "united kingdom")
  countries <- replace("scotland", "united kingdom")
  countries <- replace("wales", "united kingdom")
  countries <- replace("northern ireland", "united kingdom")
  countries <- replace("espana", "spain")
  countries <- replace("espagne", "spain")
  countries <- replace("ardoz", "spain")
  countries <- replace("french guiana france", "france")
  countries <- replace("mex", "mexico")
  countries <- replace("deutschland", "germany")
  countries <- replace("fr germany", "germany")
  countries <- replace("russian", "russia")
  countries <- replace("russian federation", "russia")
  countries <- replace("italia", "italy")
  countries <- replace("saudi", "saudi arabia")
  countries <- replace("kingdom of saudi arabia", "saudi arabia")
  countries <- replace("kingdom of saudi", "saudi arabia")
  countries <- replace("saudi arabia department of", "saudi arabia")
  countries <- replace("italia", "italy")
  countries <- replace("slovak republic", "slovakia")
  countries <- replace("czech", "czech republic")
  countries <- replace("burkina", "burkina faso")
  countries <- replace("denmark electronic", "denmark")
  countries <- replace("new zealand electronic", "new zealand")
  countries <- replace("rs brazil", "brazil")
  countries <- replace("new", "new zealand")
  countries <- replace("new zealand and", "new zealand")
  countries <- replace("zealand", "new zealand")
  countries <- replace("of china", "china")
  countries <- replace("shandong", "china")
  countries <- replace("korea", "south korea")
  countries <- replace("republic of korea", "south korea")  
  countries <- replace("singapore electronic", "singapore")
  countries <- replace("republic of singapore", "singapore")
  countries <- replace("suisse", "switzerland")
  countries <- replace("sierra", "sierra leone")
  countries <- replace("tunisie", "tunisia")
  countries <- replace("can", "canada")
  countries <- replace("france pf", "france")
  countries <- replace("republic of ireland", "ireland")
  countries <- replace("puerto", "puerto rico")
  countries <- replace("puerto rico usa", "puerto rico")
  countries <- replace("viet", "vietnam")
  countries <- replace("belgique", "belgium")
  countries <- replace("oesterreich", "austria")
  countries <- replace("polska", "poland")
  countries <- replace("poland electronic", "poland")
  countries <- replace("canada ees", "canada")
  countries <- replace("ro brazil", "brazil")
  countries <- replace("es brazil", "brazil")
  countries <- replace("juan puerto rico", "puerto rico")
  countries <- replace("huazhong", "china")
  countries <- replace("czech republic mb", "czech republic")
  countries <- replace("arab", "united arab emirates")
  countries <- replace("slovak", "slovakia")
  countries <- replace("maroc", "morocco")
  countries <- replace("schweiz", "switzerland")
  countries <- replace("french", "france")
  countries <- replace("fin", "finland")
  countries <- replace("arabia", "saudi arabia")
  countries <- replace("rica", "costa rica")
  countries <- replace("karnataka", "india")
  countries <- replace("papua", "papua new guinea")
  countries <- replace("people's", "china") # Mostly China...
  countries <- replace("people's republic", "china") 
  countries <- replace("peoples", "china") 
  countries <- replace("people's republic of", "china") 
  countries <- replace("peoples republic of", "china")
  countries <- replace("sri", "sri lanka") 
  
  countries <- grep.term("united kingdom", "united kingdom")
  countries <- grep.term("united states", "usa")
  countries <- grep.term("usa", "usa")
  countries <- grep.term("bosnia", "bosnia and herzegovina")
  countries <- grep.term("finland", "finland")
  countries <- grep.term("czech", "czech republic")
  countries <- grep.term("brasil", "brazil")
  countries <- grep.term("brazil", "brazil")
  countries <- grep.term("romania", "romania")
  countries <- grep.term("lucia", "saint lucia")
  countries <- grep.term("denmark", "denmark")
  countries <- grep.term("zealand", "new zealand")
  countries <- grep.term("russia", "russia")
  countries <- grep.term("italy", "italy")
  countries <- grep.term("imbaro", "italy")
  countries <- grep.term("emilia", "italy")
  countries <- grep.term("spain", "spain")
  countries <- grep.term("egypt", "egypt")
  countries <- grep.term("china", "china")
  countries <- grep.term("korea", "south korea")
  countries <- grep.term("hong kong", "hong kong")
  countries <- grep.term("nederland", "netherlands")
  countries <- grep.term("netherlands", "netherlands")
  countries <- grep.term("trinidad", "trinidad and tobago")
  countries <- grep.term("tobago", "trinidad and tobago")
  countries <- grep.term("cameroun", "cameroon")
  countries <- grep.term("israel", "israel")
  countries <- grep.term("anglia", "united kingdom")
  countries <- grep.term("singapore", "singapore")
  countries <- grep.term("south africa", "south africa")
  countries <- grep.term("pakistan electronic", "pakistan")
  countries <- grep.term("belgium", "belgium")
  countries <- grep.term("sweden", "sweden")
  countries <- grep.term("martinique", "france")
  countries <- grep.term("hashomer", "israel")
  countries <- grep.term("switzerland", "switzerland")
  countries <- grep.term("thailand", "thailand")
  countries <- grep.term("norway", "norway")
  countries <- grep.term("india", "india")
  countries <- grep.term("canada", "canada")
  countries <- grep.term("saudi arabia", "saudi arabia")
  countries <- grep.term("hong", "hong kong")
  countries <- grep.term("kong", "hong kong")
  countries <- grep.term("guadeloupe", "guadeloupe")
  
  countries <- replace.foreign("alg*rie", "algeria")
  countries <- replace.foreign("fran*aise", "france")
  countries <- replace.foreign("r*union", "reunion island")
  countries <- replace.foreign("r*union france", "reunion island")
  countries <- replace.foreign("la r*union", "reunion island")
  countries <- replace.foreign("la r*union france", "reunion island")
  countries <- replace.foreign("belgi*", "belgium")
  countries <- replace.foreign("b*nin", "benin")
  countries <- replace.foreign("c*te", "ivory coast")
  countries <- replace.foreign("per*", "peru")
  countries <- replace.foreign("br*sil", "brazil")
  countries <- replace.foreign("espa*a", "spain")
  countries <- replace.foreign("canad*", "canada")
  countries <- replace.foreign("m*xico", "mexico")
  countries <- replace.foreign("m*xico m*xico", "mexico", number.stars = 2)
  countries <- replace.foreign("nl m*xico", "mexico")
  countries <- replace.foreign("s*n*gal", "senegal", number.stars = 2)
  countries <- replace.foreign("gr*ce", "greece")
  countries <- replace.foreign("*tats-unis", "usa")
  countries <- replace.foreign("panam*", "panama")
  countries <- replace.foreign("guanajuato m*xico", "mexico")
  countries <- replace.foreign("panam* panam*", "panama", number.stars = 2)
  countries <- replace.foreign("estados unidos de am*rica", "usa")
  countries <- replace.foreign("brésil", "brazil")
  countries <- replace.foreign("isra*l", "israel")
  countries <- replace.foreign("rep*blica", "republica")
  countries <- replace.foreign("r*publique", "republique")
  countries <- replace.foreign("rom*nia", "romania")
  
  countries <- remove.item("institut")
  countries <- remove.item("alliance")
  countries <- remove.item("department of")
  countries <- remove.item("department")
  countries <- remove.item("republique")
  countries <- remove.item("republica")
  countries <- remove.item("orleans")
  countries <- remove.item("randwick")
  countries <- remove.item("man")
  countries <- remove.item("ont")
  countries <- remove.item("penh")
  countries <- remove.item("del monte")
  countries <- remove.item("professor")
  countries <- remove.item("d-07743")
  countries <- remove.item("abd")
  countries <- remove.item("deu")
  countries <- remove.item("nus")
  countries <- remove.item("drs")
  countries <- remove.item("tfc")
  countries <- remove.item("glb")
  countries <- remove.item("brd")
  countries <- remove.item("remote")
  countries <- remove.item("cancer")
  countries <- remove.item("catarina")
  countries <- remove.item("research")
  countries <- remove.item("federation")
  countries <- remove.item("paraskevi")
  countries <- remove.item("alta")
  countries <- remove.item("giovanni rotondo")
  
  countries
}



# Gets the unique states from all the address that have no country listed, as well as the number of times they occur
# I built this function iteratively, and tried to assign countries to the states that come up most often
# The output is a list of states and the counries they are in, and the number of times each state appears among country-less addresses
# I tried to avoid assigning a country to states that exist in multiple countries (e.g. York could be part of New York that was chopped up by the postal package parser function, or it could be York, England)
get.unique.states <- function(dirs){
  require(reshape2)
  files <- do.call("c", lapply(dirs, list.files, full.names = T))
  states <- vector(mode = "list", length = length(files))
  for(i in 1:length(files)){
    if(i %% 10 == 0) print(paste("Doing ", i, " out of ", length(files), ".", sep = ""))
    focal <- read.csv(files[i], stringsAsFactors = F)
    if(nrow(focal) > 0 & is.character(focal$state) & is.character(focal$country)){
      got.no.country <- sapply(strsplit(focal$country, split = "__"), function(x) sum(is.na(x))==0 )
      focal <- unlist(strsplit(focal$state, split = "__")[got.no.country])
      states[[i]] <- melt(table(focal[!is.na(focal)]))
      # print(head(states[[i]]))
    }
  }
  states <- do.call("rbind", states);  names(states) <- c("state", "count")
  states <- melt(with(states, tapply(count, state, sum)));   names(states) <- c("state", "count")
  states$state <- as.character(states$state)
  states <- states[order(states$count), ]
  states$country <- NA
  usa.data <- read.csv("../data/USAstates.csv", stringsAsFactors = F)
  states$country[states$state %in% tolower(usa.data[,1])] <- "usa"
  states$country[states$state %in% tolower(usa.data[,2])] <- "usa"
  states$country[states$state == "district of columbia"] <- "usa"
  india.states <- tolower(read.csv("../data/indiaStates.csv", stringsAsFactors = F)[,1])
  for(i in 1:length(india.states)) states$country[grep(india.states[i], states$state)] <- "india"
  china.states <- tolower(read.csv("../data/chinaStates.csv", stringsAsFactors = F)[,1])
  for(i in 1:length(china.states)) states$country[grep(china.states[i], states$state)] <- "china"
  states$country[states$state == "people's"] <- "china"
  states$country[states$state == "kowloon"] <- "hong kong"
  states$country[states$state == "ontario"] <- "canada"
  states$country[states$state == "on"] <- "canada"
  states$country[states$state == "british columbia"] <- "canada"
  states$country[states$state == "saskatchewan"] <- "canada"
  states$country[states$state == "quebec"] <- "canada"
  states$country[states$state == "alberta"] <- "canada"
  states$country[states$state == "manitoba"] <- "canada"
  states$country[states$state == "nova scotia"] <- "canada"
  states$country[states$state == "new brunswick"] <- "canada"
  states$country[states$state == "newfoundland"] <- "canada"
  states$country[states$state == "newfoundland and labrador"] <- "canada"
  states$country[intersect(grep(glob2rx("qu*bec"), states$state), which(nchar(states$state) == 7))] <- "canada"
  states$country[intersect(grep(glob2rx("qu*bec qc"), states$state), which(nchar(states$state) == 9))] <- "canada"
  states$country[states$state == "qc"] <- "canada"
  states$country[states$state == "bc"] <- "canada"
  states$country[states$state == "musashino"] <- "japan"
  states$country[states$state == "barcelona"] <- "spain"
  states$country[states$state == "granada"] <- "spain"
  states$country[states$state == "finland"] <- "finland"
  states$country[states$state == "scotland"] <- "united kingdom"
  states$country[states$state == "england"] <- "united kingdom"
  states$country[states$state == "england uk"] <- "united kingdom"
  states$country[states$state == "wales"] <- "united kingdom"
  states$country[states$state == "northern ireland"] <- "united kingdom"
  states$country[states$state == "victoria"] <- "australia"
  states$country[states$state == "nsw"] <- "australia"
  states$country[states$state == "vic"] <- "australia"
  states$country[states$state == "qld"] <- "australia"
  states$country[states$state == "act"] <- "australia"
  states$country[states$state == "new south wales"] <- "australia"
  states$country[states$state == "south wales"] <- "australia"
  states$country[states$state == "queensland"] <- "australia"
  states$country[states$state == "northern territory"] <- "australia"
  states$country[states$state == "tasmania"] <- "australia"
  states$country[grep("australia", states$state)] <- "australia"
  states$country[grep("people's republic", states$state)] <- "china"
  states$country[grep("kuala lumpur", states$state)] <- "malaysia"
  states$country[grep("selangor", states$state)] <- "malaysia"
  states$country[grep("istanbul", states$state)] <- "turkey"
  states$country[intersect(grep(glob2rx("*stanbul"), states$state), which(nchar(states$state) == 9))] <- "turkey"
  
  states$country[grep("catalonia", states$state)] <- "spain"
  states$country[grep("saudi", states$state)] <- "saudi arabia"
  states$country[grep("groningen", states$state)] <- "netherlands"
  states$country[grep("crete", states$state)] <- "greece"
  states$country[grep("vietnam", states$state)] <- "vietnam"
  states$country[grep("portugal", states$state)] <- "portugal"
  states$country[grep("argentina", states$state)] <- "argentina"
  states$country[grep("usa", states$state)] <- "usa"
  states$country[grep("malaysia", states$state)] <- "malaysia"
  states$country[grep("china", states$state)] <- "china"
  states$country[grep("brazil", states$state)] <- "brazil"
  states$country[grep("australia", states$state)] <- "australia"
  states$country[grep("jerusalem", states$state)] <- "israel"
  states$country[grep("finland", states$state)] <- "finland"
  
  states$country[intersect(grep("new york", states$state), grep("guangdong", states$state))] <- NA # one address lists both of these!
  states$country[intersect(grep("saskatchewan", states$state), grep("guangdong", states$state))] <- NA
  
  states$country[intersect(grep(glob2rx("s*o paulo sp"), states$state), which(nchar(states$state) == 13))] <- "brazil"
  states$country[intersect(grep(glob2rx("paran*"), states$state), which(nchar(states$state) == 7))] <- "brazil"
  
  states$country[states$state == "minas gerais"] <- "brazil"
  states$country[states$state == "gerais"] <- "brazil"
  states$country[states$state == "paulo"] <- "brazil"
  states$country[states$state == "lisboa"] <- "portugal"
  states$country[states$state == "dc"] <- "usa"
  states$country[states$state == "czech"] <- "czech republic"
  states$country[states$state == "buenos aires"] <- "argentina"
  states$country[states$state == "buenos"] <- "argentina"
  
  states[!is.na(states$state), ]
}


# Gets the unique cities from all the address that have no state or country listed, as well as the number of times they occur
# I built this function iteratively, and tried to assign countries to the cities that come up most often
# The output is a list of cities and the countries they are in, and the number of times each state appears among country-less addresses
# I tried to avoid assigning a country to cities that exist in multiple countries
get.unique.cities <- function(dirs, known.states = states$state[!is.na(states$country)]){
  require(reshape2)
  files <- do.call("c", lapply(dirs, list.files, full.names = T))
  cities <- vector(mode = "list", length = length(files))
  for(i in 1:length(files)){
    if(i %% 10 == 0) print(paste("Doing ", i, " out of ", length(files), ".", sep = ""))
    focal <- read.csv(files[i], stringsAsFactors = F)
    if(nrow(focal) > 0 & is.character(focal$state) & is.character(focal$country) & is.character(focal$city)){
      got.no.country <- sapply(strsplit(focal$country, split = "__"), function(x) sum(is.na(x))==0 )
      got.no.state <- sapply(strsplit(focal$state, split = "__"), function(x) sum(is.na(x))==0 )
      got.neither <- got.no.country & got.no.state
      if(sum(got.neither) > 0){
        focal <- unlist(strsplit(focal$city, split = "__")[got.neither])
        cities[[i]] <- melt(table(focal[!is.na(focal)]))
        cities[[i]][,1] <- as.character(cities[[i]][,1])
        cities[[i]] <- cities[[i]][nchar(cities[[i]][,1]) > 2, ]
        # print(head(cities[[i]]))
      }
    }
  }
  cities <- do.call("rbind", cities);  names(cities) <- c("city", "count")
  cities <- melt(with(cities, tapply(count, city, sum)));   names(cities) <- c("city", "count")
  cities$city <- as.character(cities$city)
  cities <- cities[order(cities$count), ]
  cities <- cities[!is.na(cities$city), ]
  
  cities <- cities[cities$count > 10000,] # restrict to most common cities to save my workload
  
  # Try to guess the countries from the cities, using data I found online
  xx <- read.csv("../data/worldcities.csv", stringsAsFactors = F)
  xx <- xx[!is.na(xx$ISO.3166.1.country.code), ]
  xx <- xx[xx$ISO.3166.1.country.code != "", ]
  xx <- xx[xx$ISO.3166.1.country.code != "-11.704167", ]
  xx <- rbind(xx[xx$ISO.3166.1.country.code == "US", ], # Major countries first (this stops it assuming we mean the philadelphia in Jamaica etc) 
              xx[xx$ISO.3166.1.country.code == "GB", ],
              xx[xx$ISO.3166.1.country.code == "CN", ],
              xx[!(xx$ISO.3166.1.country.code %in% c("US","GB","CN")), ]
  )
  xx <- xx[!duplicated(xx$name), ]
  cities$country <- as.character(xx$ISO.3166.1.country.code[match(cities$city, tolower(xx$name))])
  xx <- read.csv("../data/country codes.csv", stringsAsFactors = F)
  cities$country[!is.na(cities$country)] <- tolower(xx$country[match(cities$country[!is.na(cities$country)], xx$code)])
  
  
  # Manually correct these:
  cities$country[grep("korea, south", cities$country)] <- "south korea"
  cities$country[grep("united states", cities$country)] <- "usa"
  
  cities$country[grep("new york", cities$city)] <- "usa"
  cities$country[grep("louisville", cities$city)] <- "usa"
  cities$country[grep("new orleans", cities$city)] <- "usa"
  cities$country[grep("new haven", cities$city)] <- "usa"
  cities$country[grep("newark", cities$city)] <- "usa"
  cities$country[grep("portland", cities$city)] <- "usa"
  
  cities$country[grep("new delhi", cities$city)] <- "india"
  cities$country[grep("toronto", cities$city)] <- "canada"
  cities$country[grep("edmonton", cities$city)] <- "canada"
  cities$country[grep("vancouver", cities$city)] <- "canada"
  cities$country[grep("saskatoon", cities$city)] <- "canada"
  cities$country[grep("munich", cities$city)] <- "germany"
  cities$country[grep("berlin", cities$city)] <- "germany"
  cities$country[grep("hamburg", cities$city)] <- "germany"
  cities$country[grep("hannover", cities$city)] <- "germany"
  cities$country[grep("bonn", cities$city)] <- "germany"
  cities$country[grep("cologne", cities$city)] <- "germany"
  cities$country[grep("dresden", cities$city)] <- "germany"
  cities$country[grep("jena", cities$city)] <- "germany"
  cities$country[grep("essen", cities$city)] <- "germany"
  cities$country[grep("ulm", cities$city)] <- "germany"
  cities$country[grep("heidelberg", cities$city)] <- "germany"
  cities$country[grep("montreal", cities$city)] <- "canada"
  cities$country[grep("stockholm", cities$city)] <- "sweden"
  cities$country[grep("lund", cities$city)] <- "sweden"
  cities$country[grep("gothenburg", cities$city)] <- "sweden"
  cities$country[grep("montpellier", cities$city)] <- "france"
  cities$country[grep("paris", cities$city)] <- "france"
  cities$country[grep("lyon", cities$city)] <- "france"
  cities$country[grep("grenoble", cities$city)] <- "france"
  cities$country[grep("lille", cities$city)] <- "france"
  cities$country[grep("bordeaux", cities$city)] <- "france"
  cities$country[grep("strasbourg", cities$city)] <- "france"
  cities$country[grep("villejuif", cities$city)] <- "france"
  cities$country[grep("ghent", cities$city)] <- "belgium"
  cities$country[grep("belgrade", cities$city)] <- "serbia"
  cities$country[grep("dublin", cities$city)] <- "ireland"
  cities$country[grep("copenhagen", cities$city)] <- "denmark"
  cities$country[grep("porto", cities$city)] <- "portugal"
  cities$country[grep("nijmegen", cities$city)] <- "netherlands"
  cities$country[grep("rotterdam", cities$city)] <- "netherlands"
  cities$country[grep("sheffield", cities$city)] <- "united kingdom"
  cities$country[grep("manchester", cities$city)] <- "united kingdom"
  cities$country[grep("oxford", cities$city)] <- "united kingdom"
  cities$country[cities$city == "york"] <- "united kingdom"
  cities$country[grep("edinburgh", cities$city)] <- "united kingdom"
  cities$country[grep("glasgow", cities$city)] <- "united kingdom"
  cities$country[grep("bristol", cities$city)] <- "united kingdom"
  cities$country[grep("liverpool", cities$city)] <- "united kingdom"
  cities$country[grep("southampton", cities$city)] <- "united kingdom"
  cities$country[grep("nottingham", cities$city)] <- "united kingdom"
  cities$country[grep("london", cities$city)] <- "united kingdom"
  cities$country[grep("cairo", cities$city)] <- "egypt"
  cities$country[grep("shanghai", cities$city)] <- "china"
  cities$country[grep("changsha", cities$city)] <- "china"
  cities$country[grep("taiwan", cities$city)] <- "china"
  cities$country[grep("china", cities$city)] <- "china"
  cities$country[grep("melbourne", cities$city)] <- "australia"
  cities$country[grep("perth", cities$city)] <- "australia"
  cities$country[grep("sydney", cities$city)] <- "australia"
  cities$country[grep("vienna", cities$city)] <- "austria"
  cities$country[grep("brisbane", cities$city)] <- "austria"
  cities$country[grep("oslo", cities$city)] <- "norway"
  cities$country[grep("bergen", cities$city)] <- "norway"
  cities$country[grep("milan", cities$city)] <- "italy"
  cities$country[grep("naples", cities$city)] <- "italy"
  cities$country[grep("milano", cities$city)] <- "italy"
  cities$country[grep("rome", cities$city)] <- "italy"
  cities$country[grep("turin", cities$city)] <- "italy"
  cities$country[grep("torino", cities$city)] <- "italy"
  cities$country[grep("pavia", cities$city)] <- "italy"
  cities$country[grep("florence", cities$city)] <- "italy"
  cities$country[grep("israel", cities$city)] <- "israel"
  cities$country[grep("santiago", cities$city)] <- "chile"
  cities$country[grep("budapest", cities$city)] <- "hungary"
  cities$country[grep("amsterdam", cities$city)] <- "netherlands"
  cities$country[grep("zurich", cities$city)] <- "switzerland"
  cities$country[grep("geneva", cities$city)] <- "switzerland"
  cities$country[grep("brussels", cities$city)] <- "belgium"
  cities$country[grep("ottawa", cities$city)] <- "canada"
  cities$country[grep("calgary", cities$city)] <- "canada"
  cities$country[grep("winnipeg", cities$city)] <- "canada"
  cities$country[grep("halifax", cities$city)] <- "canada"
  cities$country[grep("warsaw", cities$city)] <- "poland"
  cities$country[grep("bangkok", cities$city)] <- "thailand"
  cities$country[grep("madrid spain", cities$city)] <- "spain"
  cities$country[grep("madrid", cities$city)] <- "spain"
  cities$country[grep("valencia", cities$city)] <- "spain"
  cities$country[grep("italy", cities$city)] <- "italy"
  cities$country[grep("korea", cities$city)] <- "south korea"
  cities$country[grep("beijing", cities$city)] <- "china"
  cities$country[intersect(grep(glob2rx("s*o paulo"), cities$city), which(nchar(cities$city) == 10))] <- "brazil"
  cities$country[intersect(grep(glob2rx("montr*al"), cities$city), which(nchar(cities$city) == 9))] <- "canada"
  cities$country[grep("aukland", cities$city)] <- "new zealand"
  cities$country[grep("cambridge", cities$city)] <- NA # Could be USA or UK
  cities$country[grep("durham", cities$city)] <- NA # Could be USA or UK
  cities$country[grep("birmingham", cities$city)] <- NA # Could be USA or UK
  cities$country[grep("worcester", cities$city)] <- NA # Could be USA or UK
  cities$country[grep("republic", cities$city)] <- NA
  cities$country[grep("louis", cities$city)] <- NA
  cities$country[grep("hamilton", cities$city)] <- NA
  cities$country[grep("athens", cities$city)] <- NA
  cities$country[grep("parkville", cities$city)] <- NA
  cities$country[grep("richmond", cities$city)] <- NA
  cities$country[grep("kingston", cities$city)] <- NA
  cities$country[cities$city == "new"] <- NA
  
  
  cities[order(cities$country, cities$city),]
}



# Function used to look in the "house" field for information that could reveal the country of origin, where city, state and country fields all failed.
# This is where the parse_addr function dumps all the information when it cannot parse the address properly. It totally fails for Japan, hence my focus on Japan here.
# Returns a dataframe with two columns: PMID, and the country(s) that was inferred by this function from the house field
parse.house.field <- function(houses,  matches, matches.with.spaces){
  countries <- lapply(houses, function(complete.house) 
  {
    # For convenient coding, split 'matches' into two data.frames - one for multi-word search terms, one for single word search terms (e.g. 'new york' vs 'york')
    matches.without.spaces <- matches[-matches.with.spaces, ]
    matches.with.spaces <- matches[matches.with.spaces, ]
    
    # First analyse the complete string from the house field, and look for search terms that have 2 or more words in them (like "new york" - which we don't want to confuse with "york")
    multi.word.search.hits <- with(matches.with.spaces, search.term[str_detect(complete.house, search.term)])
    country1 <- NULL

    if(length(multi.word.search.hits) > 0) { 
      country1 <- with(matches.with.spaces, country[search.term %in% multi.word.search.hits]) # save the coutries matching the serach hits in 'country1'
      for(i in 1:length(multi.word.search.hits)) complete.house <- gsub(multi.word.search.hits[i], "", complete.house) # Now remove any search term hits so they are not searched again at the single-word level
    }
    # Now split the house field into individual words, using the spaces or hyphens (these are what is left after parsing by parse_addr - no commas etc). There are lots of "uni nebraska-lincoln" etc
    focal.house <- strsplit(complete.house, split = "[ -]")[[1]] 
    
    # for every individual word in the house field, try to do a match on the single-word search terms 
    country2 <- unlist(lapply(focal.house, function(individual.word) with(matches.without.spaces, country[search.term %in% individual.word])))
    country <- unique(c(country1, country2)) # Combine the unique search hits, and check how many there are.
    if(length(country) != 1) country <- NA  # Only keep the country if we identified only a single country from the address - this is conservative but safe. It will throw away stuff like "The US-China genome consortium" though.
    country
  })
  countries
}


# Function that takes the address of a parsed affiliation file, and returns it even more parsed! That is, it parses the country column into standard country names, and then tries to fill in the gaps using the state and city information (where available)
parse.one.address.file <- function(file,  matches, matches.with.spaces){
  focal <- read.csv(file, stringsAsFactors = F)
  if(nrow(focal) == 0 | !is.character(focal$city)| !is.character(focal$state)| !is.character(focal$country)) return(data.frame(pmid=numeric(0), country=numeric(0)))
  
  # First parse the countries, and add parsed values back to the dataframe in "parsed.country" column
  split <- melt(strsplit(focal$country, split = "__"))
  split$country <- parse.countries(as.character(split$value))
  focal$parsed.country <- as.character(with(split, tapply(country, L1, paste0, collapse = "_"))) 
  
  # For papers for which we have no country information, add country based on the state where available
  has.no.country <- which(nchar(focal$parsed.country) - (str_count(focal$parsed.country, "NA")*2 + str_count(focal$parsed.country, "_")) == 0)
  if(length(has.no.country) > 0){
    split <- melt(strsplit(focal$state[has.no.country], split = "__"))
    split$country <- state.countries$country[match(as.character(split$value), state.countries$state)] 
    focal$parsed.country[has.no.country] <- as.character(with(split, tapply(country, L1, paste0, collapse = "_"))) 
  }
  
  # For papers for which we STILL have no country information, add country based on the city, where available
  has.no.country <- which(nchar(focal$parsed.country) - (str_count(focal$parsed.country, "NA")*2 + str_count(focal$parsed.country, "_")) == 0)
  if(length(has.no.country) > 0){
    split <- melt(strsplit(focal$city[has.no.country], split = "__"))
    split$country <- city.countries$country[match(as.character(split$value), city.countries$city)] 
    focal$parsed.country[has.no.country] <- as.character(with(split, tapply(country, L1, paste0, collapse = "_"))) 
  }
  
  # If we STILLLLL have no country information, look in the 'house' field. This one is slow but stands a good chance to recover the country as well.
  has.no.country <- which(nchar(focal$parsed.country) - (str_count(focal$parsed.country, "NA")*2 + str_count(focal$parsed.country, "_")) == 0)
  if(length(has.no.country) > 0){
    split <- melt(strsplit(focal$house[has.no.country], split = "__"))
    split$country <- parse.house.field(as.character(split$value), matches, matches.with.spaces)
    focal$parsed.country[has.no.country] <- as.character(with(split, tapply(country, L1, paste0, collapse = "_"))) 
  }
  
  focal # return the parsed address data file
}

# Wrapper to call parse.one.address.file on all the files in "input.dirs"
make.address.datafile <- function(input.dirs, over.write=F, acceptable.countries, state.countries, city.countries, focal.chunk = NULL, num.chunks = 6){
  
  # First we need to make a big list of all the search terms (i.e. cities, states, country pseudonyms), and the countries they match with. 
  # Needed for the parse.house.field function, and we save time by defining it outside of the loop below
  japan.cities <- read.csv("../data/japanCities.csv", stringsAsFactors = F) %>% mutate(city = tolower(city), country = "japan")
  matches <- rbind(japan.cities,             
                   state.countries %>% rename(city = state), 
                   city.countries, 
                   data.frame(city = acceptable.countries, country = acceptable.countries, stringsAsFactors = F),
                   as.data.frame(rbind(
                     c("states", "usa"),
                     c("estados", "usa"),
                     c("federal", "usa"),
                     c("harvard", "usa"),
                     c("yale", "usa"),
                     c("johns hopkins", "usa"),
                     c("uk", "united kingdom"),
                     c("reino unido", "united kingdom"),
                     c("england", "united kingdom"),
                     c("scotland", "united kingdom"),
                     c("wales", "united kingdom"),
                     c("northern ireland", "united kingdom"),
                     c("espana", "spain"),
                     c("espagne", "spain"),
                     c("ardoz", "spain"),
                     c("brasil", "brazil"),
                     c("indian", "india"),
                     c("nihon", "japan"),
                     c("french guiana", "france"),
                     c("deutschland", "germany"),
                     c("russian", "russia"),
                     c("italia", "italy"),
                     c("saudi", "saudi arabia"),
                     c("italia", "italy"),
                     c("slovak", "slovakia"),
                     c("czech", "czech republic"),
                     c("burkina", "burkina faso"),
                     c("zealand", "new zealand"),
                     c("korea", "south korea"),
                     c("suisse", "switzerland"),
                     c("tunisie", "tunisia"),
                     c("can", "canada"),
                     c("guadeloupe", "france"),
                     c("viet", "vietnam"),
                     c("belgique", "belgium"),
                     c("oesterreich", "austria"),
                     c("polska", "poland"),
                     c("juan puerto rico", "puerto rico"),
                     c("huazhong", "china"),
                     c("czech republic mb", "czech republic"),
                     c("arab", "united arab emirates"),
                     c("slovak", "slovakia"),
                     c("maroc", "morocco"),
                     c("schweiz", "switzerland"),
                     c("french", "france"),
                     c("fin", "finland"),
                     c("arabia", "saudi arabia"),
                     c("rica", "costa rica"),
                     c("papua", "papua new guinea"))) %>% rename(city = V1, country = V2))  %>% rename(search.term = city) 
  
  # use this to check for duplicates when bug checking: matches[matches$search.term %in% matches$search.term[duplicated(matches$search.term)], ] %>% arrange(search.term)
  
  matches$country[matches$search.term == "usa"] <- "usa" # Now deal with the duplicates that map to 2 countries
  matches$country[matches$search.term == "georgia"] <- "usa" # sorry to the country Georgia, but I check and it seems most of these hits are USA
  matches$country[matches$search.term == "taiwan"] <- "taiwan"
  matches$country[matches$search.term == "new brunswick"] <- "canada"
  matches$country[matches$search.term == "guadeloupe"] <- "guadeloupe"
  matches$country[matches$search.term == "groningen"] <- "netherlands"
  matches$country[matches$search.term == "musashino"] <- "japan"
  matches <- matches[!duplicated(matches), ] # remove the duplicate search terms
  matches <- matches[nchar(matches$search.term) > 2 & matches$search.term != "uk", ] # remove the US state 2-character search terms (too easy to get false hits), but leave UK in there as it's very commonly used
  matches.with.spaces <- grep(" ", matches$search.term)
  
  # Now we're ready to go! Loop over the files and write them to the "temp" directory
  files <- do.call("c", lapply(input.dirs, list.files, full.names = T)) # List all the files to be parsed
  nFiles <-  length(files) # count em
  
  if(!is.null(focal.chunk)){ # Specify the focal chunk of files
    chunk.size <- floor(nFiles / num.chunks)
    if(focal.chunk > num.chunks) return("Check your inputs, they're wrong!")
    
    # Make a sequence of files. Spread evenly over the available ones
    sequence <- seq(from = focal.chunk, to = (num.chunks*chunk.size)+focal.chunk-1, by = num.chunks)
    #if final chunk, add on the bonus files
    if(focal.chunk == num.chunks) sequence <- c(sequence, (1:nFiles)[(1:nFiles) > max(sequence)])
    
    files <- files[sequence]
    nFiles <- length(files) # re-count the number of files in the subset
  }
  
  done.output.files <- list.files("../temp", full.names = T)     # Get a list of what's already done
  
  for(i in 1:nFiles){
    print(paste("Doing ", i, " out of ", nFiles, ".", sep = ""))
    focal.file <-  tail(strsplit(files[i], split = "/")[[1]],1)
    temp.file.path <- file.path("../temp", focal.file)
    if(!(temp.file.path %in% done.output.files) | over.write){
      focal <- parse.one.address.file(files[i], matches, matches.with.spaces)
      print(head(focal, 100))
      if(nrow(focal) > 0){
        test.NA <- c(which(is.na(focal$parsed.country)),  which(nchar(focal$parsed.country) - (str_count(focal$parsed.country, "NA")*2 + str_count(focal$parsed.country, "_")) == 0))
        focal <- focal[-test.NA, names(focal) %in% c("pmid", "parsed.country")]
        if(nrow(focal) > 0){
          names(focal) <- c("pmid", "country")
          write.csv(focal, file = temp.file.path, row.names = F)
        }
      }
    }
  }
}


#### FUNCTION ADDED POST-HOC TO GO BACK AND GET THE TITLE OF EACH PAPER FROM THE PUBMED XML FILES

parse.many.xml.titles <- function(input.dir, output.dir, focal.chunk=NULL, num.chunks = 8, over.write = F){
  
  parse.xml.titles <- function(filepath){ # Internal function that handles a single file
    
    get.title <- function(x){ # Sub-function that handles a single PMID within a single file
      
      pmid.local <- try(x$PMID$text) # Try to get the PMID
      x <- unlist(x)
      field.names <- names(x)   # Look for a field that has "ArticleTitle" in its name
      
      title.local <- try(as.character(x[grep("ArticleTitle", field.names)]))
      title.local <- try(title.local[nchar(title.local) > 0]) # Sometimes title field is present but empty
      
      if(exists("pmid.local") & exists("title.local")) { # Assuming that we got something for both of them...
        if(length(pmid.local) == 1 & length(title.local) >= 1){ # Check if exactly one PMID and exactly one title was found
          return(paste0(c(pmid.local, title.local), collapse = "~~"))
        }}
      
      return(NULL) # If previous if-statements failed, return NULL
    }
    
    # This function is called on each individual paper in the XML file. It gets out the data we want, in an ugly text string. I separate fields with two tildas for easy cutting later.
    branch.fun <- function(x) print(get.title(xmlToList(x)))
    
    # Run the branch function on every paper in the XML file to produce a rough and ready synopsis of each paper. Use capture.output() + print() as this seems the easiest way to deal with the XML package (handlers?? WTF?)
    rough.data <- capture.output(xmlEventParse(file = filepath, 
                                               handlers = NULL, 
                                               branches = list(MedlineCitation = branch.fun)))
    # Internal function that cleans up the ugly data into a sensible dataframe, ready to be returned by the main function.
    clean.rough.data <- function(rough.data){
      cleaned <- substr(rough.data, 1, nchar(rough.data)-1)
      cleaned <- substr(cleaned, 6, nchar(cleaned))
      cleaned <- as.data.frame(do.call("rbind", strsplit(cleaned, split = "~~")))
      names(cleaned) <- c("pmid", "title")
      cleaned$pmid <- as.character(cleaned$pmid)
      cleaned$title <- as.character(cleaned$title)
      cleaned
    }
    output <- clean.rough.data(rough.data)
    
    end.in.a.dot <- substr(output$title, nchar(output$title), nchar(output$title)) == "." # cut periods off the end of the titles
    output$title[end.in.a.dot] <- substr(output$title[end.in.a.dot], 1, nchar(output$title[end.in.a.dot]) - 1)
    bounded.by.brackets <- substr(output$title, 1, 1) == "[" &   # Some titles are [bounded by brackets like this] - cut them off
      substr(output$title, nchar(output$title), nchar(output$title)) == "]"
    output$title[bounded.by.brackets] <- substr(output$title[bounded.by.brackets], 2, nchar(output$title[end.in.a.dot]) - 1)
    output
  }
  
  # Now we can begin the parsing a chunk of XML files
  xml.files <- list.files(input.dir, full.names = T) # list all the files in the input
  xml.files <- xml.files[grep("medline",xml.files)] # Make sure they are all Medline files
  num.files <- length(xml.files) # count them
  
  if(!is.null(focal.chunk)){ # Specify the focal chunk of files
    chunk.size <- floor(num.files / num.chunks)
    if(focal.chunk > num.chunks) return("Check your inputs, they're wrong!")
    
    # Make a sequence of files. Spread evenly over the available ones
    sequence <- seq(from = focal.chunk, to = (num.chunks*chunk.size)+focal.chunk-1, by = num.chunks)
    #if final chunk, add on the bonus files
    if(focal.chunk == num.chunks) sequence <- c(sequence, (1:num.files)[(1:num.files) > max(sequence)])
    
    xml.files <- xml.files[sequence]
    num.files <- length(xml.files) # re-count the number of files in the subset
  }
  
  if(!over.write) completed.files <- list.files(output.dir, full.names = T) # list of file paths that were processed already
  
  # Loop over the focal files and do them in order, following overwrite rules. Write parsed files to disk in the output directory
  for(i in 1:num.files){ 
    print(paste("Doing file ", i, " out of ", num.files, ".", sep=""))
    output.file.path <- paste(output.dir, "/", gsub("[.]xml[.]gz","",tail(strsplit(xml.files[i], split = "/")[[1]],1)), ".csv", sep="")
    if(!(output.file.path %in% completed.files) | over.write) write.csv(parse.xml.titles(xml.files[i]), file = output.file.path, row.names = F)
  }
}

#### END







# Add some 95% confidence limits to every single proportion in the summary data sheets, using the default method from binom.test - this takes a while to calculate
add.CIs.to.summary <- function(df){
  df <- df %>% mutate(lCI.overall = NA, uCI.overall = NA, lCI.first = NA, uCI.first = NA, lCI.last = NA, uCI.last = NA, lCI.single = NA, uCI.single = NA)
  for(i in 1:nrow(df)) {
    try(df[i,grep("CI.overall",names(df))] <- with(df[i, ], as.numeric(binom.test(nFemales.overall, nFemales.overall+nMales.overall)$conf)), silent = T)
    try(df[i,grep("CI.first",names(df))] <- with(df[i, ], as.numeric(binom.test(nFemales.first, nFemales.first + nMales.first)$conf)), silent = T)
    try(df[i,grep("CI.last",names(df))] <- with(df[i, ], as.numeric(binom.test(nFemales.last, nFemales.last + nMales.last)$conf)), silent = T)
    try(df[i,grep("CI.single",names(df))] <- with(df[i, ], as.numeric(binom.test(nFemales.single, nFemales.single + nMales.single)$conf)), silent = T)
  }
  return(df)
}
