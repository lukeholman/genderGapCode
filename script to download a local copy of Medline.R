# Function which writes a bash script to download all the data from Medline/Pubmed, so you can manipulate it locally rather than on their servers (which takes weeks and generates lots of web traffic). Call this bash script from the terminal using the 'source' command.

#  Note that there are two types of Medline data files - normal files, and update files. Each year they add all the updates to the main files, scrap the updates, and then add more as the year progresses. The updates take precedence over duplicates in the main files, allowing people to correct errors in Medline or add new inforation. Thus, we need to keep them separate and in correct order (the newest update PMIDs are the ones to use)

# load this library
library(RCurl)

# Basically the funciton reads the names of everything that we want on the FTP server, and makes a script that can be used to download all those files from the Terminal
make.bash.script <- function(url, output.file){
  ftp_contents <- unlist(strsplit(getURL(url),split=" "))
  ftp_contents <- ftp_contents[grep("medline", ftp_contents)]
  ftp_contents <- ftp_contents[-c(grep("md5", ftp_contents),grep("stats", ftp_contents),grep("notes", ftp_contents))]
  ftp_contents <- paste("curl -O ", url, gsub("\n-rw-r--r--", "", ftp_contents),sep="")
  write.csv(ftp_contents, output.file, row.names = F)
}

make.bash.script("ftp://ftp.nlm.nih.gov/nlmdata/.medleasebaseline/gz/", "~/Pubmed xml files/Bash script to get main pubmed files.txt")
make.bash.script("ftp://ftp.nlm.nih.gov/nlmdata/.medlease/gz/", "~/Pubmed xml files/Bash script to get update pubmed files.txt")

# You will have to edit the script a tiny bit manually to make it run - open it and look for the minor errors in it that I didn't bother to fix (e.g. R gives the column of commands a title, which we don't need, and it puts double quotes around each command - you can remove these easily using search and replace).