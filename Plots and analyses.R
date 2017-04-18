setwd("~/Desktop/Gender Pubmed mining project/R scripts")

library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(grid)
library(gridExtra)
suppressPackageStartupMessages(library(dplyr))
library(stringr)
library(tidyr)
library(jsonlite)
source("Plot and analysis functions.R")


######################################################################################################
# SET UP THE DATA FOR ANALYSIS AND PLOTTING
######################################################################################################  

# Load up the database connection - this allows quick access to the whole PubMed dataset without loading the whole thing to memory
pubmed_db <- src_sqlite("../data for analysis/pubmed_gender_data.sqlite3", create = F)
pubmed_sqlite <- tbl(pubmed_db, "papers")
journals_sqlite <- tbl(pubmed_db, "journals")


######################################################################################################################################
# GET SOME BASIC INFORMATION ABOUT THE PUBMED DATASET (note: takes a moment to run, as it accesses the database)
######################################################################################################################################

# Sample size information (used to make supplementary table)
length(unique.journals) # There are 6471 journals in the dataset (each of these has at least 100 papers, as I discarded those with fewer than 100)
pubmed_sqlite %>% summarise(n = n()) # There are 9,885,392 papers in the full PubMed dataset (note that for some of these, we could not ascertain any authors' genders) 
gender.column <- as.data.frame(collect(pubmed_sqlite %>% select(gender), n = Inf))[,1] # Get the gender and gender.95 columns from the database
gender95.column <- as.data.frame(collect(pubmed_sqlite %>% select(gender.95), n = Inf))[,1]

sum(str_detect(gender.column, "[MF]")) # There are 9,514,554 papers where we know at least one author's gender with some degree of confidence
sum(str_detect(gender95.column, "[MF]")) # There are 9,154,630 papers where we know at least one author's gender with >=95% confidence
sum(nchar(gender.column)) # There are 48,716,542 authors in total
sum(nchar(gsub("U", "", gender.column))) # There are 39,226,243 authors whose gender is known
sum(nchar(gsub("U", "", gender95.column))) # There are 35,499,685 authors whose gender is known to >95% certainty
sum(str_count(gender95.column, "F")) / sum(str_count(gender95.column, "[MF]")) # The overall gender ratio of the entirety of PubMed (using the >95% certain gender data) is 35.1632 % female
rm(list = c("gender.column", "gender95.column"))


################################################################################################################
# INSPECT THE SUCCESS RATE AT WHICH GENDER WAS ASSIGNED TO AUTHOR NAMES FROM AROUND THE WORLD
# It's important that we check if the gender assignment success rate is similar across countries 
# This figure shows that gender assignment was mostly very good, but did not work well for e.g. China and Taiwan
################################################################################################################

mean.success <- with(read.csv("../outputs/genderizing.success.rate.by.country.csv"), 100 * sum(success.rate * n / sum(n)))
genderisation.rate.fig <-  ggplot(read.csv("../outputs/genderizing.success.rate.by.country.csv"), aes(x = factor(country, levels = country), y = 100 * success.rate, fill = log10(n))) + geom_hline(yintercept = mean.success, colour = "darkgrey", linetype = 2) + geom_bar(stat = "identity", colour = "black") + coord_flip() + scale_fill_gradientn(colours = brewer.pal(3, "Blues"), name = "Sample\nsize", breaks = 2:6, labels=c(expression(10^2), expression(10^3), expression(10^4), expression(10^5), expression(10^6))) + ylab("% authors to whom we could assign gender with >95% confidence") + xlab(NULL) + scale_y_continuous(expand = c(0,0), limits = c(0,103), breaks = seq(0,100,by=20))  + theme(legend.position = c(0.9,0.2)) + annotate("text", x = 4, y = mean.success+1, label = paste("Overall: ", round(mean.success, 1), "%", sep = ""), hjust = 0)
ggsave(genderisation.rate.fig, file = "../figures/genderisation.rate.fig.pdf", width = 8.5, height = 12.3)


##############################################################################################################################################
# Double-check that omitting the names whose gender score is between 0.05 and 0.95 does not consistently remove male- or female-type names
# For example, if Chris and Robin are mostly male and mostly female names respectively, and there are more Chrises than Robins, 
# we would be adding bias to the estiamtes of gender ratio by omitting the gender-ambiguous names.
# The test shows that there are roughly equal numbers of authors with a mostly-male name as there are with a mostly-female name, so we are ok
##############################################################################################################################################

genders <- pubmed_sqlite %>% select(gender, gender.score) %>% collect(n = Inf) %>% as.data.frame()
genders <- data.frame(gender = unlist(strsplit(genders$gender, split = "")), score = unlist(strsplit(genders$gender.score, split = "_")), stringsAsFactors = F)
genders$score[genders$score == "NA"] <- NA
genders$score <- as.numeric(genders$score)
genders %>% filter(score > 0.05 & score < 0.95) %>% nrow() # 3,726,558 authors eliminated by using gender.95 over gender
genders %>% filter(score < 0.05 | score > 0.95) %>% nrow() # 35,292,864 authors left in
genders %>% nrow() # 48716542 total authors 
1 - (genders %>% filter(!is.na(score)) %>% filter(score > 0.05 & score < 0.95) %>% nrow())  / (genders %>% filter(!is.na(score)) %>% nrow()) # 0.9049983 (proportion authors remaining, or those where there was at least some gender info)
1 - (genders %>% filter(score > 0.05 & score < 0.95) %>% nrow())  / (genders %>% nrow()) # 0.9235053 (proportion authors remaining, or those where there was at least some gender info)
genders %>% filter(score > 0.05 & score < 0.95) %>% ggplot(aes(x=1, y=score)) + geom_boxplot() + ylab("Femaleness of the name") + xlab(" ")
genders %>% filter(score > 0.05 & score < 0.95) %>% summarise(mean.score.of.rejected.names = mean(score), median = median(score),  se = sd(score)/sqrt(length(score)))
#   mean.score.of.rejected.names median           se
#                       0.471781   0.49 0.0001417435
rm(genders)
# The average gender score of the names that were scored as U in the gender.95 column BUT NOT the gender column is: 0.472 +- 0.00014
# This means we removed very slightly more male-ish names than we did female-ish names. Thus, the gender ratios in the paper are probably very slightly more female-biased than the "real" ones
# Thus, journals/disciplines are very slightly more male-biased than they appear in the paper


###########################################################################################################################################
# FIT SOME LOGISTIC MODELS TO FIND THE GENDER RATIO, ITS RATE OF CHANGE, AND YEARS-TO-PARITY ON PUBMED DISCIPLINES AND JOURNALS, PLUS ARXIV
# The function 'make.figure.one.dataset' fits a logistic curve to the gender-time relationship, and uses this curve to estimate the
# current and future gender ratio. It also bootstraps the data to get 95% CIs on these values (so this analysis takes some time)
###########################################################################################################################################

# List all the unique disciplines and journals in the PubMed data
unique.discs <- (pubmed_sqlite %>% select(discipline) %>% distinct() %>% collect(n = Inf) %>% as.data.frame())[,1]
unique.discs <- unique.discs[!is.na(unique.discs)]
unique.journals <- (pubmed_sqlite %>% select(journal) %>% distinct() %>% collect(n = Inf) %>% as.data.frame())[,1]

# Script to make the datasets to plot Figures 1 and 2, and get comparable info for every journal. 
# Run these pairs in 6 different cores (I just opened extra R windows via Terminal)
write.csv(make.figure.one.dataset(unique.discs[1:18], "disc", chunk.size = 10, nChunks = 100), file = "../outputs/Regressions for each discipline1.csv", row.names = F)
write.csv(make.figure.one.dataset(unique.journals[1:1078], "journal", chunk.size = 10, nChunks = 100), file = "../outputs/Regressions for each journal1.csv", row.names = F)

write.csv(make.figure.one.dataset(unique.discs[19:36], "disc", chunk.size = 10, nChunks = 100), file = "../outputs/Regressions for each discipline2.csv", row.names = F)
write.csv(make.figure.one.dataset(unique.journals[1079:2156], "journal", chunk.size = 100, nChunks = 10), file = "../outputs/Regressions for each journal2.csv", row.names = F)

write.csv(make.figure.one.dataset(unique.discs[37:54], "disc", chunk.size = 10, nChunks = 100), file = "../outputs/Regressions for each discipline3.csv", row.names = F)
write.csv(make.figure.one.dataset(unique.journals[2157:3234], "journal", chunk.size = 100, nChunks = 10), file = "../outputs/Regressions for each journal3.csv", row.names = F)

write.csv(make.figure.one.dataset(unique.discs[55:72], "disc", chunk.size = 10, nChunks = 100), file = "../outputs/Regressions for each discipline4.csv", row.names = F)
write.csv(make.figure.one.dataset(unique.journals[3235:4312], "journal", chunk.size = 100, nChunks = 10), file = "../outputs/Regressions for each journal4.csv", row.names = F)

write.csv(make.figure.one.dataset(unique.discs[73:90], "disc", chunk.size = 10, nChunks = 100), file = "../outputs/Regressions for each discipline5.csv", row.names = F)
write.csv(make.figure.one.dataset(unique.journals[4313:5390], "journal", chunk.size = 100, nChunks = 10), file = "../outputs/Regressions for each journal5.csv", row.names = F)

write.csv(make.figure.one.dataset(unique.discs[91:108], "disc", chunk.size = 10, nChunks = 100), file = "../outputs/Regressions for each discipline6.csv", row.names = F)
write.csv(make.figure.one.dataset(unique.journals[5391:6471], "journal", chunk.size = 100, nChunks = 10), file = "../outputs/Regressions for each journal6.csv", row.names = F)

disc.stats <- rbind(read.csv("../outputs/Regressions for each discipline1.csv", stringsAsFactors = F), read.csv("../outputs/Regressions for each discipline2.csv", stringsAsFactors = F), read.csv("../outputs/Regressions for each discipline3.csv", stringsAsFactors = F), read.csv("../outputs/Regressions for each discipline4.csv", stringsAsFactors = F), read.csv("../outputs/Regressions for each discipline5.csv", stringsAsFactors = F), read.csv("../outputs/Regressions for each discipline6.csv", stringsAsFactors = F))

# For some reason, several combos failed in the wrapper function "make.figure.one.dataset" even though they work fine when I queue them directly into the gender.stats() function. Let's re-run them now
partially.done <- disc.stats %>% filter(discipline %in% names(table(disc.stats$discipline)[table(disc.stats$discipline) < 4])) %>% select(1:2)
missing.entries <- do.call("rbind", (c(sapply(unique.discs, function(x) paste0(x, c("_First", "_Last", "_Single", "_Overall"))))[!(c(sapply(unique.discs, function(x) paste0(x, c("_First", "_Last", "_Single", "_Overall")))) %in% paste(disc.stats$discipline, disc.stats$position, sep = "_"))] %>% strsplit(split = "_"))) %>% as.data.frame() %>% rename(item = V1, position = V2) %>% mutate(position = tolower(position), item = as.character(item)) %>% filter(item %in% unique(partially.done$discipline)) %>% filter(!(paste(item, position) %in% paste(partially.done$discipline, partially.done$position)))
answers <- list()
for(i in 1:nrow(missing.entries)){
  foc <- gender.stats(item = missing.entries$item[i], filter.type = "disc", position = missing.entries$position[i], data.source = pubmed_sqlite, chunk.size=100, nChunks=10, plot = F, print.each.one = T, export.JSON = F, run.sim = T)
  if(!is.null(foc)) answers[[i]] <- foc
}
#write.csv(do.call("rbind", answers), file = "../outputs/Regressions for each discipline - extras.csv", row.names = F)

# Append all the files together and write them to the "data for analysis" directory
disc.stats <- rbind(read.csv("../outputs/Regressions for each discipline1.csv", stringsAsFactors = F), read.csv("../outputs/Regressions for each discipline2.csv", stringsAsFactors = F), read.csv("../outputs/Regressions for each discipline3.csv", stringsAsFactors = F), read.csv("../outputs/Regressions for each discipline4.csv", stringsAsFactors = F), read.csv("../outputs/Regressions for each discipline5.csv", stringsAsFactors = F), read.csv("../outputs/Regressions for each discipline6.csv", stringsAsFactors = F), read.csv("../outputs/Regressions for each discipline - extras.csv", stringsAsFactors = F))
# write.csv(disc.stats, file = "../data for analysis/Statistics for each discipline.csv", row.names = F)

journal.stats <- rbind(read.csv("../outputs/Regressions for each journal1.csv", stringsAsFactors = F), read.csv("../outputs/Regressions for each journal2.csv", stringsAsFactors = F), read.csv("../outputs/Regressions for each journal3.csv", stringsAsFactors = F), read.csv("../outputs/Regressions for each journal4.csv", stringsAsFactors = F), read.csv("../outputs/Regressions for each journal5.csv", stringsAsFactors = F), read.csv("../outputs/Regressions for each journal6.csv", stringsAsFactors = F), read.csv("../outputs/Regressions for each journal - extras.csv", stringsAsFactors = F))
#write.csv(journal.stats, file = "../data for analysis/Statistics for each journal.csv", row.names = F)

# Now do the ArXiv major and minor categories too, for Figure 1 and suppl material
arxiv.data <- read.csv("../data for analysis/arxiv data.csv", stringsAsFactors = F) # Load ArXiv data and change the date format to be the same as PubMed
if(substr(arxiv.data$date[1], 3, 3) != "_") arxiv.data$date <- paste(substr(arxiv.data$date, 9, 10), substr(arxiv.data$date, 6, 7), substr(arxiv.data$date, 1, 4), sep = "_") 
write.csv(make.figure.one.dataset(items = unique(arxiv.data$category), filter.type = "arxiv.cat", data.source = arxiv.data, chunk.size = 10, nChunks = 100), file = "../data for analysis/Statistics for ArXiv major categories.csv", row.names = F)

# Do the ArXiv sub-categories too
write.csv(make.figure.one.dataset(items = unique(paste(arxiv.data$category, arxiv.data$sub.category, sep = ": ")), filter.type = "arxiv.subcat", data.source = arxiv.data, chunk.size = 100, nChunks = 10), file = "../data for analysis/Statistics for ArXiv sub-categories.csv", row.names = F)


########################################################################################################
# CALCULATE THE DATA POINTS AND CURVE SHAPES NEEDED TO MAKE THE d3 WEB APP USING LOGISTIC MODELS
########################################################################################################

unique.discs <- (pubmed_sqlite %>% select(discipline) %>% distinct() %>% collect(n = Inf) %>% as.data.frame())[,1]
unique.discs <- unique.discs[!is.na(unique.discs)]
unique.journals <- (pubmed_sqlite %>% select(journal) %>% distinct() %>% collect(n = Inf) %>% as.data.frame())[,1]
unique.countries <- read.csv("../outputs/genderizing.success.rate.by.country.csv", stringsAsFactors = F) %>% filter(success.rate * n > 999)  # only bother doing countries where we measured the gender of at least 1000 authors
unique.countries <-  tolower(unique.countries$country)

# Get all possible combinations of journal/disc/country/author position, which we will need to make the data for the web app
arguments <- make.list.of.combos.for.web.app(unique.discs, unique.journals, unique.countries)
write.csv(arguments, file = "../outputs/List of arguments for making the web app.csv", row.names = F)

# Run these in 5 different cores
arguments <- read.csv("../outputs/List of arguments for making the web app.csv", stringsAsFactors = F)
starts <- seq(1, 1 + floor(nrow(arguments)/5)*5, by = floor(nrow(arguments)/5))
ends <- seq(floor(nrow(arguments)/5), floor(nrow(arguments)/5)*5, by= floor(nrow(arguments)/5))
ends[length(ends)] <- nrow(arguments)

core <- 1; arguments <- arguments[starts[core]:ends[core], ]; lapply(1:nrow(arguments), make.web.app.data, arguments=arguments)
core <- 2; arguments <- arguments[starts[core]:ends[core], ]; lapply(1:nrow(arguments), make.web.app.data, arguments=arguments)
core <- 3; arguments <- arguments[starts[core]:ends[core], ]; lapply(1:nrow(arguments), make.web.app.data, arguments=arguments)
core <- 4; arguments <- arguments[starts[core]:ends[core], ]; lapply(1:nrow(arguments), make.web.app.data, arguments=arguments)
core <- 5; arguments <- arguments[starts[core]:ends[core], ]; lapply(1:nrow(arguments), make.web.app.data, arguments=arguments)


# Double check that there are no trailing spaces in journal names, and
# that there is a discipline properly assigned, and that General Surgery is correctly named
lapply(list.files("../outputs/json files for web app/", full.names = T), function(f){
  focal <- fromJSON(f)
  nc <- nchar(focal$Journal)
  if(substr(focal$Journal, nc, nc) == " ") focal$Journal <- substr(focal$Journal, 1, nc-1)
  if(is.na(focal$Discipline)) focal$Discipline <- (journals_sqlite %>% filter(short.title == focal$Journal) %>% select(Discipline) %>% as.data.frame())[,1][1]
  if(focal$Discipline == "Surgery") focal$Discipline <- "General Surgery"
  write(toJSON(focal, pretty = F), f)
} )

# Append together all the JSON files into one, and write it to disk
write(lapply(list.files("../outputs/json files for web app/", full.names = T), function(f) fromJSON(f)) %>% toJSON(pretty=F), file = "../outputs/data_for_web_app.json")

# Make 100% sure the ones that erroneously list "Kingdom" as the country are deleted
json.data <- jsonlite::fromJSON(txt = "../outputs/data_for_web_app.json") %>% filter(Country != "Kingdom")
json.data <- json.data %>% select(Discipline, Country, Journal, Curve, Points, Position) # re-order the columns in the manner Errol asked for
for(i in c(1:3, 6)) json.data[,i] <- unlist(json.data[,i]) # Unlist these ones, to get rid of square brackets

write(toJSON(json.data, pretty = F), file = "../Errol's web app/genderGap/data_for_web_app_unparsed.json")


########################################################################################################
# MAKE FIGURES 1 AND 2
########################################################################################################

disc.stats <- clean.up.data.for.figures(rbind(read.csv("../data for analysis/Statistics for each discipline.csv", stringsAsFactors = F) %>% select(-journal) %>% mutate(origin = "PubMed"), read.csv("../data for analysis/Statistics for ArXiv major categories.csv", stringsAsFactors = F) %>% select(-sub.category) %>% rename(discipline = category) %>% mutate(origin = "ArXiv")), "disc")
sub.cats <- clean.up.data.for.figures(read.csv("../data for analysis/Statistics for ArXiv sub-categories.csv", stringsAsFactors = F) %>% mutate(discipline = sub.category), filter.type = "disc")
complete.data.set1 <- disc.stats %>% filter(discipline %in% unique(disc.stats$discipline)[1:(-1+ceiling(length(unique(disc.stats$discipline)) / 2))])
complete.data.set2 <- disc.stats %>% filter(!(discipline %in% complete.data.set1$discipline)])

# Here is a nice small summary plot. I made this in case the journal says that Figures 1 and 2 are too big
library(ggbeeswarm)
grid.arrange(
  ggplot(disc.stats, aes(x=position, y=gender.ratio.at.present, fill=gender.ratio.at.present)) + geom_hline(yintercept = 50, linetype = 2)  + geom_quasirandom(pch=21, colour = "darkgrey") + theme_classic() + theme(legend.position = "none") + xlab(NULL) + ylab("% women authors in 2016") + scale_fill_distiller(type='div', palette = 7),
  ggplot(disc.stats, aes(x=position, y=current.rate.of.change, fill=gender.ratio.at.present)) + geom_hline(yintercept = 0, linetype = 2) + geom_quasirandom(pch=21, colour = "darkgrey") + theme_classic() + theme(legend.position = "none") + xlab(NULL) + ylab("Change in % women authors over the next year") + scale_fill_distiller(type='div', palette = 7),
  ggplot(disc.stats %>% mutate(years.to.parity = replace(years.to.parity, is.na(years.to.parity), 350)), aes(x=position, y=years.to.parity, fill=gender.ratio.at.present)) + geom_hline(yintercept = 0, linetype = 2) + geom_quasirandom(pch=21, colour = "darkgrey") + theme_classic() + theme(legend.position = "none")  + xlab(NULL) + ylab("Number of years until gender parity") + scale_fill_distiller(type='div', palette = 7),
  ncol=3, bottom = "Authorship position")




# Plot figures 1 and 2, and export to PDF files
main.figure(complete.data.set1, c(7,50))
main.figure(complete.data.set2, c(10,90))
ggsave(main.figure(complete.data.set1, c(7,50)), file = "../figures/Figure 1.pdf", height = 8.4, width = 8)
ggsave(main.figure(complete.data.set2, c(10,90)), file = "../figures/Figure 2.pdf", height = 8.4, width = 8)

# Accompanying text:
length(unique(complete.data.set1$discipline)) + length(unique(complete.data.set2$discipline)) # There are 115 disciplines in the figure
sum(disc.stats$up.CI.1[disc.stats$position == "Overall"] < 45) # 87 of them are significantly <50% women
sum(disc.stats$low.CI.1[disc.stats$position == "Overall"] > 55) # 8 are >55% women

# Make equivalent supplementary figures for the ArXiv sub-categories
ggsave(main.figure(sub.cats %>% filter(category == "Physics"), c(0,40), c(-1,3), add.density3 = F), filename = "../figures/arxiv.subcats1.pdf", width = 8, height = 8)
ggsave(main.figure(sub.cats %>% filter(category == "Mathematics"), c(0,33), c(-1,3)), filename = "../figures/arxiv.subcats2.pdf", width = 8, height = 8)
ggsave(main.figure(sub.cats %>% filter(category == "Computer Science"), c(0,35), c(-2,4), add.density3 = F), filename = "../figures/arxiv.subcats3.pdf", width = 8, height = 8)
ggsave(main.figure(sub.cats %>% filter(category == "Quantitative Biology"), c(10,35), c(-1.5, 1.8), add.density = F, single.present = F), filename = "../figures/arxiv.subcats4.pdf", width = 8, height = 6)
ggsave(main.figure(sub.cats %>% filter(category == "Astrophysics"), c(0,35), c(-2, 3), add.density = F), filename = "../figures/arxiv.subcats5.pdf", width = 8, height = 4)
ggsave(main.figure(sub.cats %>% filter(category == "Statistics"), c(5, 35), add.density = F), filename = "../figures/arxiv.subcats6.pdf", width = 8, height = 4)
ggsave(main.figure(sub.cats %>% filter(category == "Nonlinear Sciences"), c(10, 25), c(-0.5, 1.5), add.density = F, single.present = F), filename = "../figures/arxiv.subcats7.pdf", width = 8, height = 4)

# Some results for the text to go with Figures 1 and 2
disc.stats %>% filter(position == "Overall" & !is.na(qualitative.result) & grepl("Essentially at parity", qualitative.result)) # 28 out of
disc.stats %>% filter(position == "Overall") %>% nrow()                                                               # 123 disciplines are currently within 5% of parity
disc.stats %>% filter(position == "Overall" & !grepl("Essentially at parity", qualitative.result) & years.to.parity <= 20) %>% nrow() # 35 are not at parity, but will get there within the next 20 years
disc.stats %>% filter(position == "Overall" & !grepl("Essentially at parity", qualitative.result) & years.to.parity > 20) %>% nrow() # 48 will reach parity in >20 years
disc.stats %>% filter(position == "Overall" & !is.na(qualitative.result) & grepl("biased", qualitative.result)) # These ones are getting more biased


##########################################################################################################
# MAKE SUPPLEMENTARY FIGURE SHOWING THAT VERY MALE-BIASED FIELDS ARE IMPROVING SLOWER THAN NEAR-EQUAL ONES
##########################################################################################################

# The relationship between current % women and the rate of change, for male-biased disciplines
hump.figure <- disc.stats %>% filter(gender.ratio.at.present < 50) %>% 
  ggplot(aes(x = gender.ratio.at.present, y= current.rate.of.change, weight = n.authors))  + geom_hline(yintercept = 0, linetype = 2) + stat_smooth(method = "gam", formula = y ~ s(x, k = 3), colour = "black") + geom_point(alpha=0.7, aes(colour = position), size = 0.8) + scale_colour_manual(values = c(brewer.pal(3, "Set1"), "black")) + xlab("% women authors in 2016") + ylab("Change in % women authors per year") + facet_wrap(~position, labeller = as_labeller(c(`First` = "First authors", `Last` = "Last authors", `Single` = "Single authors", `Overall` = "All authors"))) 
ggsave(hump.figure, file = "../figures/hump figure.pdf", height = 7, width = 7)

#Correlating the gap between first or last authors, and overall, with the gender gap. Smaller difference in male-biased ones
gap.data <- disc.stats %>% filter(gender.ratio.at.present < 50) %>% group_by(discipline) %>% 
  summarise(overall.gr = gender.ratio.at.present[1],
            n.authors = n.authors[1],
            diff1 = gender.ratio.at.present[2] - gender.ratio.at.present[1],
            diff2 = gender.ratio.at.present[3] - gender.ratio.at.present[1])

gap.figure <- grid.arrange(
  ggplot(gap.data, aes(x = overall.gr, y = diff1, weight = n.authors))+ geom_hline(yintercept = 0, linetype = 2) + stat_smooth(method = "gam", formula = y ~ s(x, k = 3), colour = "black") + geom_point(alpha=0.7, size = 0.8, colour = rgb(178,56,71,maxColorValue = 255)) + xlab(NULL) + ylab("Relative % women first authors") + scale_y_continuous(limits = c(-17,11)) + theme_classic(),
  ggplot(gap.data, aes(x = overall.gr, y = diff2, weight = n.authors))+ geom_hline(yintercept = 0, linetype = 2) + stat_smooth(method = "gam", formula = y ~ s(x, k = 3), colour = "black") + geom_point(alpha=0.7, size = 0.8, colour = rgb(66,124,171,maxColorValue = 255)) + xlab(NULL) + ylab("Relative % women last authors")+ scale_y_continuous(limits = c(-17,11)) + theme_classic(),
ncol=2, bottom = "% women authors in 2016")
ggsave(gap.figure, file = "../figures/gap figure.pdf", height = 3.7, width = 7)


########################################################################################################         
# WHAT TYPES OF JOURNAL HAVE THE MOST MEN/WOMEN? LOOKING AT IF, REVIEW JOURNALS, AND OA
# Make the heat map figure, and do the accompanying stats
########################################################################################################        

journal.data <- left_join(read.csv("../data for analysis/Statistics for each journal.csv", stringsAsFactors = F), journals_sqlite %>% rename(journal = short.title) %>% select(journal, IF, ISSN), copy = T)
OA.data <-  read.csv("../data/DOAJ open access data.csv", stringsAsFactors = F) %>% select(Journal.title, Journal.ISSN..print.version.,  Journal.EISSN..online.version., Journal.license)
journal.data$is.OA <- "No"      # Add in the Open Access data
journal.data$is.OA[journal.data$ISSN %in% OA.data$Journal.ISSN..print.version.] <- "Yes"
journal.data$is.OA[journal.data$ISSN %in% OA.data$Journal.EISSN..online.version.] <- "Yes"
journal.data$is.OA[grep("BMC", journal.data$journal)] <- "Yes"    # Some well-known OA journals did not get added (e.g. because publications can have several ISSNs, or database is incomplete), so do it manually
journal.data$is.OA[grep("plos ", tolower(journal.data$journal))] <- "Yes"
journal.data$is.OA[grep("PeerJ", journal.data$journal)] <- "Yes"
journal.data$is.OA[grep("F1000Res", journal.data$journal)] <- "Yes"
journal.data$Review <- "No"    # Review journals are defined as those with "Rev" or "Trends" in the title
journal.data$Review[unique(c(grep("Rev", journal.data$journal),grep("Trends", journal.data$journal)))] <- "Yes"

# Now throw out all the data where the journal discipline is listed as "Unclassified", and where we do not know the impact factor. Also restrict to "Overall" measurements
journal.data <- journal.data %>% filter(discipline != "Unclassified" & !is.na(IF) & position == "overall")

# Only include journals where we observed the gender of at least 500 authors (thus, the estimated gender ratios are accurate)
journal.data <- journal.data %>% filter(n.authors >= 500)  

# Calculate the impact factor of each journal relative to other journals in its discipline (using residuals from a random effects model)
journal.data$logIF <- log10(journal.data$IF) # Take the log10 of the impact factor
journal.data$predictedIF <- resid(lmer(logIF ~ (1|discipline), data = journal.data)) # Calculate residual log10 IF for each journal (i.e. IF relative to the field that it is in)

# Define journals as "top tier" if their absolute IF is in the top 25% for their discipline
benchmark <- melt(tapply(journal.data$IF, journal.data$discipline, function(x) quantile(x, probs = 0.75))) %>% rename(discipline = Var1, cutoff = value) # find 75% IF quantile for each discipline
benchmark <- benchmark$cutoff[match(journal.data$discipline, benchmark$discipline)] # make a column giving the appropriate benchmark for each journal
journal.data$IF.category <- "Standard"
journal.data$IF.category[journal.data$IF >= benchmark] <- "High impact"

# Calculate "adjusted % women" for each journal. A journal in a discpline where the average % women authors is 30% would have a score of +5 if it had 35% women, and -10 if it had 20% women.
benchmark <- melt(tapply(journal.data$gender.ratio.at.present, journal.data$discipline, mean)) %>% rename(discipline = Var1, cutoff = value) # find mean IF for each discipline
benchmark <- benchmark$cutoff[match(journal.data$discipline, benchmark$discipline)] # make a column giving the appropriate benchmark for each journal
journal.data$adjusted.gender.ratio <- journal.data$gender.ratio.at.present - benchmark # subtract the appropriate mean from each one to get adjusted % women

heat.map.data <- journal.data %>% group_by(IF.category, Review, is.OA) %>% summarise(mean = mean(adjusted.gender.ratio), n.authors = sum(n.authors), n.papers = sum(n.papers), n.journals = n()) %>% mutate(is.OA = replace(is.OA, is.OA=="No", "Not Open Access"), is.OA = replace(is.OA, is.OA=="Yes", "Open Access")) %>% as.data.frame()  %>% mutate(IF.category = factor(IF.category, levels = c("Standard", "High impact")))
sample.size.labels <- heat.map.data %>% add_row(IF.category="High impact", Review="Yes", is.OA= "Open Access", mean=NA, n.authors=0, n.papers=0, n.journals = 0)

heat.map <- heat.map.data %>% ggplot(aes(x = IF.category, y = Review, fill = mean)) + geom_tile(colour = "white", size = 2.7) + facet_wrap(~is.OA) + scale_fill_gradient2(low = brewer.pal(7, "Blues")[6], mid = "white", high = brewer.pal(7, "Reds")[4], name = "Relative\n% women\nauthors") + scale_x_discrete(expand = c(0,0)) + scale_y_discrete(expand = c(0,0)) + xlab("Impact factor") + ylab("Review journal") + theme(panel.grid.major = element_blank(), panel.border = element_rect(colour = "white", size =2.7, fill = NA), strip.background = element_rect(colour = NULL, size =1, fill = "white"), axis.ticks = element_blank()) + geom_text(data = sample.size.labels, aes(label = n.journals))

ggsave(heat.map, file = "../figures/Figure 4 - heat map.pdf", width = 7.7, height = 3.76) # save the heat map


# Accompanying stats for the heat map
m1 <- (lm(adjusted.gender.ratio ~ predictedIF * Review * is.OA, data = journal.data)) # Define the full model 
dredge.output  <- with(options(na.action = "na.fail"), MuMIn::dredge(m1)) # Rank all the component models by AICc
cbind(coef(MuMIn::model.avg(dredge.output)), confint(MuMIn::model.avg(dredge.output))) # Get the model averaged parameter values
MuMIn::importance(dredge.output)
summary(lm(adjusted.gender.ratio ~ predictedIF + Review * is.OA, data = journal.data)) # This is the top model, as ranked by AICc
nrow(journal.data) # There are 2721 journals that we are able to use in this analysis (i.e. where we have at least 500 authors gendered, plus IF information)
# Table for export to LaTex
xtable::xtable(car::Anova(lm(adjusted.gender.ratio ~ predictedIF + Review * is.OA, data = journal.data), type = "III"))




########################################################################################################
# CALCULATE THE GENDER RATIO, RATE OF CHANGE, AND YEARS TO PARITY FOR ALL COMBOS OF COUNTRY AND DISCIPLINE
# Only includes combinations where we have 250 authors and 100 papers overall, and at least 5 different 
# years in which 50 authors' genders could be identified
########################################################################################################

# The data we made earlier for the web app tells us which combos of country and discipline have a high enough sample size to predict the gender ratio accurately
# Reminder: we specified that the combination must have data for 100 papers and 250 authors, and data from at least 5 different years in which at least 50 authors were gendered
# So, let's load up the combinations that we need to re-visit and analyse using "run.sim = T" in order to fit a curve and bootstrap some CIs for the present gender ratio
web.app.data <- as_tibble(jsonlite::fromJSON(txt = "../Errol's web app/genderGap/data_for_web_app_unparsed.json")) 
web.app.data <- web.app.data %>% filter(Position == "Overall", Journal == "allJournals", Country != "allCountries") %>% select(Discipline, Country) 

# Now run the model for each of these combinations (n=2654), and get the data we need to make nice graphs split by country
get.gender.by.disc.and.country <- function(to.do){ # Expects a single value for discipline and country, like c("Medicine", "Italy")
  output <- "not done yet"
  attempts.left <- 5   # Very occasionally fails when the curve is almost totally flat, due to negative numbers in the likelihood or something. Just re-run when this happens
  FT <- "disc"
  if(to.do$Discipline == "allDisciplines") FT <- "everything.for.one.country"
  while(attempts.left > 0) {
    try(
      output <- gender.stats(item = to.do$Discipline, country = tolower(to.do$Country), filter.type = FT, position = "overall", data.source = pubmed_sqlite, plot = F, chunk.size = 50, nChunks = 20, export.JSON = F, run.sim = T, print.each.one = T)
    )
    if(!is.character(output)) attempts.left <- 0 # If the function worked, "output" will have been overwritten with a dataframe or NULL, so we can quit already
    attempts.left <- attempts.left - 1
  }
  if(!is.character(output)) {
    if(is.null(output)) return(output)
    else return(data.frame(country = to.do$Country, output)) # If it worked, return the dataframe or NULL
  }
  return(NULL) # If it didn't work after 5 attempts, return NULL
}
unique.countries <- tolower(unique(web.app.data$Country))
disc.and.country.gender.ratios <- do.call("rbind", lapply(1:nrow(web.app.data), function(x) get.gender.by.disc.and.country(web.app.data[x,])))
# write.csv(disc.and.country.gender.ratios, file = "../data for analysis/discipline.by.country.gender.data.csv", row.names = F)



########################################################################################################
# Do the same thing again, but for each unique country across all disciplines. Needed for the country figure.
########################################################################################################

# only bother doing countries where we measured the gender of at least 1000 authors
unique.countries <- (read.csv("../outputs/genderizing.success.rate.by.country.csv", stringsAsFactors = F) %>% filter(success.rate * n > 999, country != "United", country != "Kingdom"))[,1] 

all.discs.for.each.country <- do.call("rbind", lapply(1:length(unique.countries), function(x) {
  print(unique.countries[x])
  output <- NULL
  try(output <- gender.stats(item = "all", filter.type = "everything.for.one.country", position = "overall", country = tolower(unique.countries[x]), data.source = pubmed_sqlite, chunk.size = 20, nChunks = 50, plot = F, print.each.one = T, run.sim = T, export.JSON = F, verbose = T))
  return(output)
}))
all.discs.for.each.country <- data.frame(country = unique.countries, all.discs.for.each.country, stringsAsFactors = F)
write.csv(all.discs.for.each.country, file = "../data for analysis/all.discs.for.each.country.csv", row.names = F)



#############################################################################################################################
# USING THAT DATA, MAKE FIGURES SHOWING AUTHOR GENDER RATIOS SPLIT BY COUNTY (AND COUNTRY+DISCIPLINE) - USES PUBMED DATA ONLY
#############################################################################################################################

country.rankings <- read.csv("../data for analysis/discipline.by.country.gender.data.csv", stringsAsFactors = F)[,c(1,3,5:9)]
all.discs.for.each.country <- read.csv("../data for analysis/all.discs.for.each.country.csv", stringsAsFactors = F)[,c(1,5:9)]

# Arg: Which elements of unique.discs to plot (or just give specific discipline names to disc argument)
# Arg: How many countries to plot (picks the most productive ones)
# Arg: minimum.n excludes all the gender ratio estimates below this threshold. Must measure at least n authors' genders to get in
country.discipline.plot <- function(disc.numbers = 1:12, discs = NULL, nCountries = 50, minimum.n = 50, xlim = c(0, 100), breaks = c(0, 25, 50, 75, 100), nrow = 2, xlab = "Estimated % women authors in 2016 \u00B1 95% CIs", free.x = T){
  
  # Remove data points where the CIs are > 20% (clearly we didn't manage to estimate the gender gap very well in these cases)
  country.rankings <- country.rankings %>% filter(up.CI.1 - low.CI.1 < 20)
  
  unique.discs <- table(country.rankings$discipline)[table(country.rankings$discipline) > 9] %>% names %>% sort # These are the discs where we have data for at least 5 countries
  unique.discs <- unique.discs[unique.discs != "allDisciplines"]
  country.rankings <- country.rankings %>% rename(n = n.authors, prop.female = gender.ratio.at.present) %>% mutate(col = log10(n)) 

  if(is.null(discs)) country.plot.data <- country.rankings %>% filter(discipline %in% unique.discs[disc.numbers])
  else country.plot.data <- country.rankings %>% filter(discipline %in% discs) 
  
  # Find the overall gender ratio within each discipline (i.e. across all countries), by calculating from all the papers in the disc-country dataset
  disc.averages <- country.rankings %>% group_by(discipline) %>% summarise(mean.prop = mean((prop.female*n)/n , na.rm = T)) %>% as.data.frame
  disc.averages <- disc.averages %>% filter(as.character(discipline) %in% as.character(country.plot.data$discipline))
  
  name.mappings <- data.frame(discipline = c("Physical and Rehabilitation Medicine", "Chemistry Techniques, Analytical", "Clinical Laboratory Techniques", "Sexually Transmitted Diseases", "Substance-Related Disorders", "Speech-Language Pathology", "Laboratory Animal Science", "Health Services Research", "Complementary Therapies", "Biomedical Engineering", "Allergy and Immunology", "Antineoplastic Agents", "Occupational Medicine", "Anti-Infective Agents", "Reproductive Medicine", "Computational Biology", "Communicable Diseases", "Nutritional Sciences", "Environmental Health", "History of Medicine", "Statistics as Topic", "Behavioral Sciences", "Primary Health Care", "Medical Informatics", "Veterinary Medicine", "Emergency Medicine"),
                              short = c("Phys. Rehabil.", "Anal. Chem. Tech.", "Clin. Lab. Tech.", "STDs", "Subs. Related Disord.", "Speech-Lang. Pathol.", "Lab. Animal Sci.", "Health Serv. Res.", "Complement. Ther.", "Biomed. Eng.", "Allergy & Immunol.", "Antineoplastics", "Occupational Med.", "Anti-Infect. Agents", "Reproductive Med.", "Comput. Biology", "Commun. Disease", "Nutri. Sci.", "Enviro. Health", "History of Med.", "Stat. as Topic", "Behav. Sci.", "Prim. Health Care", "Med. Informatics", "Veterinary Med.", "Emergency Med."),
                              stringsAsFactors = F)
  country.plot.data <- suppressMessages(left_join(country.plot.data, name.mappings))
  country.plot.data$discipline[!is.na(country.plot.data$short)] <- country.plot.data$short[!is.na(country.plot.data$short)] 
  disc.averages <- suppressMessages(left_join(disc.averages, name.mappings))
  disc.averages$discipline[!is.na(disc.averages$short)] <- disc.averages$short[!is.na(disc.averages$short)] 
  
  # Countries ranked by sample size
  all.countries.n <- (country.rankings %>% group_by(country) %>% summarise(n = sum(n)) %>% arrange(-n) %>% select(country) %>% as.data.frame)[,1]
  countries <- all.countries.n[1:nCountries] # The top n countries, by sample size, are the ones that get plotted
  
  # Countries ranked by gender ratio
  all.countries.ratio <- (all.discs.for.each.country %>% arrange(gender.ratio.at.present) %>% select(country) %>% as.data.frame)[,1]
  country.levels <- all.countries.ratio[all.countries.ratio %in% countries] # Make the levels, for ordering the countries by gender ratio in the plot
  country.plot.data <- country.plot.data %>% filter(country %in% countries) %>% arrange(discipline, prop.female) %>% mutate(country = factor(country, levels = country.levels)) %>% filter(!is.na(country))

  p <- ggplot(country.plot.data, aes(y = country, x = prop.female)) + geom_vline(xintercept = 50, colour = "pink") + geom_vline(data = disc.averages, aes(xintercept = mean.prop), linetype=2, colour = "darkgrey") + geom_errorbarh(aes(xmin = low.CI.1, xmax = up.CI.1), height = 0) + geom_point(size = 1) + xlab(xlab) + ylab(NULL) + theme(strip.text = element_text(size=8)) 
  if(free.x) return(p + facet_wrap(~discipline, nrow=nrow, scales = "free_x"))
  else return(p + facet_wrap(~discipline, nrow=nrow) + scale_x_continuous(breaks = breaks) + coord_cartesian(xlim = xlim))
}

# Here is the country plot from the main paper
side.panels <- country.discipline.plot(discs = c("Dermatology", "Cardiology", "Chemistry", "Biophysics"), xlim = c(10,70), breaks = c(10,20,30,40,50,60,70), xlab = NULL, nCountries = 50, free.x = FALSE)

# Find the gender ratio for each country (Across disciplines), and for all of PubMed, for this country-disc dataset
pubmed.overall <- all.discs.for.each.country %>% summarise(mean.prop = mean((gender.ratio.at.present*n.authors)/n.authors)) %>% as.numeric
countries.overall <- left_join(all.discs.for.each.country,  read.csv("../outputs/genderizing.success.rate.by.country.csv", stringsAsFactors = F) %>% select(country, success.rate)) %>% arrange(gender.ratio.at.present) %>% mutate(country = factor(country, levels = country))

country.plot <- countries.overall %>% ggplot(aes(x = country, y = gender.ratio.at.present, fill = 100*success.rate)) + geom_hline(yintercept = pubmed.overall, linetype = 2) + geom_hline(yintercept = 50, colour = "pink") + geom_errorbar(width = 0, aes(ymin = low.CI.1, ymax = up.CI.1)) + geom_bar(stat="identity", colour = "black", size = 0.5) + scale_y_continuous(expand = c(0,0), breaks = c(0,10,20,30,40,50,60), limits = c(0,63)) + scale_fill_gradientn(colours = brewer.pal(3, "Blues"), name = "% authors\nwith known\ngender", limits = c(0,100)) + coord_flip() + ylab(NULL) + xlab(NULL) + labs(fill = "Sample\nsize") + theme(legend.position = c(0.84,0.2))
grid.arrange(country.plot, side.panels, ncol = 2, bottom = "Estimated % women authors in 2016 (\u00B1 95% CIs)")
ggsave(grid.arrange(country.plot, side.panels, ncol = 2, bottom = "Estimated % women authors in 2016 (\u00B1 95% CIs)"), file = "../figures/Figure 3 - country info.pdf", height = 11, width = 10)

# Here are all the supplmentary plots of all the discipline-country combinations where we have good estiamtes of the gender ratio of 10+ countries
ggsave(country.discipline.plot(1:12), file = "../figures/country.disciplines1.pdf", width = 8, height = 11.6)
ggsave(country.discipline.plot(13:24), file = "../figures/country.disciplines2.pdf", width = 8, height = 11.6)
ggsave(country.discipline.plot(25:36), file = "../figures/country.disciplines3.pdf", width = 8, height = 11.6)
ggsave(country.discipline.plot(37:48), file = "../figures/country.disciplines4.pdf", width = 8, height = 11.6)
ggsave(country.discipline.plot(49:60), file = "../figures/country.disciplines5.pdf", width = 8, height = 11.6)
ggsave(country.discipline.plot(72:81), file = "../figures/country.disciplines6.pdf", width = 8, height = 11.6)


##########################################################################################
# Statistical analysis searching for country-level predictors of author gender ratio
##########################################################################################

# Load the human development and gender inequality data from the UN
un.data <- left_join(read.csv("../data/UNDP HDI data.csv", stringsAsFactors = F)[,1:7], read.csv("../data/UNDP GII data.csv", stringsAsFactors = F)[,1:9])

# Rename the variables to R-friendly format
un.data <- with(un.data, data.frame(country = Country,
                                    HDI = Human.Development.Index..HDI.,
                                    life.expectancy = Life.expectancy.at.birth,
                                    years.of.school = Mean.years.of.schooling, 
                                    GNI = Gross.national.income..GNI..per.capita,
                                    GII = Gender.Inequality.Index, 
                                    maternal.mortality = Maternal.mortality.ratio..deaths.per.100.000.live.births.,
                                    adolescent.births = Adolescent.birth.rate..births.per.1.000.women.ages.15.19.,
                                    seats.in.parliament = Share.of.seats.in.parliament....held.by.women.,
                                    education.gap = Men.with.at.least.some.secondary.education....ages.25.and.older. - Women.with.at.least.some.secondary.education....ages.25.and.older.,
                                    labour.gap = Men.labour.force.participation.rate....ages.15.and.older. - Women.labour.force.participation.rate....ages.15.and.older.,
                                    stringsAsFactors = F
)) %>% filter(country != "")

# Scale the UN predictors so that all have mean 0, variance 1 (makes it easier to compare their relative effect sizes on gender ratio)
un.data <- data.frame(country = un.data$country, apply(un.data[,2:ncol(un.data)], 2, function(x) str_replace(x,",","") %>% # The UN hard-coded the commas into their numbers (like 20,000), so remove them
                                                         as.numeric %>%
                                                         scale), stringsAsFactors = F)

# Define countries where women are MORE abundant in education and labour as having no gap in these areas (i.e. there can be no 'negative gap')
un.data$education.gap[un.data$education.gap < 0] <- 0
un.data$labour.gap[un.data$labour.gap < 0] <- 0

# Rename some countries in the UN data for proper matching with our own data
un.data$country[un.data$country == "United States"] <- "USA"
un.data$country[un.data$country == "Russian Federation"] <- "Russia"
un.data$country[un.data$country == "Viet Nam"] <- "Vietnam"
un.data$country[un.data$country == "Brunei Darussalam"] <- "Brunei"
un.data$country[un.data$country == "Côte d'Ivoir"] <- "Ivory Coast"

merged <- left_join(read.csv("../data for analysis/Author counts for journal country and position.csv", stringsAsFactors = F) %>% 
                      spread(key = gender, value = n) %>% 
                      rename(nFemales = F, nMales = M, nUnknown = U) %>%
                      mutate(year = as.numeric(scale(year))), #scale year to have the same mean and variance as the UN predictors (0 and 1)
                    un.data)
merged$nFemales[is.na(merged$nFemales)] <- 0
merged$nMales[is.na(merged$nMales)] <- 0
merged$nUnknown[is.na(merged$nUnknown)] <- 0
merged <- merged %>% filter(country != "Taiwan" & country != "Unknown") # Taiwan is not in the UN data separately 
# These places do not have complete UN data and will be removed
# "Hong Kong"      "Nigeria"        "Puerto Rico"    "Madagascar"     "Brunei"         "Uzbekistan"     "Guadeloupe"     "Ivory Coast"    "Reunion Island" "Bermuda"        "Saint Lucia"
merged <- merged[complete.cases(merged), ] # Remove countries with incomplete UN data
merged <- left_join(merged, journals_sqlite %>% select(short.title, Discipline) %>% rename(journal = short.title, discipline = Discipline), copy = TRUE) # Add the disciplines
merged <- merged %>% filter(nFemales + nMales >= 10) # throw out combinations with fewer than 10 authors measured
merged$gender.ratio <- with(merged, nFemales / (nFemales + nMales)) # Compute the gender ratio, in formats that work with lmer and glmer

# Sample size for the following analysis (number of authors)
sum(c(merged$nFemales, merged$nMales))

# First try to run a binomial GLMM. 
# It never converges properly, no matter the optimiser used. The data are just not amenable to GLMM I guess.
library(afex) # This package has the 'allFit' function, which runs lots of optimisers on the same model to look for one that works
library(optimx) # This packages has extra optimiser algorithms
resp <- with(merged, cbind(nFemales, nMales)) # Make a two-column response vector of number of men/women for each combination
models1 <- allFit(glmer(resp ~  year + position + life.expectancy + years.of.school + GNI + adolescent.births + seats.in.parliament + education.gap + labour.gap + (1|journal) + (1|discipline) + (1|country), data = merged, family = "binomial"))  # Lots of failures to coverge, irrespective of optimiser used


# Fit the data as a linear model instead of a GLMM (with % women as the response variable)
full.model1 <- lmerTest::lmer(100*gender.ratio ~ year + position + life.expectancy + years.of.school + GNI + adolescent.births + seats.in.parliament + education.gap + labour.gap + (year|journal) + (year|discipline) + (year|country), data = merged)
summary.full <- summary(full.model1)
anova.full <- car::Anova(full.model1, type = "III")
MuMIn::r.squaredGLMM(full.model1) # r squared information

# Get the random and fixed effects tables in LaTex format
print(xtable::xtable((VarCorr(full.model1) %>% as.data.frame() %>% select(grp,var1,vcov))[c(1,2,4,5,7,8,10),] %>% mutate(prop = 100*vcov/sum(vcov)), digits = c(1,1,1,2,1)), include.rownames = F)
xtable::xtable(summary.full$coefficients, digits = c(1,2,2,1,2,4))
xtable::xtable(anova.full)



# Make a supplementary figure of the fixed effects coefficients
correlates.plot <- summary.full$coefficients %>% as.data.frame() %>% mutate(low = Estimate - 1.96 * `Std. Error`, up = Estimate + 1.96 * `Std. Error`, Parameter = rownames(summary.full$coefficients)) %>% filter(Parameter != "(Intercept)") %>% 
  mutate(Parameter = c("Publication date", "Being last author", "Being middle author", "Being a single author", "Life expectancy", "Mean years of school \n(both genders)", "Gross national income \n(per capita)", "Adolescent birth rate", "% women in parliament", "Gender gap in education", "Gender gap in labour force")) %>%
  arrange(Estimate) %>% mutate(Parameter = factor(Parameter, levels =Parameter)) %>% ggplot(aes(x = Parameter, y = Estimate)) + geom_hline(yintercept = 0, linetype = 2) + geom_errorbar(aes(ymin =low, ymax = up), width = 0) + geom_point(size = 2) + coord_flip() + xlab(NULL) + ylab("Effect on % women authors ")
ggsave(correlates.plot, width = 5, height = 4, file = "../figures/correlates.plot.pdf")


################################################################################################################
# Does my authorship data correlate with the NSF's data on the % women working in different fields?
# This serves as a sanity check on my data, and gives some insight into the question: "Do fields with fewer
# female authors have fewer female researchers, or does this represent lower publication rates by women (or both)?"
################################################################################################################

# First get all the gender ratios for the disciplines for USA only - for all 4 author positions
combos <- expand.grid(position = c("first", "last", "single"), discipline = (read.csv("../data for analysis/discipline.by.country.gender.data.csv") %>% filter(country == "USA") %>% select(discipline))[,1] %>% as.character()) %>% mutate(position = as.character(position), discipline = as.character(discipline))
usa.gender.ratio <- do.call("rbind", lapply(1:nrow(combos), function(x) {
  output <- NULL
  try(output <- gender.stats(item = combos$discipline[x], filter.type = "disc", position = combos$position[x], country = "usa", data.source = pubmed_sqlite, chunk.size = 20, nChunks = 50, plot = F, print.each.one = T, run.sim = T, export.JSON = F, verbose = T))
  return(output)
}))
# write.csv(usa.gender.ratio, file = "../data for analysis/usa.gender.ratio.csv", row.names = F)

usa.gender.ratio <- rbind(read.csv("../data for analysis/usa.gender.ratio.csv", stringsAsFactors = F)[,c(2,3:8)],
                          read.csv("../data for analysis/discipline.by.country.gender.data.csv", stringsAsFactors = F)[,c(1,3,5:9)] %>% filter(country == "USA" & discipline != "allDisciplines") %>% mutate(position = "overall") %>% select(discipline, position, n.authors, n.papers, gender.ratio.at.present, low.CI.1,  up.CI.1))
usa.gender.ratio$position[usa.gender.ratio$position == "first"] <- "First"
usa.gender.ratio$position[usa.gender.ratio$position == "last"] <- "Last"
usa.gender.ratio$position[usa.gender.ratio$position == "single"] <- "Single"
usa.gender.ratio$position[usa.gender.ratio$position == "overall"] <- "Overall"
usa.gender.ratio$position <- factor(usa.gender.ratio$position, levels = c("First", "Last", "Single", "Overall"))

# Load data from the NSF which give the % female PhD students and postdocs in various specialties in 2014 in the USA
# There are found here: https://ncsesdata.nsf.gov/gradpostdoc/2014/index.html -> Tables 21, 22 and 38
# Merge them with the PubMed data set, using a manually-assigned mapping of NSF categories to PubMed categories (I did this blind to the data)
# Some NSF categories (e.g. Ocean Science) do not have obvious analogs on PubMed, so I just left them unmatched
usa.gender.ratio <- left_join(usa.gender.ratio, read.csv("../data/NSF data - female phd students and postdocs by subfield in 2014.csv", stringsAsFactors = F) %>% mutate(`PhD students` = Female.phd.students/(Female.phd.students + Male.phd.students), Postdocs = Female.postdocs/(Female.postdocs + Male.postdocs)) %>% select(PubMed.discipline, `PhD students`, Postdocs) %>% rename(discipline = PubMed.discipline)) %>% filter(!is.na(`PhD students`) | !is.na(Postdocs)) %>% distinct(paste(discipline, position), .keep_all = T)
usa.gender.ratio$`PhD students`[is.nan(usa.gender.ratio$`PhD students`)] <- NA
usa.gender.ratio<-usa.gender.ratio[,1:9]
usa.gender.ratio <- gather(usa.gender.ratio, key = level, value = nsf.ratio, `PhD students`,  Postdocs)  %>% arrange(discipline, position, level)


# Add the R2 and slope statistics (m = slope)
make.annotation <- function(level, pos, statistic){
  dat <-  usa.gender.ratio %>% filter(position == pos)
  if(grepl("PhD", level)) dat <- dat[grep("PhD", dat$level), ] %>% filter(!is.na(nsf.ratio) & !is.nan(nsf.ratio))
  else dat <- dat[grep("Post", dat$level), ]
  model <- lm(gender.ratio.at.present ~ nsf.ratio, data = dat)
  slope <- round(summary(model)$coefficients[2,1], 2)
  slope.cis <- round(confint(model)[2,], 2)
  if(statistic == "slope") return(paste("italic(m) == ", slope, "~ (", slope.cis[1], "-", slope.cis[2], ")", sep = ""))
  paste("R^2 ==", round(summary(model)$adj.r.squared, 2))
}
annotations <- expand.grid(level = unique(usa.gender.ratio$level), position = unique(usa.gender.ratio$position), slope = NA, r2 = NA, x = 0, y = 70) %>% 
  mutate(level = as.character(level), author.position = as.character(position))
for(i in 1:nrow(annotations)) {
  annotations$slope[i] <- with(annotations, make.annotation(level[i], position[i], "slope"))
  annotations$r2[i] <- with(annotations, make.annotation(level[i], position[i], "r2"))
}
annotations$position <- factor(annotations$position, levels = c("First", "Last", "Single", "Overall"))

# Plot the data. Very clear positive slope, but it's lower than 1:1 - this mean that the % young female scientists in a field is an overestimate of the % female authors
# This suggests A) women publish less often than men, and/or B) the gender ratio at levels above PhD/postdoc is more male-biased than at these levels (as I strong suspect it is, but the NSF does not have detailed data on science faculty in these detailed sub-disciplines, as far as I can tell - they just have data on faculty gender ratio for very broad categories like "Biology and Medicine")
nsf.figure <- usa.gender.ratio %>% ggplot(aes(x = 100*nsf.ratio, y = gender.ratio.at.present)) + stat_smooth(method = "lm") + geom_point(alpha = 0.9, size = 0.9) + geom_abline(slope = 1, intercept = 0, linetype = 2) + facet_grid(position~level) + xlab("% women early career researchers in the NSF survey") + ylab("% women USA-based authors on PubMed") + scale_y_continuous(limits = c(0,75), breaks = c(0,25,50,75)) + scale_x_continuous(limits = c(0,100), breaks = c(0,25,50,75,100)) 
nsf.figure <- nsf.figure + geom_text(data = annotations, aes(x=x, y=y, label = slope), inherit.aes=FALSE, hjust = 0, parse=T) + geom_text(data = annotations, aes(x=x, y=y-6, label = r2), inherit.aes=FALSE, hjust = 0, parse=T) 

ggsave(nsf.figure, file = "../figures/Comparison with NSF.pdf", width = 7, height = 10) 


##########################################################################################################################       
# Export some data giving the gender ratio, rate of change, and years-to-parity for all journals and disciplines
# These are uploaded as .csv files in the supplementary material, as they could be useful for reference and future studies
##########################################################################################################################  

pubmed.data.for.export <- rbind(
  read.csv("../data for analysis/Statistics for each discipline.csv", stringsAsFactors = F),
  read.csv("../data for analysis/Statistics for each journal.csv", stringsAsFactors = F)
) %>% neaten.the.data

country.specific.pubmed.data.for.export <- read.csv("../data for analysis/discipline.by.country.gender.data.csv", stringsAsFactors = F) %>% neaten.the.data

arxiv.data.for.export <- rbind(
  read.csv("../data for analysis/Statistics for ArXiv major categories.csv", stringsAsFactors = F),
  read.csv("../data for analysis/Statistics for ArXiv sub-categories.csv", stringsAsFactors = F)
) %>% neaten.the.data

journal.data.for.export <- read.csv("../data for analysis/Statistics for each journal.csv") %>% neaten.the.data

write.csv(pubmed.data.for.export, file = "../data for analysis/Supplementary dataset 1 - PubMed gender information.csv", row.names = F)
write.csv(country.specific.pubmed.data.for.export, file = "../data for analysis/Supplementary dataset 2 - PubMed gender information - split by country.csv", row.names = F)
write.csv(arxiv.data.for.export, file = "../data for analysis/Supplementary dataset 3 - ArXiv gender information.csv", row.names = F)
write.csv(journals_sqlite %>% as.data.frame() %>% select(short.title, Discipline, title) %>% rename(journal = short.title, discipline = Discipline, fullTitle = title), file = "../genderGapCode/journal_disciplines.csv", row.names = F)


################################################################################################################
# Do papers with a senior female (male) author have more females or males on them than expected? (YES, for sure)
# I decided to do this separately for individual journals, because the 'disciplines' tend to pool across 
# sub-specialties which may have different gender ratios, which screws up the test. 
# See note by Bergstrom  et al here: http://eigenfactor.org/gender/assortativity/measuring_homophily.pdf
################################################################################################################

# This function calculates Bergstrom et al's "coefficient of homophily", for each journal in 'journals'
# I chose to restrict the analysis to the last year (with 'nYears' argument), because if the gender ratio changes over time, it will generate spurious homophily
# minimum.n refers to the number of papers - we only use papers where all authors' genders are known for this analysis
gender.grouping.test.by.journal <- function(journals, minimum.n, nSimulations, nYears = 1, db = pubmed_sqlite){
  # Internal function to estimate alpha for a sample of author lists, using bootstrap sampling
  # alpha is defined as p - q, and is termed "the coefficient of homophily" (comments section, http://www.michaeleisen.org/blog/?p=1931)
  # p is the probability that a randomly-chosen co-author of a male author is a man
  # p is the probability that a randomly-chosen co-author of a female author is a man
  # Thus, positive alpha means that there are more male-male and female-female associations than expected by change. Negative alpha means the fewer than expected by chance.
  # One could probably do this using ANOVA or something, but this method is probably more robust as it's non-parametric. 
  # Bergstrom et al's note here (http://eigenfactor.org/gender/assortativity/measuring_homophily.pdf). It seems to assume that all papers have 2 authors, which I guess makes things easier to do analytically
  find.alpha <- function(counts, nBoots = 100, return.cis = FALSE){
    
    # Add vector giving the probability of getting a male author when picking the first random author
    p.first.pick.is.male <- counts$nMales / (counts$nMales + counts$nFemales)
    
    # Add vector giving the probability of getting a male co-author, after picking a MALE the first time
    p.second.pick.is.male.M <- (counts$nMales - 1) / (counts$nMales + counts$nFemales - 1)
    p.second.pick.is.male.M[p.first.pick.is.male == 0] <- 0 # If it was zero the first time, make sure it's zero now
    
    # Add vector giving the probability of getting a male co-author, after picking a FEMALE the first time
    p.second.pick.is.male.F <- (counts$nMales) / (counts$nMales + counts$nFemales - 1)
    
    num.papers <- nrow(counts)
    
    # Flip some coins to determine whether we pick a male or a female as the focal author, for each paper. 
    # Make columns to hold the co-author we'll pick next, and to identify which of the 'nBoots' replicates this is for
    results <- data.frame(first.pick.is.male = rbinom(num.papers * nBoots, 1, rep(p.first.pick.is.male, nBoots)), 
                          second.pick.is.male = NA, 
                          boot.replicate = rep(1:nBoots, each = num.papers))
    
    # Now that we have picked the first random authors, we know which and how many of these probabilities we will need
    p.second.pick.is.male.M <- rep(p.second.pick.is.male.M, nBoots)
    p.second.pick.is.male.M <- p.second.pick.is.male.M[results$first.pick.is.male == 1]   # CHECK THIS WEIRD STUFF
    
    # Same again
    p.second.pick.is.male.F <- rep(p.second.pick.is.male.F, nBoots)
    p.second.pick.is.male.F <- p.second.pick.is.male.F[results$first.pick.is.male == 0]
    
    results$second.pick.is.male[results$first.pick.is.male == 1] <- rbinom(sum(results$first.pick.is.male == 1), 1, p.second.pick.is.male.M)  # CHECK TOO
    results$second.pick.is.male[results$first.pick.is.male == 0] <- rbinom(sum(results$first.pick.is.male == 0), 1, p.second.pick.is.male.F)
    
    p.and.q.estimates <- results %>% group_by(boot.replicate, first.pick.is.male) %>% summarise(answer = sum(second.pick.is.male) / length(second.pick.is.male))
    alpha.estimates <- p.and.q.estimates$answer[p.and.q.estimates$first.pick.is.male == 1] - p.and.q.estimates$answer[p.and.q.estimates$first.pick.is.male == 0] # WARNINGS HERE
    
    if(return.cis) return(c(median(alpha.estimates), as.numeric(quantile(alpha.estimates, probs = c(0.025, 0.975)))))
    return(median(alpha.estimates)) # The median of the bootstrap samples is our best guess of alpha, given the data. Usually the median is close to the mean, but median should be more robust
  }
  
  output <- data.frame(journal = journals, n.useable.papers = NA, n.authors = NA, gender.ratio.of.sample = NA, observed.alpha = NA, low.CI.alpha = NA, up.CI.alpha = NA, two.tail.p.value = NA, stringsAsFactors = F)
  
  # Now, loop over all the journals...
  for(i in 1:length(journals)){
    focal <- db %>% filter(journal == journals[i]) %>% select(gender.95, date) %>% collect(n = Inf) %>% filter(nchar(gender.95) > 1) # only get papers with >1 author
    focal$date <- convert.dates(focal$date, type = "before2016")
    focal <- focal %>% filter(date > nYears * -1) %>% select(-date) # Restrict the analysis to "nYears" before the present date. My logic is that if the gender ratio changes over time, it will give the false impression of homophily, as explained in Bergstrom's comment here: http://www.michaeleisen.org/blog/?p=1931
    if(nrow(focal) > minimum.n){ # Proceed if the sample size is big enough
      focal <- (focal %>% as.data.frame())[,1] # Get the gender column as a character vector
      focal <- focal[!str_detect(focal, "U")] # to keep it simple, let's only include papers where we know all the authors' genders
      
      if(length(focal) > minimum.n){ # Proceed with test if the sample size is still big enough
        
        # First, make a summary of the number of M and F authors on each paper
        counts <- data.frame(nMales = str_count(focal, "M"), nFemales = str_count(focal, "F"))
        total.males <- sum(counts$nMales)
        total.females <- sum(counts$nFemales)
        
        num.authors <- sum(counts)
        num.papers <- nrow(counts)
        
        pasted.sex <- c(rep(1, total.males), rep(0, total.females)) # Males encoded by 1s here, females by 0s
        simulated.data <- data.frame(sex = pasted.sex[unlist(unname(tapply(runif(num.authors * nSimulations), rep(1:nSimulations, each = num.authors), rank)))],
                                     paperID = rep(unlist(mapply(rep, 1:num.papers, each = rowSums(counts))), nSimulations),
                                     replicate = rep(1:nSimulations, each = num.authors))
        
        simulated.data <- simulated.data %>% group_by(replicate, paperID) %>% summarise(nFemales = sum(sex == 0), nMales = sum(sex == 1)) %>% select(-paperID)
        
        simulated.alpha.under.null <- sapply(1:nSimulations, function(i) find.alpha(simulated.data[simulated.data$replicate == i,2:3]))
        observed.alpha <- find.alpha(counts, return.cis = TRUE)
        
        # The 2-tail p value is the proportion of null-simulated alpha values with a value at least as far from zero as the observed alpha
        # Thus, a small p-value means that the observed alpha is larger or smaller than expected
        two.tail.p.value <- sum(abs(simulated.alpha.under.null) >= abs(observed.alpha[1])) / length(simulated.alpha.under.null)
        
        out <- data.frame(n.useable.papers = num.papers,
                          n.authors = num.authors,
                          gender.ratio.of.sample = round(100 * total.females/(total.females + total.males), 2),
                          observed.alpha = observed.alpha[1],
                          low.CI.alpha = observed.alpha[2],
                          up.CI.alpha = observed.alpha[3],
                          two.tail.p.value = two.tail.p.value)
        output[i, 2:ncol(output)] <- out
        print(output[i, ])
      }
    }
  }
  return(output)
}

unique.journals <- (pubmed_sqlite %>% select(journal) %>% distinct() %>% collect(n = Inf) %>% as.data.frame())[,1]
gender.homophily.results <- gender.grouping.test.by.journal(journals = unique.journals, minimum.n = 50, nSimulations = 500, nYears = 1, db = pubmed_sqlite)
gender.homophily.results <- left_join(gender.homophily.results, journals_sqlite %>% rename(journal = short.title, discipline = Discipline, country = journal.country) %>% select(journal, discipline, country) %>% collect(n=Inf) %>% as.data.frame %>% distinct(journal, .keep_all = T)) # add the disciplines
gender.homophily.results <- gender.homophily.results %>% select(discipline, journal, country, n.useable.papers, n.authors, gender.ratio.of.sample, observed.alpha, low.CI.alpha, up.CI.alpha, two.tail.p.value) %>% arrange(-n.useable.papers)
gender.homophily.results <- gender.homophily.results %>% mutate(adjusted.p.value = p.adjust(two.tail.p.value, method = "BH")) # Add another column of p-values after controlling the False Discovery Rate using Benjamini-Hochberg method

# There is plenty of homophily, and no heterophily at all. The disciplines that contain more diverse topics (e.g. "Multidisciplinary") show most homophily, as expected. 
# Still, there is homophily elsewhere too
disc.summary <- gender.homophily.results %>% filter(!is.na(observed.alpha) & !is.na(discipline)) %>% 
  group_by(discipline) %>% 
  summarise(n.journals = n(), mean.alpha = mean(observed.alpha), se = sd(observed.alpha)/sqrt(length(observed.alpha)), 
            median = median(observed.alpha), 
            number.significant.homophily = sum(observed.alpha > 0 & two.tail.p.value < 0.05), 
            number.significant.heterophily = sum(observed.alpha < 0 & two.tail.p.value < 0.05),
            number.very.significant.homophily = sum(observed.alpha > 0 & adjusted.p.value < 0.05), 
            number.very.significant.heterophily = sum(observed.alpha < 0 & adjusted.p.value < 0.05)) %>% 
  arrange(-number.very.significant.homophily / n.journals) %>% as.data.frame()

disc.summary %>% filter(n.journals > 9) %>% arrange(-mean.alpha)

# Slightly MORE homophilly in journals with a more even gender ratio. Positive assortment is not strongest in fields where women are rare, as one might guess
gender.homophily.results %>% filter(!is.na(adjusted.p.value)) %>% mutate(sig = adjusted.p.value < 0.05) %>% ggplot(aes(x = gender.ratio.of.sample, y = observed.alpha)) + stat_smooth(method = "gam", formula = y ~ s(x, k=3)) + geom_point(alpha=0.3, aes(colour = sig, size = log10(n.useable.papers))) + xlab("% female authors in the journal") + ylab("Coefficient of homophilly (a)")

# Grouping by the country of the journal's publisher (e.g. Saudi journals have more gender segregation)
gender.homophily.results %>% filter(!is.na(adjusted.p.value)) %>% ggplot(aes(x = country, y = observed.alpha)) + geom_boxplot() + geom_jitter(alpha=0.2) + coord_flip()
