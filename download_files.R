#---- Script to download ------------------------------------------------------
#---- Ethan Roubenoff 
#---- Downloads CDC files and organizes them
library(rvest)
library(stringr)
library(tidyverse)
library(magrittr)

setwd("~/kriging_PAA")
if (!dir.exists("./data")) {dir.create("./data")}


url <- "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/Datasets/NVSS/USALEEP/CSV/"
files.list <- url %>% read_html() %>% html_nodes("a")
files.list %<>% str_match("\"(.*?)\"")
files.list %<>% as.data.frame() %>% select(V2) %>% mutate(V2 = as.character(V2))

files.list <- files.list[-1, ]

for (i in 1:104) {
  s <- files.list[i]
  s <- strsplit(s, "/")[[1]][9]
  download.file(
    paste0(url, s),
    paste0("./data/", s)
  )
}
