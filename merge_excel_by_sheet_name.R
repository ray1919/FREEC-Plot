#!/usr/bin/env Rscript
# Author: Ryan Zhao
# Date: 2021-01-28
# Purpose: merge excel by sheet name

library(getopt)
library(openxlsx)
suppressMessages(library(dplyr))
COLOR_LST = c("#264653","#2A9D8F","#E9C46A","#F4A261","#E76F51", "#a94442")

opt_spec <- matrix(c('help', 'h', 0, 'logical', 'help manual',
                     'dir', 'd', 1, 'character', 'dir path for .xlsx',
                     'pattern', 'p', 1, 'character', 'file pattern for .xlsx in target dir'),
                   byrow=TRUE, ncol=5)
opt = getopt(opt_spec, commandArgs(TRUE))

if (!is.null(opt$help)){
  cat(getopt(opt_spec, usage=TRUE))
  q(save='no', status=1)
}
stopifnot(dir.exists(opt$dir))
# opt <- list()
# opt$pattern = "chr8"
# opt$dir <- "chr8_plots"

filelist <- dir(path = opt$dir, pattern = paste0(opt$pattern, ".*xlsx$"))

sheetnames <- tibble()
for (filename in filelist) {
  filepath <- paste0(opt$dir, "/", filename)
  sheetnames <- bind_rows(sheetnames,
                          tibble(filepath,
                                 sample = sub("_.*", "", filename),
                                 sheetname = getSheetNames(filepath)))
}
sheetnames <- sheetnames %>% mutate(sheet_new = sub(".* - ", "", sheetname))

wb <- createWorkbook(creator = "Ryan")
sheetnames_uniq <- unique(sheetnames$sheet_new)
for (i in 1:length(sheetnames_uniq)) {
  cat(paste("Add worksheet:", sheetnames_uniq[i], "\n"))
  addWorksheet(wb, sheetName = sheetnames_uniq[i], tabColour = COLOR_LST[i])
  dt <- tibble()
  filepath_sheet <- sheetnames %>% filter(sheet_new == sheetnames_uniq[i])
  for (j in 1:nrow(filepath_sheet)) {
    cat(paste("read Excel:", filepath_sheet$filepath[j], "\n"))
    excel <- read.xlsx(filepath_sheet$filepath[j], sheet = filepath_sheet$sheetname[j])
    dt <- bind_rows(dt, excel %>% mutate(Sample = filepath_sheet$sample[j]) %>%
                      relocate(Sample))
  }
  writeData(wb, sheet = sheetnames_uniq[i], x = dt, withFilter = F)
}

saveWorkbook(wb, paste0(opt$pattern, ".xlsx"), overwrite = T)
