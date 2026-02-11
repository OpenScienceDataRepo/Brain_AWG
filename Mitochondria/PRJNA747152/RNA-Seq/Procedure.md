# Procedure

## QC and DESeq2

The raw counts were not available. When I reviewed the original paper, however, I found that Supplementary Table 1-5 contained DEGs for the following conditions:
- 0.4Gy LDR + 10Gy irradiated vs. unirradiated larval male brains
- 0.4Gy LDR + 10Gy vs. 10Gy irradiated larval male brains
- 0.4Gy LDR + 10Gy vs. 0.4Gy LDR irradiated larval male brains
- 10Gy irradiated vs. unirradiated larval male brains
- 10Gy vs. 0.4Gy LDR irradiated larval male brains

This meant that the QC and DESeq2 was already completed by the original researchers, making the data ready for the rest of the workflow.

*Porrazzo A, Cipressa F, De Gregorio A, De Pittà C, Sales G, Ciapponi L, Morciano P, Esposito G, Tabocchini MA, Cenci G. Low dose rate γ-irradiation protects fruit fly chromosomes from double strand breaks and telomere fusions by reducing the esi-RNA biogenesis factor Loquacious. Commun Biol. 2022 Sep 3;5(1):905. doi: 10.1038/s42003-022-03885-w. Erratum in: Commun Biol. 2022 Sep 29;5(1):1033. doi: 10.1038/s42003-022-03984-8. PMID: 36057690; PMCID: PMC9440893.*

The tables were extracted from the document using the following code:
```
setwd("/Volumes/Marians SSD/ADBR Mito/PRJNA747152")

# 1. Install + Load Packages

if (!requireNamespace("pdftools", quietly = TRUE)) install.packages("pdftools")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(pdftools)
library(stringr)
library(dplyr)

# 2. Load PDF

file <- "Supplementary_Data.pdf"
pdf_text_data <- pdf_text(file)

cat("Total pages in PDF:", length(pdf_text_data), "\n")

# 3. Find Pages with Supplementary Tables 1–5

table_pages <- c()

for (i in seq_along(pdf_text_data)) {
  if (str_detect(pdf_text_data[i], "Supplementary Table 1") |
      str_detect(pdf_text_data[i], "Supplementary Table 2") |
      str_detect(pdf_text_data[i], "Supplementary Table 3") |
      str_detect(pdf_text_data[i], "Supplementary Table 4") |
      str_detect(pdf_text_data[i], "Supplementary Table 5")) {
    
    cat("Found table on page:", i, "\n")
    table_pages <- c(table_pages, i)
  }
}

# 4. Extract Tables 1–5

tables_list <- list()

for (i in seq_along(table_pages)) {
  
  raw_text <- pdf_text_data[table_pages[i]]
  
  # Split into lines
  lines <- unlist(str_split(raw_text, "\n"))
  
  # Remove empty lines
  lines <- lines[lines != ""]
  
  # Remove header lines mentioning Supplementary
  lines <- lines[!str_detect(lines, "Supplementary")]
  
  # Convert to dataframe (split on multiple spaces)
  df <- read.table(
    text = lines,
    header = TRUE,
    sep = "",
    fill = TRUE,
    stringsAsFactors = FALSE
  )
  
  rownames(df) <- NULL
  
  tables_list[[i]] <- df
  
  cat("Extracted Table", i, "\n")
}

# 5. Export Tables Individually

for (i in seq_along(tables_list)) {
  write.csv(
    tables_list[[i]],
    file = paste0("Supp_Table_", i, ".csv"),
    row.names = FALSE
  )
}

cat("✅ Tables 1–5 exported as CSV files in working directory.\n")
```

## Visualization
```
