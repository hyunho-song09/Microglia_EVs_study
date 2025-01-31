library(xlsx)

load_data <- function(filename, sheet) {
  filepath <- paste0(data_path, filename)
  read.xlsx(filepath, sheetName = sheet)
}

# Load EVs datasets
input.BC.df <- load_data("example_data_01.xlsx", "input_fig_BC")
input.D.df <- load_data("example_data_01.xlsx", "input_fig_D")
input.D.ex <- input.D.df[c(1:3,7:9),]
input.D.mi <- input.D.df[c(4:6,10:12),]

# Load Metabolomic datasets
cell.df <- load_data("example_data_02.xlsx",  "ex2")
media.df <- load_data("example_data_03.xlsx", "ex3")
EV.df <- load_data("example_data_04.xlsx", "ex4")