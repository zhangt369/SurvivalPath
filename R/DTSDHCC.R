#'@title  SYSUCC Hepatocellular Carcinoma Dynamic Time Series Data
#'@name DTSDHCC
#'@description Time series dataset of 2360 patients with intermediate stage
#' hepatocellular carcinoma (HCC), with each time point observation included data of 12
#' clinical variables and 2 survival outcome variables.
#'
#'Format
#'
#'A dataframe with 11684 observations and 14 variables. The data of each patient
#' at each time point is sorted into a separate row. The variable "ID" refers to
#' individual patient identification numbers. The variable "Date" refers to the time point of each observation.
#' A total of 12 clinical variables were arranged sequentially, including 1 demographic variable ("Age"),
#' 8 observational variables ("Amount of Hepatic Lesions", "Largest Diameter of Hepatic Lesions",
#' "New Lesion", "Vascular Invasion", "Local Lymph Node Metastasis", "Distant Metastasis",
#' "Child_pugh_score" and "AFP"), 3 treatment variables ("TargetedTherapy", "Embolization", "Resection")
#' and 2 outcome variables ("Status_of_death","OStime_day"). The missing values of the original dataset have been
#' filled using Random forest regression method.

#'@author Lujun Shen,Tao Zhang
NULL
