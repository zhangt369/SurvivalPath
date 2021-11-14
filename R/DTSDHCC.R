#'@title  SYSUCC Hepatocellular Carcinoma Dynamic Time Series Data
#' @docType Dynamic time series data of hepatocellular carcinoman (DTSDHCC)
#' @name DTSDHCC
#' @description Time series dataset of 2365 patients with intermediate stage
#' hepatocellular carcinoma (HCC), with each time point included data of 12
#' clinical variables and 3 survival outcome variables.
#'
#'Format
#'
#'A dataframe with 11707 observations and 15 variables. The data of each patient
#' at each time point is sorted into a separate row. The variable "ID" refers to
#' individual patient. The variable "Date" refers to the time point of each observation.
#' A total of 12 clinical variables were arranged sequentially, including 1 demographic variable ("Age"),
#' 8 observational variables ("Amount of Hepatic Lesions", "Largest Diameter of Hepatic Lesions (mm)",
#' "New Lesion", "Vascular Invasion", "Local Lymph Node Metastasis", "Distant Metastasis",
#' "Child_pugh_score" and "AFP"), 3 treatment variables ("TargetedTherapy", "Embolization", "Resection")
#' and 3 outcome variables ("Status_of_death","OStime_day", "Last_follow_up").

#'@author Lujun Shen, Jinqing Mo, Tao Zhang
NULL
