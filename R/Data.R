#' Gene expression profiling of chronic lung disease for the Lung Genomics Research Consortium (GSE47460)
#' Microarray data from total RNA extracted from whole lung homogenate from subjects undergoing thoracic surgery.
#' @format a list of four components.
#' \itemize{
#' \item{X: }{Expression matrix with 319 samples(Rows) and 2000 genes(Cols) after preprocessing and centering}
#' \item{Y: }{7 outcomes (Cols) of 319 samples (Rows): "fev1pd1a", "fvcprd1", "bode", "ratiopre", "RV", "WBCDIFF1", "WBCDIFF4"}
#' \item{Z: }{intercept and covariates of age, gender and BMI for 319 samples}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47460}
"LungData"

#' Gene expression profiling of METARBIC breast cancer dataset, 271 triple negative breast cancer samples are picked.
#' @format a list of five components.
#' \itemize{
#' \item{X: }{Expression matrix with 271 samples(Cols) and 2000  genes(Rows) after preprocessing and centering}
#' \item{Y: }{3 outcomes (Cols) of 271 samples (Rows): "Tumor_size", "lymph_nodes", "OS"(survival) }
#' \item{Y.ind: }{censoring indicators for survival outcome where 1 means events and 0 means censoring}
#' \item{Z: }{intercept and covariates of age and Chemotherapy for 271 samples}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47460}
"TNBC_Metabric"

#' Gene expression profiling of SCAN-B breast cancer dataset, 141 triple negative breast cancer samples are picked.
#' @format a list of five components.
#' \itemize{
#' \item{X: }{Expression matrix with 141 samples(Cols) and 2000  genes(Rows) after preprocessing and centering}
#' \item{Y: }{3 outcomes (Cols) of 141 samples (Rows): "Tumor_size", "lymph_nodes", "OS"(survival) }
#' \item{Y.ind: }{censoring indicators for survival outcome where 1 means events and 0 means censoring}
#' \item{Z: }{intercept and covariates of age and Chemotherapy for 141 samples}
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60789}
"TNBC_ScanB"
