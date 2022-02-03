#' @title Omics data downloading and preprocessing function
#' @param cohort Cohort abbreviation from Xenahub.
#' @param datatypes Datatypes to be downloaded and merged including "clinical" "molecular" "mRNA".
#' @param cancertype Specific tumor TCGA abbreviation (e.g., "LAML", "BRCA").
#' @importFrom magrittr "%>%"
#' @importFrom dplyr mutate
#' @importFrom stats setNames
#' @importFrom UCSCXenaTools XenaFilter
#' @importFrom UCSCXenaTools XenaGenerate
#' @importFrom UCSCXenaTools XenaDownload
#' @importFrom UCSCXenaTools XenaQuery
#' @importFrom UCSCXenaTools XenaPrepare
#' @importFrom RCurl getURL
#' @import dplyr
#' @export
#' @examples data <- Xenaprep(cancertype = "LAML")
NULL

Xenaprep <- function(cohort = "TCGA",
                      datatypes = c("clinical", "molecular", "mRNA"),
                          cancertype = "LAML"){

  XenaHostNames <- NULL
  subset_index <- UCSCXenaTools::XenaData$XenaHostNames == "pancanAtlasHub"
  Xena_host_PANCAN <- XenaGenerate(subset = subset_index) %>%
    XenaFilter(filterDatasets = "TCGASubtype.20170308.tsv") -> dir_molecsubtype_PANCAN

  #Download the molecular subtype dataset TCGASubtype.20170308.tsv
   suppressMessages(XenaQuery(dir_molecsubtype_PANCAN)) %>%
       XenaDownload() -> molecsubtype_PANCAN

   #Prepare the downloaded molecular-subtype data from TCGA PANCAN cohort for R.
  molecsubtype_PANCAN <- XenaPrepare(molecsubtype_PANCAN)
  if(identical("molecular", datatypes)){
    return(molecsubtype_PANCAN)
  }

  #Check if the code can access to xenahub database https://xenabrowser.net/
  if (is.character(getURL("https://xenabrowser.net/"))) {
    cat("Connecting\n")
  } else {
    cat("Cannot access to the Xenahub database.\n
          Check your internet connection and https://xenabrowser.net/\n")
  }

  create.list <- function(num.of.objects, object.names){
    list <- vector(mode = "list", length = num.of.objects)  %>%
      setNames(object.names)
    return(list)
  }
  message("Downloading the following datasets:")
  message(paste(paste(datatypes, sep  = " "), "data from",
                cancertype, "TCGA cohort // ", sep = " "))

  sampleID <- NULL #visible binding fix

  #Download other clinical covariates (check https://xenabrowser.net/")
  XenaGenerate(subset = XenaHostNames=="tcgaHub") %>%
    XenaFilter(filterDatasets = paste("survival/", cancertype, "_survival.txt",
                                      sep = '', collapse = NULL)) -> df_todo
  suppressMessages(XenaQuery(df_todo)) %>%
    XenaDownload() -> xe_download
  XenaPrepare(xe_download) %>%
    select_if(~sum(!is.na(.)) > 0) %>%
    mutate(sampleID = sample) %>%
    relocate(sampleID, .before = NULL, .after = NULL) -> clinical_TCGA

    #Download and prepare the first clinical information file
    XenaGenerate(subset = XenaHostNames=="tcgaHub") %>%
    XenaFilter(filterDatasets =  paste("TCGA.", cancertype, ".sampleMap/",
                                       cancertype, "_clinicalMatrix",
                                       sep = '', collapse = NULL)) -> df_todo
  XenaQuery(df_todo) %>%
    XenaDownload() -> xe_download
  XenaPrepare(xe_download) %>%
    select_if(~sum(!is.na(.)) > 0) -> clinical_TCGA2


  #Merge clinical and molecular subtype datasets (see "https://xenabrowser.net/").
  clinical_TCGA <- merge(clinical_TCGA, clinical_TCGA2, by = "sampleID")

  if("clinical" %in% datatypes &
     !("mRNA" %in% datatypes) &
     !("molecular" %in% datatypes)){
    return(clinical_TCGA)
  }

  merged_clinic.molecsubtype <- merge(clinical_TCGA, molecsubtype_PANCAN, by = "sampleID")
  #Download mean-normalized (per gene) mRNAseq PANCAN data normalized across all TCGA cohorts.
  XenaGenerate(subset = XenaHostNames=="tcgaHub") %>%
    XenaFilter(filterDatasets = paste("TCGA.", cancertype,
                                      ".sampleMap/HiSeqV2_PANCAN", sep = '', collapse = NULL)) -> query_data_PANCAN

  suppressMessages(XenaQuery(query_data_PANCAN)) %>%
    XenaDownload() -> download_data_PANCAN
  TmRNA_data_PANCANnorm <- XenaPrepare(download_data_PANCAN);

    #Prepare PANCAN normalized mRNA data for R analysis
  TmRNA_data_PANCANnorm <- as.data.frame(TmRNA_data_PANCANnorm)

  #These datasets need to be transposed to get samples to rows and genes to columns
  rownames(TmRNA_data_PANCANnorm) <- TmRNA_data_PANCANnorm[,1]
  TmRNA_data_PANCANnorm[,1] <- NULL
  mRNA_data_PANCANnorm <- as.data.frame(t(as.matrix(TmRNA_data_PANCANnorm)))
  mRNA_data_PANCANnorm %>%
    mutate(sampleID = rownames(mRNA_data_PANCANnorm),
           sample = rownames(mRNA_data_PANCANnorm)) %>%
    relocate(sample, .before = NULL, .after = NULL) %>%
    relocate(sampleID, .before = NULL, .after = NULL) -> mRNA_data_PANCANnorm
  rm(TmRNA_data_PANCANnorm)


  merged_clinic.molecsubtype.mRNA <- merge(merged_clinic.molecsubtype,
                                           mRNA_data_PANCANnorm, by = "sampleID")
  merged_clinic.mRNA <- merge(clinical_TCGA, mRNA_data_PANCANnorm, by = "sampleID")


  if("clinical" %in% datatypes & !("mRNA" %in% datatypes) &
     "molecular" %in% datatypes){
    return(merged_clinic.molecsubtype)

  } else if (("clinical" %in% datatypes & "mRNA" %in% datatypes &
              "molecular" %in% datatypes) |("clinical" %in% datatypes &
                                            "mRNA" %in% datatypes &
                                            !("molecular" %in% datatypes)) ){

    #Return a list consisting of merged datasets
    output_datasets <- create.list(3, c("clin", "clin.molsubtype", "clin.molsubtype.mRNA"))

    output_datasets$clin <- clinical_TCGA
    output_datasets$clin.molsubtype <- merged_clinic.molecsubtype
    output_datasets$clin.molsubtype.mRNA <- merged_clinic.molecsubtype.mRNA
  }

  return(output_datasets)

}
