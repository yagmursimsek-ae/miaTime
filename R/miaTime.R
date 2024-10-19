#' \code{miaTime} Package.
#'
#' \code{miaTime} implements time series methods in \link[mia]{mia} ecosystem.
#' @name miaTime-package
#' @aliases miaTime
#' @docType package
#' @seealso \link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment} class
NULL

#' @import mia
#' @import TreeSummarizedExperiment
NULL

#' @title HITChip Atlas with 1006 Western Adults
#' @description This data set contains genus-level microbiota profiling with
#' HITChip for 1006 western adults with no reported health complications,
#' reported in Lahti et al. (2014).
#' @name hitchip1006
#' @details The data is also available for download from the Data Dryad
#' \url{http://doi.org/10.5061/dryad.pk75d} and was initially released
#' as a phyloseq object under the name atlas1006 in \pkg{microbiome} R package.
#' The \pkg{miaTime} package provides this as
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object.
#' Some of the subjects also have short time series.
#' @docType data
#' @author Leo Lahti
#' @return Loads the data set in R.
#' @references
#' Lahti et al. Tipping elements of the human intestinal ecosystem.
#' Nature Communications 5:4344, 2014. \url{https://doi.org/10.1038/ncomms5344}.
#' @usage data(hitchip1006)
#' @format The data set in
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' format.
#' @keywords data
#'
NULL

#' @title Gut Microbiome Profiling 20 Belgian Women
#' @description This time-series data set contains absolute temporal variation
#' for most major gut genera, as well as diversity indicators for both
#' relative and quantitative abundance profiles of healthy women living in an
#' industrialized country.
#' @name temporalMicrobiome20
#' @details The data is available also on
#' \url{https://www.nature.com/articles/s41467-021-27098-7#Sec43}
#' The data contains 713 fecal samples from 20 Belgian women collected
#' over six weeks.
#' The current implementation in \pkg{miaTime} provides this as
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object.
#' @docType data
#' @author Doris Vandeputte
#' @return Loads the data set in R.
#' @references
#' Vandeputte, D., De Commer, L., Tito, R.Y. et al. Temporal variability
#' in quantitative human gut microbiome profiles and implications for clinical
#' research.
#' Nat Commun 12, 6740 (2021). https://doi.org/10.1038/s41467-021-27098-7
#' @usage data(temporalMicrobiome20)
#' @format The data set in
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' format.
#' @keywords data
#'
NULL

#' @title Human Gut Minimal Microbiome Profiling Data
#' @description This time-series data set contains absolute temporal abundances
#' for 16 human gut species that were assembled in an in vitro gut system. These
#' were subjected to a variety of disturbances over a period of 20 days. The
#' sample data includes measurements for Acetate, Butyrate, Propionate, and
#' optical density.
#' @name minimalgut
#' @details The data is available also on
#' \url{https://github.com/microsud/syncomR}
#' The data contains 183 samples from 3 in vitro gut system, 61 per bioreactor
#' collected three times daily.
#' The current implementation in \pkg{miaTime} provides this as
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' object.
#' @docType data
#' @author Sudarshan A. Shetty
#' @return Loads the data set in R.
#' @references
#' Shetty, S.A., Kostopoulos, I., Geerlings, S.Y. et al. Dynamic metabolic
#' interactions and trophic roles of human gut microbes identified using a
#' minimal microbiome exhibiting ecological properties.
#' ISME, 16, 2144–2159 (2022). https://doi.org/10.1038/s41396-022-01255-2.
#' @usage data(minimalgut)
#' @format The data set in
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' format.
#' @keywords data
#'
NULL

#' @title SilvermanAGutData
#' @description The SilvermanAGutData dataset contains 16S rRNA gene based high-throughput
#' profiling of 4 in vitro artificial gut models. The sampling was done hourly
#' and daily to capture sub-daily dynamics of microbial community originating
#' from human feces. The data consists of 413 taxa from 639 samples.
#' The data set can be used to investigate longitudinal dynamics of microbial
#' community in a controlled environment.
#'
#' Column metadata includes the days of sampling, vessel identifier, sampling
#' frequency pre-post challenge with Bacteroides ovatus.
#'
#' The row metadata of the microbiome data contains taxonomic information on
#' the Kingdom, Phylum, Class, Order, Family and Genus and Species level.
#'
#' The row tree consists of a phylogenetic tree build using sequence
#' information of 413 taxa.
#'
#' As reference sequences the ASV are provided.
#'
#' All data are downloaded from ExperimentHub and cached for local re-use
#'
#' @name SilvermanAGutData
#' @docType data
#' @author Sudarshan A. Shetty and Felix G.M. Ernst
#' @return Loads the data set in R.
#' @references
#' Silveman J.D et al. (2018): Dynamic linear models guide design and analysis
#' of microbiota studies within artificial human guts.
#' Microbiome 6:202 \url{https://doi.org/10.1186/s40168-018-0584-3}
#' @usage data(SilvermanAGutData)
#' @format The data set in
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' format.
#' @keywords data
#'
NULL
