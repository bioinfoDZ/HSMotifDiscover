#' @title  HSMotifDiscover
#'
#' @description  Discover motif in heparan sulphate sequences or any type of sequnces.
#' @param input_HSseq_file A heparin sulphate sequence (or any other seqeunce) file. The file contains header and sequence information of the samples similar to fasta file, but sequence character counts per line are not constrained.
#' @param motifLenVec Describe a vector of the motif lenghts to be discovered.
#' @param charGrpFile (Optional parameter) File having two columns , \emph{One}- dimers/trimers/tetramers that should be considered as single characters to discover motif,  \emph{Two}- character group information in numeric form. In heparin sulphate dimer and trimer occupy alternate positions are. So, dimers are grouped in one group and trimers are grouped in other group. If this file is not given then each character will be considered independently to discover motif as other motif discovery tools such as MEME, but will work for any other type of characters.
#' @param seq_weight_file  (Optional parameter) File having column of sequence header and sequence weight in motif discovery. If this file is not selected then all sequences have equal weight in motif discovery.
#' @param numCores  (Optional parameter) the number of cores to be used. This feature is useful when motif is discovered in large range. Multiple cores feature will not work on widows computer or RStudio environment.
#' @param affinity_threshold  (Optional parameter) The sequences with weight greater than the threshold are used for motif discovery. This is required only of weight file is given.
#' @param itr  (Optional parameter) Number of iterations for gibbs sampling optimisation. Higher itrations may improve the results but at the cost of time.
#' @return \enumerate{
#' \item \emph{list} A list of following parameters of all the discovered motifs in the given range.
#' \enumerate{
#'   \item  \emph{PSSM:} Position specific scoring matrix of the discovered motif.
#'   \item  \emph{MotifEntropy:} Average entropy of discovered motif.
#'   \item  \emph{IC:} Information content of the discovered motif.
#'   \item  \emph{bkgProb_vec:} Probability vector of the characters in the background sequences.
#'   \item  \emph{motif_Pval:} P-Value of the discovered motif.
#'   \item  \emph{CharDataf:} Summary of the input sequences characters.
#'   \item  \emph{orgChar_discoveredMotifs:} Sequences of the discovered motif.
#'   \item  \emph{MotifLogo:} Discovered motif seqlogo.
#'   \item  \emph{resultsTable_df:} Results of the discovered motif are summarised in the data frame. Where columns have following information
#'   \itemize{
#'      \item \emph{header:} header of the sequences used for motif discovery.
#'      \item \emph{sequence:} Input sequences.
#'      \item \emph{seqLen:} sequence length.
#'      \item \emph{if sequence weight is greater than the threshold:} boolean (TRUE/FALSE). This column is present only of weight of the sequences is given as input.
#'      \item \emph{max_score:} Maximum likelihood score of the of the PSSM in the sequnce.
#'      \item \emph{P_value:} The probability that a random sequence (with the same length and conforming to the background) would have position p-values such that the product is smaller or equal to the value calculated for the sequence under test.
#'      \item \emph{MotifLoc:} Start point of discovered motif. (This is the location of max_score in the sequence )
#'      \item \emph{seq_weights:} Input weights of the sequnces. This column is present only if weight of the sequences is given as input.
#'      \item \emph{score:} Liklihood score of the PSSM at different sequence locations.
#'  }
#'  }
#'  \item \emph{MotifSummary_runTime_*.txt:} Summary of motifs in the text file.
#'  \item \emph{log_runTime_*.txt:} log file of the run.
#'  }
#' @details Discover motif in heparan sulphate sequences or any type of sequnces.
#' @export
#' @examples
#' HSSeq_file=system.file("extdata", "Ex_simulated_HSseqs_new_M1c_and_M3c_S100E.txt",
#' package = "HSMotifDiscover", mustWork = TRUE)
#' charGrpFile=system.file("extdata", "Chars.txt", package = "HSMotifDiscover", mustWork = TRUE)
#' seq_weight_file=system.file("extdata", "motif1c_Weight.txt", package = "HSMotifDiscover", mustWork = TRUE)
#' motifLenVec=c(5,7)   # motif length range
#'
#' out=HSMotifDiscover(input_HSseq_file=HSSeq_file, motifLenVec=motifLenVec, charGrpFile=charGrpFile,
#' seq_weight_file=seq_weight_file, numCores=1,  affinity_threshold=0,  itr=40000)
#'
#'

#
# @seealso \code{\link{vcfToSNV}} to generate \emph{snv} dataframe, and \url{https://genome.ucsc.edu/FAQ/FAQformat.html} for BED file format.
#




HSMotifDiscover <- function(input_HSseq_file, motifLenVec,  charGrpFile, seq_weight_file, numCores=1,  affinity_threshold=0,  itr=40000)
{
  RunFileExt=paste0(format(Sys.time(), "%d%b%Y%H%M%S"),'.txt')

  sink(paste0(getwd(),'/','log_runTime_',RunFileExt))

    print(' in > HSMotifDiscover')



    isRStudio <- Sys.getenv("RSTUDIO") == "1"

    if((Sys.info()[['sysname']]== 'Linux') | (Sys.info()[['sysname']]== 'Darwin')){
      os_Flag= 'good'
    }

    #motifLenVec=seq(motif_range[1], motif_range[2])

      if(missing(input_HSseq_file) & missing(motifLenVec)){
        geterrmessage("input sequnce file and motif length/ range can't be left empty")
      }
      if (missing(seq_weight_file) & missing(charGrpFile))
      {
        print('>>>C1')
        if((os_Flag=='good') & (numCores > 1) & (isRStudio== FALSE)){
          motif_range_data_results=mclapply(motifLenVec , function(x) CompileResults(input_HSseq_file=input_HSseq_file, motif_length=x, affinity_threshold=affinity_threshold,  itr=itr), mc.cores = getOption("mc.cores", 2L))
        }else(
          motif_range_data_results=lapply(motifLenVec , function(x) CompileResults(input_HSseq_file=input_HSseq_file, motif_length=x,  affinity_threshold=affinity_threshold,  itr=itr))
        )
      }
      if(missing(seq_weight_file) & !missing(charGrpFile))
      {
        print('>>>C2')
        if((os_Flag=='good') & (numCores > 1) & (isRStudio== FALSE)){
          motif_range_data_results=mclapply(motifLenVec , function(x) CompileResults(input_HSseq_file=input_HSseq_file, motif_length=x, charGrpFile=charGrpFile, affinity_threshold=affinity_threshold, itr=itr),mc.cores = getOption("mc.cores", 2L) )
        }else(
          motif_range_data_results=lapply(motifLenVec , function(x) CompileResults(input_HSseq_file=input_HSseq_file, motif_length=x,  charGrpFile=charGrpFile, affinity_threshold=affinity_threshold,  itr=itr))
        )
      }
      if(!missing(seq_weight_file) & missing(charGrpFile))
      {
        print('>>>C3')
        if((os_Flag=='good') & (numCores > 1) & (isRStudio== FALSE)){
          motif_range_data_results=mclapply(motifLenVec , function(x) CompileResults(input_HSseq_file=input_HSseq_file, motif_length=x,  seq_weight_file=seq_weight_file, affinity_threshold=affinity_threshold,  itr=itr), mc.cores = getOption("mc.cores", 2L))
        }else(
          motif_range_data_results=lapply(motifLenVec , function(x) CompileResults(input_HSseq_file=input_HSseq_file, motif_length=x, seq_weight_file=seq_weight_file, affinity_threshold=affinity_threshold,  itr=itr))
        )
      }
      if(!missing(seq_weight_file) & !missing(charGrpFile))
      {
        print('>>>C4')
        if((os_Flag=='good') & (numCores > 1) & (isRStudio== FALSE)){
          motif_range_data_results=mclapply(motifLenVec , function(x) CompileResults(input_HSseq_file=input_HSseq_file, motif_length=x,  charGrpFile=charGrpFile,  seq_weight_file=seq_weight_file, affinity_threshold=affinity_threshold, itr=itr),  mc.cores = getOption("mc.cores", 2L))
        }else(
          motif_range_data_results=lapply(motifLenVec , function(x) CompileResults(input_HSseq_file=input_HSseq_file, motif_length=x, charGrpFile=charGrpFile,  seq_weight_file=seq_weight_file, affinity_threshold=affinity_threshold,  itr=itr))
        )
      }

      print(' out > HSMotifDiscover')

    sink()


    sink(paste0('MotifSummary_runTime_',RunFileExt))

      for(j in 1:length(motifLenVec))
      {
        motif_range_data_results[[j]]=motif_range_data_results[[j]]$motifData
        print(paste0('Motif Length: ',motifLenVec[j] ))
        cat('\n=================\n')
        cat('Position Specific Scoring Matrix: \n')
        print(motif_range_data_results[[j]]$PSSM)
        cat('\nMotif statistics: \n')
        cat(paste0('Motif Entropy: ', motif_range_data_results[[j]]$MotifEntropy,', Information Content: ',motif_range_data_results[[j]]$IC,', Motif P-value: ',signif(motif_range_data_results[[j]]$motif_Pval, 100)))
        cat('\n\n')
        cat('\n')

      }

      names(motif_range_data_results)=paste0('MotifLength_',motifLenVec)

    sink()




  return(motif_range_data_results)
}
