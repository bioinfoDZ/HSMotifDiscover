

HSSeq_file=system.file("extdata", "Ex_simulated_HSseqs_new_M1c_and_M3c_S100E.txt",
                       package = "HSMotifDiscover", mustWork = TRUE)
charGrpFile=system.file("extdata", "Chars.txt", package = "HSMotifDiscover", mustWork = TRUE)
seq_weight_file=system.file("extdata", "motif1c_Weight.txt", package = "HSMotifDiscover", mustWork = TRUE)
motifLenVec=c(5,7)   # motif length range

#out=HSMotifDiscover(input_HSseq_file=HSSeq_file, motifLenVec=motifLenVec, charGrpFile=charGrpFile,
#                    seq_weight_file=seq_weight_file, numCores=1,  affinity_threshold=0,  itr=40000)
