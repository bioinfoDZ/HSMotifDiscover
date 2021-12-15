summary_motif_data <- function(motif_data, if_revCompStrand)
{

    print('in > summary_motif_data')
    counts=table(nchar(motif_data$orgChar_discoveredMotifs))
    most_freq_char=names(counts)[which.max(counts)]
    des_motifs_ind=which(nchar(motif_data$orgChar_discoveredMotifs)==most_freq_char)

    str1=paste0(motif_data$CharData$orginal_char,collapse = '')
    uniqchars <- function(x) unique(strsplit(x, "")[[1]])

    UC=uniqchars(str1)
    #rand_colors <- randomColor(length(UC), luminosity="bright")  #library(randomcoloR)
    c28 <- c(
        "green4","blue1", "#E31A1C", # red
        "gold1","black",
        "#6A3D9A", # purple
        "#FF7F00", # orange
         "#555555",
        "maroon", "orchid1", "deeppink1", "dodgerblue2", "steelblue4",
        "darkturquoise", "green1", "yellow4", "yellow3","#225555",
        "darkorange4", "brown",
        "#364B9A",
        "skyblue2", "#FB9A99", # lt pink
        "palegreen2",
        "#CAB2D6", # lt purple
        "#FDBF6F", # lt orange
        "gray70", "khaki2"
    )
   # pie(rep(1, 28), col = c28)

    rand_colors=c28[seq(1,length(UC))]
    cs1 = make_col_scheme(chars=UC,  cols=rand_colors)
    p1= ggseqlogo(motif_data$orgChar_discoveredMotifs[des_motifs_ind], namespace=UC,  method='b', col_scheme=cs1) + ylim(0,log2(length(UC)))
    
   
    motif_data$motif_DesiredSeqsName=motif_data$seqs[des_motifs_ind]
    seqsT=BStringSet(motif_data$orgChar_discoveredMotifs[des_motifs_ind])
    unique_letters=uniqueLetters(seqsT)
    FreqMat=sapply(seq(1:width(seqsT)[1]),function(x) colSums(letterFrequency(subseq(seqsT, start=x, width=1),letters=unique_letters)))
    PSSM=sweep(FreqMat,2,colSums(FreqMat),'/')

    motif_data$PSSM=PSSM
    motif_data$MotifEntropy=Entropy(PSSM)     #sum(sapply(seq(1,ncol(PSSM)), function(x) Entropy(PSSM[,x])))
    motif_data$IC=signif(log2(length(unique_letters))-mean(sapply(seq(1,ncol(PSSM)), function(x) Entropy(PSSM[,x]))),3)
    motif_data$MotifLogo=p1
    motif_data$ReverseComplement_PSSM=reverseComplement(PSSM)
    
    if(if_revCompStrand==TRUE){
        rc_p1= ggseqlogo(as.vector(reverseComplement(DNAStringSet(motif_data$orgChar_discoveredMotifs[des_motifs_ind]))), namespace=UC,  method='b', col_scheme=cs1)+ ylim(0,2)
        motif_data$ReverseComplement_MotifLogo=rc_p1
    }
    

    mapC_p1=ggseqlogo(motif_data$mapChar_discoveredMotifs[des_motifs_ind], namespace=motif_data$CharData$mapped_char, method='b')
    mapC_seqsT=BStringSet(motif_data$mapChar_discoveredMotifs[des_motifs_ind])
    mapC_unique_letters=uniqueLetters(mapC_seqsT)
    mapC_FreqMat=sapply(seq(1:width(mapC_seqsT)[1]),function(x) colSums(letterFrequency(subseq(mapC_seqsT, start=x, width=1),letters=mapC_unique_letters)))
    mapC_PSSM=sweep(mapC_FreqMat,2,colSums(mapC_FreqMat),'/')


    motif_data$mapC_PSSM=mapC_PSSM
    motif_data$mapC_MotifEntropy=Entropy(mapC_PSSM)     #sum(sapply(seq(1,ncol(PSSM)), function(x) Entropy(PSSM[,x])))

    motif_data$mapC_IC=signif(log2(length(mapC_unique_letters))-mean(sapply(seq(1,ncol(mapC_PSSM)), function(x) Entropy(mapC_PSSM[,x]))),3)     #sum(sapply(seq(1,ncol(PSSM)), function(x) Entropy(PSSM[,x])))
    motif_data$mapC_MotifLogo=mapC_p1
    
    motif_data$motif_Obj=UVmotif_Obj(motif_data$mapC_PSSM, motif_data$final_bkg_vec, CharDataf=motif_data$CharData )

    print('out > summary_motif_data')
    
    return(motif_data)
}



