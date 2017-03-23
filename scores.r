rm(list=ls())
library(Biostrings)
nRegion <- 600
library(stringr)


#-------------------------------------read matrices ------------------------------------------


pfm_ARF5 <- read.table("m_ARF5.txt",header=TRUE,sep="\t",skip=1)
pfm_ARF5 <- rbind(pfm_ARF5,0,0) 
pfm_ARF5 <- round((t(as.matrix(pfm_ARF5)))*nRegion)+1 ;pfm_ARF5
maxi_ARF5 <- apply(pfm_ARF5,FUN=max, 2)
maxi_ARF5 <- matrix(nrow=4, rep(maxi_ARF5,4),byrow=TRUE)
pwm_ARF5 <- log(pfm_ARF5/maxi_ARF5)
pwm_ARF5_rev <- pwm_ARF5 
#write.table(t(pwm_ARF5_rev),"PWM_ARF2.txt",quote=FALSE,sep='\t',row.names=FALSE)
pwm_ARF5 <-  reverseComplement(pwm_ARF5_rev) ; pwm_ARF5



pfm_ARF2 <- read.table("m_ARF2.txt",header=TRUE,sep="\t",skip=1)
pfm_ARF2 <- rbind(0,pfm_ARF2,0) 
pfm_ARF2 <- round((t(as.matrix(pfm_ARF2)))*nRegion)+1 ;pfm_ARF2
maxi_ARF2 <- apply(pfm_ARF2,FUN=max, 2)
maxi_ARF2 <- matrix(nrow=4, rep(maxi_ARF2,4),byrow=TRUE)
pwm_ARF2 <- log(pfm_ARF2/maxi_ARF2)
pwm_ARF2_rev <- pwm_ARF2 
#write.table(t(pwm_ARF2_rev),"PWM_ARF2.txt",quote=FALSE,sep='\t',row.names=FALSE)
pwm_ARF2 <-  reverseComplement(pwm_ARF2_rev) ; pwm_ARF2


#-------------------------------------read fasta files-----------------------------------------


ARF5_pos <- readDNAStringSet('promoters_dofs.fasta')
width_pos <- width(ARF5_pos)
seq_pos <- as.character(ARF5_pos)



#-------------------------------------Compute Scores-----------------------------------------
k <- 0
for (PWM in list(pwm_ARF5,pwm_ARF2))
{
    PWM_rev <- reverseComplement(PWM)
    k <- k + 1
    ARF <- list("ARF5","ARF2")
    match_ARF5 <- sapply(FUN=matchPWM,seq_pos,pwm=PWM,min.score=-12)
    scores_ARF5 <- mapply(seq_pos,FUN=PWMscoreStartingAt,starting.at=sapply(FUN=start,match_ARF5[]),SIMPLIFY=FALSE,MoreArgs=list(pwm=PWM))
    
    match_ARF5_rev<- sapply(FUN=matchPWM,seq_pos,pwm=PWM_rev,min.score=-12)
    scores_ARF5_rev<- mapply(seq_pos,FUN=PWMscoreStartingAt,starting.at=sapply(FUN=start,match_ARF5_rev[]),SIMPLIFY=FALSE,MoreArgs=list(pwm=PWM_rev))


    promoter <- "pDOF58"

    sites <- start(match_ARF5$pDOF58)
    sites_rev<- start(match_ARF5_rev$pDOF58)

    DR <- NULL
    i <- 0
    for (elt1 in sites)   
    {
        i <- i+1
        j <- 0
        for (elt2 in sites)
        {
            j <- j+1  
            if ((elt2 - elt1) < 20 && (elt2 - elt1) >0)
            {
                DR <- rbind(DR,c((elt2 - elt1), sites[i],scores_ARF5$pDOF58[i],scores_ARF5$pDOF58[j]))
            }
        }
    }
    rownames(DR) <- rep("DR",dim(DR)[1])
    colnames(DR) <- c("spacing","position","score1","score2")
    
    DR_rev <- NULL
    i <- 0
    for (elt1 in sites_rev)   
    {
        i <- i+1
        j <- 0
        for (elt2 in sites_rev)
        {
            j <- j+1  
            if ((elt2 - elt1) < 20 && (elt2 - elt1) >0)
            {
                DR_rev <- rbind(DR_rev,c((elt2 - elt1), sites_rev[i],scores_ARF5_rev$pDOF58[i],scores_ARF5_rev$pDOF58[j]))
            }
        }
    }
    rownames(DR_rev) <- rep("DR rev",dim(DR_rev)[1])
    colnames(DR_rev) <- c("spacing","position","score1","score2")
    
    
    ER <- NULL
    i <- 0
    for (elt1 in sites)   
    {
        i <- i+1
        j <- 0
        for (elt2 in sites_rev)
        {
            j <- j+1  
            if ((elt2 - elt1) < 20 && (elt2 - elt1) >0)
            {
                ER <- rbind(ER,c((elt2 - elt1), sites[i],scores_ARF5$pDOF58[i],scores_ARF5_rev$pDOF58[j]))
            }
        }
    }
    rownames(ER) <- rep("ER",dim(ER)[1])
    colnames(ER) <- c("spacing","position","score1","score2")


    IR <- NULL
    i <- 0
    for (elt1 in sites_rev)   
    {
        i <- i+1
        j <- 0
        for (elt2 in sites)
        {
            j <- j+1  
            if ((elt2 - elt1) < 20 && (elt2 - elt1) >0)
            {
                IR <- rbind(IR,c((elt2 - elt1), sites_rev[i],scores_ARF5_rev$pDOF58[i],scores_ARF5$pDOF58[j]))
            }
        }
    }
    rownames(IR) <- rep("IR",dim(IR)[1])
    colnames(IR) <- c("spacing","position","score1","score2")
    
    tab <- rbind(DR,DR_rev,ER,IR)
    tab[,1] <- tab[,1] - 6 
    tab[,2] <- tab[,2] + 3 
    tab <- tab[tab[,1]> 0 , ]
    write.table(round(tab,5),paste("Interdistances_",ARF[k],"_DOF58",".csv",sep=""),col.names=NA,quote=FALSE,sep="\t")
}
