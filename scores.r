rm(list=ls())
library(Biostrings)
nRegion <- 600
library(stringr)
args = commandArgs(trailingOnly=TRUE)


#-------------------------------------read matrices ------------------------------------------


pfm_ARF5 <- read.table("m_ARF5.txt",header=TRUE,sep="\t",skip=1)
pfm_ARF5 <- rbind(0,0,pfm_ARF5) 
pfm_ARF5 <- round((t(as.matrix(pfm_ARF5)))*nRegion)+1 #;pfm_ARF5
maxi_ARF5 <- apply(pfm_ARF5,FUN=max, 2)
maxi_ARF5 <- matrix(nrow=4, rep(maxi_ARF5,4),byrow=TRUE)
pwm_ARF5 <- log(pfm_ARF5/maxi_ARF5)
pwm_ARF5_rev <- pwm_ARF5 
#write.table(t(pwm_ARF5_rev),"PWM_ARF2.txt",quote=FALSE,sep='\t',row.names=FALSE)
pwm_ARF5 <-  reverseComplement(pwm_ARF5_rev) #; pwm_ARF5



pfm_ARF2 <- read.table("m_ARF2.txt",header=TRUE,sep="\t",skip=1)
pfm_ARF2 <- rbind(0,pfm_ARF2,0) 
pfm_ARF2 <- round((t(as.matrix(pfm_ARF2)))*nRegion)+1 #;pfm_ARF2
maxi_ARF2 <- apply(pfm_ARF2,FUN=max, 2)
maxi_ARF2 <- matrix(nrow=4, rep(maxi_ARF2,4),byrow=TRUE)
pwm_ARF2 <- log(pfm_ARF2/maxi_ARF2)
pwm_ARF2_rev <- pwm_ARF2 
#write.table(t(pwm_ARF2_rev),"PWM_ARF2.txt",quote=FALSE,sep='\t',row.names=FALSE)
pwm_ARF2 <-  reverseComplement(pwm_ARF2_rev) #; pwm_ARF2


#-------------------------------------read fasta files-----------------------------------------


ARF5_pos <- readDNAStringSet(args[1])
width_pos <- width(ARF5_pos)
seq_pos <- as.character(ARF5_pos)



#-------------------------------------Compute Scores-----------------------------------------
tab_ARF2 <- NULL
tab_ARF5 <-  NULL
k <- 0
if(len(args == 3))
{
    threshold <- args[3]
}
else
{
    threshold <- -12
}

while(k < length(seq_pos))
{
    k <- k+1
    p <- 0
    for (PWM in list(pwm_ARF5,pwm_ARF2))
    {
        p <- p+1
        PWM_rev <- reverseComplement(PWM)
        ARF <- list("ARF5","ARF2")
        match_ARF5 <- sapply(FUN=matchPWM,seq_pos,pwm=PWM,min.score=-12)
        scores_ARF5 <- mapply(seq_pos,FUN=PWMscoreStartingAt,starting.at=sapply(FUN=start,match_ARF5[]),SIMPLIFY=FALSE,MoreArgs=list(pwm=PWM))
        
        match_ARF5_rev<- sapply(FUN=matchPWM,seq_pos,pwm=PWM_rev,min.score=-12)
        scores_ARF5_rev<- mapply(seq_pos,FUN=PWMscoreStartingAt,starting.at=sapply(FUN=start,match_ARF5_rev[]),SIMPLIFY=FALSE,MoreArgs=list(pwm=PWM_rev))

        sequence <- as.character(ARF5_pos[k])
        promoter <- names(sequence)

        sites <- start(match_ARF5[[k]])
        sites_rev<- start(match_ARF5_rev[[k]])

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
                    DR <- rbind(DR,c((elt2 - elt1), sites[i],scores_ARF5[[k]][i],scores_ARF5[[k]][j],str_sub(sequence,sites[i] +3, sites[i] +3 + 6 + (elt2 - elt1) -1)))
                }
            }
        }
        if(!is.null(DR))
        {
            rownames(DR) <- rep("DR",dim(DR)[1])
            colnames(DR) <- c("spacing","position","score1","score2","sequence")
        }
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
                    DR_rev <- rbind(DR_rev,c((elt2 - elt1), sites_rev[i],scores_ARF5_rev[[k]][i],scores_ARF5_rev[[k]][j],str_sub(sequence,sites_rev[i] +3, sites_rev[i] +3 + 6 + (elt2 - elt1) -1)))
                }
            }
        }
        if(!is.null(DR_rev))
        {
            rownames(DR_rev) <- rep("DR rev",dim(DR_rev)[1])
            colnames(DR_rev) <- c("spacing","position","score1","score2","sequence")
        }
        
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
                    ER <- rbind(ER,c((elt2 - elt1), sites[i],scores_ARF5[[k]][i],scores_ARF5_rev[[k]][j],str_sub(sequence,sites[i] +3, sites[i] +3 + 6 + (elt2 - elt1) -1)))
                }
            }
        }
        if(!is.null(ER))
        {
            rownames(ER) <- rep("ER",dim(ER)[1])
            colnames(ER) <- c("spacing","position","score1","score2","sequence")
        }

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
                    IR <- rbind(IR,c((elt2 - elt1), sites_rev[i],scores_ARF5_rev[[k]][i],scores_ARF5[[k]][j],str_sub(sequence,sites_rev[i] +3, sites_rev[i] +3 + 6 + (elt2 - elt1) -1)))
                }
            }
        }
        if(!is.null(IR))
        {
            rownames(IR) <- rep("IR",dim(IR)[1])
            colnames(IR) <- c("spacing","position","score1","score2","sequence")
        }
        
        tab <- rbind(DR,DR_rev,ER,IR)
        if(!is.null(tab))
        {
            tab[,1] <- as.integer(tab[,1]) - 6 
            tab[,2] <- as.integer(tab[,2]) + 3 
            tab <- tab[as.integer(tab[,1])> 0 , ,drop=FALSE]
            tab[,3] <- round(as.integer(tab[,3]))
            tab[,4] <- round(as.integer(tab[,4]))
            tab <- tab[order(as.integer(tab[,2])),]
            if(length(args)==1)
            {
                write.table(tab,paste("Interdistances_",ARF[p],"_",promoter,".csv",sep=""),col.names=NA,quote=FALSE,sep="\t")
            }
            else
            {
                if(ARF[p]=='ARF2')
                {
                   tab_ARF2 <- rbind(tab_ARF2,rep(promoter,5))
                   tab_ARF2 <- rbind(tab_ARF2,tab)
                }
                if(ARF[p]=='ARF5')
                {
                   tab_ARF5 <- rbind(tab_ARF5,rep(promoter,5))
                   tab_ARF5 <- rbind(tab_ARF5,tab)
                }
            }
        }
    }
}

if(length(args)!=1)
{
    write.table(tab_ARF2,paste("Interdistances_ARF2_",args[2],".csv",sep=""),col.names=NA,quote=FALSE,sep="\t")
    write.table(tab_ARF5,paste("Interdistances_ARF5_",args[2],".csv",sep=""),col.names=NA,quote=FALSE,sep="\t")
}
