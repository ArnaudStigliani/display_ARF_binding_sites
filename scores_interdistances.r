rm(list=ls())
library(Biostrings)
nRegion <- 600
library(stringr)

#-------------------------------------read matrices ------------------------------------------


pfm_ARF5<- read.table("m_ARF5.txt",header=TRUE,sep="\t",skip=1)
pfm_ARF5 <- round((t(as.matrix(pfm_ARF5)))*nRegion)+1 ;pfm_ARF5
pfm_ARF5 <- cbind(rep(151,4),rep(151,4),pfm_ARF5)
maxi_ARF5 <- apply(pfm_ARF5,FUN=max, 2)
maxi_ARF5 <- matrix(nrow=4, rep(maxi_ARF5,4),byrow=TRUE)
pwm_ARF5 <- log(pfm_ARF5/maxi_ARF5)
pwm_ARF5_rev <- pwm_ARF5 
pwm_ARF5 <-  reverseComplement(pwm_ARF5_rev) ; pwm_ARF5

#-------------------------------------read fasta files-----------------------------------------


ARF5_pos <- readDNAStringSet('promoters_dofs.fasta')

width_pos <- width(ARF5_pos)


seq_pos <- as.character(ARF5_pos)



#-------------------------------------Compute Scores-----------------------------------------

#pos

scores_ARF5_pos<- mapply(seq_pos,FUN=PWMscoreStartingAt,SIMPLIFY=FALSE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF5)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF5))

scores_ARF5_rev_pos <- mapply(seq_pos,FUN=PWMscoreStartingAt,SIMPLIFY=FALSE,  starting.at=mapply(seq,1,width_pos-dim(pwm_ARF5_rev)[2],SIMPLIFY=FALSE),MoreArgs=list(pwm=pwm_ARF5_rev))

#-----------------------------------Compute Scores Interdistances-----------------------------

#initialise

scores_DR_pos <- NULL
scores_DR_pos_rev <-  NULL
scores_ER_pos <- NULL
scores_IR_pos <- NULL



# Compute DR 

for(i in 1:length(ARF5_pos))
{
    scores_DR_pos[[i]] <- list()
    for (j in 7:27)
    {
        scores_DR_pos[[i]][[j-6]] <- scores_ARF5_pos[[i]][j:(width_pos[i]-dim(pwm_ARF5)[2])] + scores_ARF5_pos[[i]][1:(width_pos[i]-j+1-dim(pwm_ARF5)[2])]
    }
}

DR <- lapply(FUN=which,lapply(FUN= ">",unlist(scores_DR_pos,recursive=FALSE),-14))


# Compute DR rev

for(i in 1:length(ARF5_pos))
{
    for (j in 7:27)
    {
        scores_DR_pos_rev[[j-6]] <- scores_ARF5_rev_pos[[i]][j:(width_pos[i]-dim(pwm_ARF5_rev)[2])] + scores_ARF5_rev_pos[[i]][1:(width_pos[i]-j+1-dim(pwm_ARF5_rev)[2])]
    }
}

DR_rev <- lapply(FUN=which,lapply(FUN= ">",scores_DR_pos_rev,-20))

# Best score between DR and DR rev

#scores_DR_pos <- ifelse(scores_DR_pos > scores_DR_pos_rev, scores_DR_pos,scores_DR_pos_rev)

# Compute IR 

for(i in 1:length(ARF5_pos))
{
    for (j in 7:27)
    {
        scores_IR_pos[[j-6]] <- (scores_ARF5_pos[[i]][j:(width_pos[i]-dim(pwm_ARF5)[2])] + scores_ARF5_rev_pos[[i]][1:(width_pos[i]-j+1-dim(pwm_ARF5_rev)[2])])
    }
}

IR <- lapply(FUN=which,lapply(FUN= ">",scores_IR_pos,-20))

# Compute ER

for(i in 1:length(ARF5_pos))
{
    for (j in 7:27)
    {
        scores_ER_pos[[j-6]] <- (scores_ARF5_rev_pos[[i]][j:(width_pos[i]-dim(pwm_ARF5_rev)[2])] + scores_ARF5_pos[[i]][1:(width_pos[i]-j+1-dim(pwm_ARF5)[2])])
    }
}


ER  <- lapply(FUN=which,lapply(FUN= ">",scores_ER_pos,-20))


#------------------------Meilleurs scores--------------------------------




## pn01 <- c(rep(0,length(ARF2_pos)),rep(1,length(ARF2_neg)))
    
## lm1 <- glm(pn01~rbind(scores_ER_pos,scores_ER_neg),family=binomial)   # variable 1

## an1 <- anova(lm1,test="Chisq")
