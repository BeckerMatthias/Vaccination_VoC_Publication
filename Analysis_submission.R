#This code is used to generate the graphs used in the manuscript
#In general graphs are exported to an output folder in .pdf or .svg format, 
#from where they are further processed into the finalized figures using a vector graphics editing software such as InkScape

#R version used was version 3.6.1 (2019-07-05) -- "Action of the Toes" within RStudio 1.2.5001

#Required are an input folder at "/input" with the file containing the raw measurement data 
#"Raw_Data_for_analysis.csv" with annotations, the raw data for the parallelism and Quality Control performance 
#"Parallelism.csv.txt", "QC_IgG.csv.txt" and "QC_IgA.csv.txt"

#All export commands such as "write.table" or "pdf" have been disabled by addition of "#" at the start of the line.
#Plots will therefore appear in RStudio instead of exporting and tables will be in the environment
#Each subsection corresponds to one figure in the manuscript

#####Library#####

#used for color sets
library(RColorBrewer)

#used for specific depictions
library(beeswarm)

#####
#####Data read-in####

#MULTICOV-AB Data, this includes the quality control samples CO1,CO2,CO3 and a Blank which were processed on every assay plate
MFI <- read.csv("Input/MULTICOV-AB_Results.csv",h=T)


#MULTICOV-AB Saliva Measurement Data
S <- read.csv("Input/MULTICOV-AB_Saliva_Results.csv",h=T)


#NeutrobodyPlex (NBP) Assay Results
N <- read.csv("Input/NBP_Results.csv",h=T)


#ACE2 Neutralization Assay Results
A <- read.csv("Input/ACE2_Results.csv",h=T)


#Virus Neutralisation Assay (VNT) Results
VNT <- read.csv("Input/VNT_Results.csv",h=T)




#####
#####Data normalisation and processing####

#MULTICOV-AB Data

#Convert dT (Time between Infection or Vaccination 1/2) from factor to numeric
MFI[,"dT1"] <- as.numeric(as.character(MFI[,"dT1"]))
MFI[,"dT2"] <- as.numeric(as.character(MFI[,"dT2"]))
MFI[,"dTinf"] <- as.numeric(as.character(MFI[,"dTinf"]))

#Normalize MFI signal to control sample on plates, use CO 2 for IgG Detection and CO 3 for IgA

M <- MFI

for (j in names(M)[10:39]){
  for (i in 1:nrow(M)){
    if (M[i,"SampleID"]%in%c("CO 1", "CO 2", "CO 3", "Blank", "NC1", "NC2")){
    }else if (M[i,"Assay"]%in%c("IgG")){
      M[i,j] <- M[i,j]/mean(M[M$Plate==M[i,"Plate"]&M$SampleID=="CO 2",j])
    }else if (M[i,"Assay"]=="IgA"){
      M[i,j] <- M[i,j]/mean(M[M$Plate==M[i,"Plate"]&M$SampleID=="CO 3",j])
    }
  }
}
rm(i,j)

#Remove Control Samples
M <- M[!M$SampleID%in%c("CO 1", "CO 2", "CO 3", "Blank", "NC1", "NC2"),]
MFI <- MFI[!MFI$SampleID%in%c("CO 1", "CO 2", "CO 3", "Blank", "NC1", "NC2"),]


#MULTICOV-AB Saliva Data

#Remove longitudinal series samples from Vaccinated individuals (only one time point used in MULTICOV-AB Saliva)
S <- S[!S$Sample_id%in%c(28,79,80,82,83),]



#ACE2 and NBP Data

#Convert dT from factor to numeric
A[,"dT1"] <- as.numeric(as.character(A[,"dT1"]))
A[,"dT2"] <- as.numeric(as.character(A[,"dT2"]))
N[,"dT1"] <- as.numeric(as.character(N[,"dT1"]))
N[,"dT2"] <- as.numeric(as.character(N[,"dT2"]))
A[,"dTinf"] <- as.numeric(as.character(A[,"dTinf"]))
N[,"dTinf"] <- as.numeric(as.character(N[,"dTinf"]))


#Annotate VNT Data with Meta Data by matching sample IDs
VNT <- merge(VNT,A[,1:6],by="SampleID",all.x=T,sort=F)


#####
#####Manuscript Figure 1####


mycol <- c("#ff4000","#00bfff","#004d67","#7F7F7F")

#a) - N vs RBD for MULTICOV-AB IgG highlighting different sample groups
#svg(paste("Output/Fig1/Fig1_a.svg",sep=""),7.5,5)
par(mfrow=c(1,1),mar=c(4,4,3,1),pty="m")
tmp <- M[M$Assay=="IgG",]
plot(tmp$SARS2.RBD,tmp$SARS2.N,log="xy",cex=0,
     main="IgG",xlab="SARS-CoV-2 RBD WT (Normalized MFI)",ylab="SARS-CoV-2 N (Normalized MFI)")
points(tmp[tmp$COVID.Status=="Infected","SARS2.RBD"],
       tmp[tmp$COVID.Status=="Infected","SARS2.N"],col=mycol[1],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"SARS2.RBD"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"SARS2.N"],col=mycol[2],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"SARS2.RBD"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"SARS2.N"],col=mycol[3],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Pre-Pandemic","SARS2.RBD"],
       tmp[tmp$COVID.Status=="Pre-Pandemic","SARS2.N"],col=mycol[4],pch=1,lwd=1.5,cex=1)
legend("topleft",c("Infected","Pre 2nd vaccination","Post 2nd vaccination","Pre-Pandemic"),fill=mycol)
#dev.off()

#b) - N vs RBD for MULTICOV-AB IgA highlighting different sample groups
#svg(paste("Output/Fig1/Fig1_b.svg",sep=""),7.5,5)
par(mfrow=c(1,1),mar=c(4,4,3,1),pty="m")
tmp <- M[M$Assay=="IgA",]
plot(tmp$SARS2.RBD,tmp$SARS2.N,log="xy",cex=0,
     main="IgA",xlab="SARS-CoV-2 RBD WT (Normalized MFI)",ylab="SARS-CoV-2 N (Normalized MFI)")
points(tmp[tmp$COVID.Status=="Infected","SARS2.RBD"],
       tmp[tmp$COVID.Status=="Infected","SARS2.N"],col=mycol[1],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"SARS2.RBD"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"SARS2.N"],col=mycol[2],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"SARS2.RBD"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"SARS2.N"],col=mycol[3],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Pre-Pandemic","SARS2.RBD"],
       tmp[tmp$COVID.Status=="Pre-Pandemic","SARS2.N"],col=mycol[4],pch=1,lwd=1.5,cex=1)
legend("topleft",c("Infected","Pre 2nd vaccination","Post 2nd vaccination","Pre-Pandemic"),fill=mycol)
#dev.off()

#c) Vaccinated Sample Time course of RBD Signals - IgG
#svg(paste("Output/Fig1/Fig1_c.svg",sep=""),7.5,5)
par(mfrow=c(1,1),mar=c(4,4,3,1),pty="m")
tmp <- M[M$Assay=="IgG",]
plot(tmp$dT1,tmp$SARS2.RBD,cex=0,log="y",
     main="IgG",ylab="SARS-CoV-2 RBD WT (Normalized MFI)",xlab="Days after first vaccination")
for (i in unique(tmp$DonorID)){
  lines(tmp[tmp$DonorID==i,"dT1"],tmp[tmp$DonorID==i,"SARS2.RBD"],col=mycol[4])
}
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"dT1"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"SARS2.RBD"],col=mycol[2],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"dT1"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"SARS2.RBD"],col=mycol[3],pch=1,lwd=1.5,cex=1)
legend("topleft",c("Pre 2nd vaccination","Post 2nd vaccination"),fill=mycol[2:3])
#dev.off()


#d) Vaccinated Sample Time course of RBD Signals - IgA
#svg(paste("Output/Fig1/Fig1_d.svg",sep=""),7.5,5)
par(mfrow=c(1,1),mar=c(4,4,3,1),pty="m")
tmp <- M[M$Assay=="IgA",]
plot(tmp$dT1,tmp$SARS2.RBD,cex=0,log="y",
     main="IgA",ylab="SARS-CoV-2 RBD WT (Normalized MFI)",xlab="Days after first vaccination")
for (i in unique(tmp$DonorID)){
  lines(tmp[tmp$DonorID==i,"dT1"],tmp[tmp$DonorID==i,"SARS2.RBD"],col=mycol[4])
}
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"dT1"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"SARS2.RBD"],col=mycol[2],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"dT1"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"SARS2.RBD"],col=mycol[3],pch=1,lwd=1.5,cex=1)
legend("topleft",c("Pre 2nd vaccination","Post 2nd vaccination"),fill=mycol[2:3])
#dev.off()
#####
#####Manuscript Figure 2####

mycol <- c("#ff4000","#00bfff","#004d67","#7F7F7F")

#svg(paste("Output/Fig2/Fig2_ab.svg",sep=""),8,5)
par(mfrow=c(1,2),mar=c(4,4,3,1),pty="m")
#highlight sample type outliers
pwpch <- as.character(S$sample_type)
pwpch[pwpch%in%c("Infected","negative","vaccinated")] <- 21
pwpch[pwpch%in%c("Inf_and_Vac","vac_other")] <- 24
pwpch <- as.numeric(pwpch)
#a)
boxplot(S[S$sample_type%in%c("Infected","Inf_and_Vac"),"IgA_RBD"],
        S[S$sample_type=="negative","IgA_RBD"],
        S[S$sample_type%in%c("vaccinated","vac_other"),"IgA_RBD"],
        outline=F,lwd=0.75,col=paste(mycol[c(1,4,2)],"30",sep=""),xlim=c(0.5,3.5),xaxt="n",yaxt="n",border="#00000000",
        main="SARS-CoV-2 WT RBD - IgA",ylab="Normalized MFI")
axis(2,at=c(1,5,10),labels=c(1,5,10))
beeswarm(list(S[S$sample_type%in%c("Infected","Inf_and_Vac"),"IgA_RBD"],
              S[S$sample_type=="negative","IgA_RBD"],
              S[S$sample_type%in%c("vaccinated","vac_other"),"IgA_RBD"]),
         col="black",bg=mycol[c(1,4,2)],pwpch=pwpch,xaxt="n",add=T,cex=1)
axis(side=1,at=c(1,2,3),labels=c("Infected","Negative","Vaccinated"))
boxplot(S[S$sample_type%in%c("Infected","Inf_and_Vac"),"IgA_RBD"],
        S[S$sample_type=="negative","IgA_RBD"],
        S[S$sample_type%in%c("vaccinated","vac_other"),"IgA_RBD"],
        outline=F,lwd=0.75,col="#00000000",add=T,yaxt="n",xaxt="n")
#b)
boxplot(log10(S[S$sample_type%in%c("Infected","Inf_and_Vac"),"IgG_RBD"]),
        log10(S[S$sample_type=="negative","IgG_RBD"]),
        log10(S[S$sample_type%in%c("vaccinated","vac_other"),"IgG_RBD"]),
        outline=F,lwd=0.75,col=paste(mycol[c(1,4,2)],"30",sep=""),xlim=c(0.5,3.5),xaxt="n",yaxt="n",border="#00000000",
        main="SARS-CoV-2 WT RBD - IgG",ylab="Normalized MFI")
axis(2,at=log10(c(0.1,1,10,100)),labels=c(0.1,1,10,100))
beeswarm(list(log10(S[S$sample_type%in%c("Infected","Inf_and_Vac"),"IgG_RBD"]),
              log10(S[S$sample_type=="negative","IgG_RBD"]),
              log10(S[S$sample_type%in%c("vaccinated","vac_other"),"IgG_RBD"])),
         col="black",bg=mycol[c(1,4,2)],pwpch=pwpch,xaxt="n",add=T,cex=1)
axis(side=1,at=c(1,2,3),labels=c("Infected","Negative","Vaccinated"))
boxplot(log10(S[S$sample_type%in%c("Infected","Inf_and_Vac"),"IgG_RBD"]),
        log10(S[S$sample_type=="negative","IgG_RBD"]),
        log10(S[S$sample_type%in%c("vaccinated","vac_other"),"IgG_RBD"]),
        outline=F,lwd=0.75,col="#00000000",add=T,yaxt="n",xaxt="n")
#dev.off()


#MWU Tests for boxplots
MWU <- data.frame(matrix(0,nrow=2,ncol=3))
names(MWU) <- c("InfVsNeg","VacvsNeg","InfVsVac")
rownames(MWU) <- c("IgG","IgA")
MWU["IgA","InfVsNeg"] <- wilcox.test(S[S$sample_type%in%c("Infected","Inf_and_Vac"),"IgA_RBD"],
                                     S[S$sample_type=="negative","IgA_RBD"],
                                     alternative="two.sided")$p.value
MWU["IgA","VacvsNeg"] <- wilcox.test(S[S$sample_type%in%c("vaccinated","vac_other"),"IgA_RBD"],
                                     S[S$sample_type=="negative","IgA_RBD"],
                                     alternative="two.sided")$p.value
MWU["IgA","InfVsVac"] <- wilcox.test(S[S$sample_type%in%c("vaccinated","vac_other"),"IgA_RBD"],
                                     S[S$sample_type%in%c("Infected","Inf_and_Vac"),"IgA_RBD"],
                                     alternative="two.sided")$p.value
MWU["IgG","InfVsNeg"] <- wilcox.test(S[S$sample_type%in%c("Infected","Inf_and_Vac"),"IgG_RBD"],
                                     S[S$sample_type=="negative","IgG_RBD"],
                                     alternative="two.sided")$p.value
MWU["IgG","VacvsNeg"] <- wilcox.test(S[S$sample_type%in%c("vaccinated","vac_other"),"IgG_RBD"],
                                     S[S$sample_type=="negative","IgG_RBD"],
                                     alternative="two.sided")$p.value
MWU["IgG","InfVsVac"] <- wilcox.test(S[S$sample_type%in%c("vaccinated","vac_other"),"IgG_RBD"],
                                     S[S$sample_type%in%c("Infected","Inf_and_Vac"),"IgG_RBD"],
                                     alternative="two.sided")$p.value
#write.table(MWU,"Output/Fig2/MWU.csv",sep=",",row.names=T,col.names=NA)

#####
#####Manuscript Figure 3####

mycol <- c("#ff4000","#00bfff","#004d67")

#a) IgG Signal RBD WT vs RBD UK
#svg(paste("Output/Fig3/Fig3_a.svg",sep=""),5,5)
par(mfrow=c(1,1),mar=c(4,4,3,1),pty="s")
tmp <- MFI[MFI$Assay=="IgG"&M$COVID.Status%in%c("Vaccinated","Infected"),]
plot(tmp$Anteo.RBD.WT,tmp$Anteo.RBD.UK,cex=0,xlim=range(tmp[,"Anteo.RBD.WT"]),ylim=range(tmp[,"Anteo.RBD.WT"]),
     main="IgG WT vs UK",xlab="SARS-CoV-2 RBD WT (MFI)",ylab="SARS-CoV-2 RBD UK (MFI)")
points(tmp[tmp$COVID.Status=="Infected","Anteo.RBD.WT"],
       tmp[tmp$COVID.Status=="Infected","Anteo.RBD.UK"],col=mycol[1],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"Anteo.RBD.WT"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"Anteo.RBD.UK"],col=mycol[2],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"Anteo.RBD.WT"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"Anteo.RBD.UK"],col=mycol[3],pch=1,lwd=1.5,cex=1)
legend("topleft",c("Infected","Pre 2nd vaccination","Post 2nd vaccination"),fill=mycol)
abline(a=0,b=1,lty=2,col="gray30")
legend("bottomright",paste("Kendall's tau\n",round(cor(tmp[,"Anteo.RBD.WT"],tmp[,"Anteo.RBD.UK"],method="kendall"),3)),bty="n")
#dev.off()

#b) IgG Signal RBD WT vs RBD ZAF
#svg(paste("Output/Fig3/Fig3_b.svg",sep=""),5,5)
par(mfrow=c(1,1),mar=c(4,4,3,1),pty="s")
tmp <- MFI[MFI$Assay=="IgG"&MFI$COVID.Status%in%c("Vaccinated","Infected"),]
plot(tmp$Anteo.RBD.WT,tmp$Anteo.RBD.ZAF,cex=0,xlim=range(tmp[,"Anteo.RBD.WT"]),ylim=range(tmp[,"Anteo.RBD.WT"]),
     main="IgG WT vs ZAF",xlab="SARS-CoV-2 RBD WT (MFI)",ylab="SARS-CoV-2 RBD ZAF (MFI)")
points(tmp[tmp$COVID.Status=="Infected","Anteo.RBD.WT"],
       tmp[tmp$COVID.Status=="Infected","Anteo.RBD.ZAF"],col=mycol[1],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"Anteo.RBD.WT"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"Anteo.RBD.ZAF"],col=mycol[2],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"Anteo.RBD.WT"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"Anteo.RBD.ZAF"],col=mycol[3],pch=1,lwd=1.5,cex=1)
legend("topleft",c("Infected","Pre 2nd vaccination","Post 2nd vaccination"),fill=mycol)
abline(a=0,b=1,lty=2,col="gray30")
legend("bottomright",paste("Kendall's tau\n",round(cor(tmp[,"Anteo.RBD.WT"],tmp[,"Anteo.RBD.ZAF"],method="kendall"),3)),bty="n")
#dev.off()



#####
#####Manuscript Figure 4####

mycol <- c("#ff4000","#00bfff","#004d67","#7F7F7F")


#a) IVNT IC50 correlation

pwcol <- rep(1,nrow(VNT))
pwcol[which(VNT$COVID.Status=="Vaccinated"&VNT$dT2<=0)] <- 2
pwcol[which(VNT$COVID.Status=="Vaccinated"&VNT$dT2>0)] <- 3
pwcol[which(VNT$COVID.Status=="Pre-Pandemic")] <- 4
pwcol1 <- pwcol[VNT$COVID.Status=="Vaccinated"]
pwcol2 <- pwcol[VNT$COVID.Status%in%c("Infected","Pre-Pandemic")]
#svg(paste("Output/Fig4/Fig4_a.svg",sep=""),7.5,5)
par(mfrow=c(1,1),mar=c(4,4,3,1),pty="m")
plot(1,1,cex=0,ylim=range(VNT$VNT50_WT),xlim=c(0.5,6.25),xaxt="n",log="y",xlab="",ylab=expression("VNT"["50"]),yaxt="n")
beeswarm(list(VNT[VNT$COVID.Status=="Vaccinated","VNT50_WT"],
              VNT[VNT$COVID.Status=="Vaccinated","VNT50_SA"]),
         col="black",pwbg=mycol[rep(pwcol1,2)],pch=21,xaxt="n",cex=2,at=c(1,2.5),add=T)
beeswarm(list(VNT[VNT$COVID.Status%in%c("Infected","Pre-Pandemic"),"VNT50_WT"],
              VNT[VNT$COVID.Status%in%c("Infected","Pre-Pandemic"),"VNT50_SA"]),
         col="black",pwbg=mycol[rep(pwcol2,2)],pch=21,xaxt="n",cex=2,at=c(4.25,5.75),add=T)
for (i in which(VNT$COVID.Status=="Vaccinated")){
  lines(c(1,2.5),c(VNT[i,"VNT50_WT"],VNT[i,"VNT50_SA"]),col="gray30")
}
for (i in which(VNT$COVID.Status%in%c("Infected","Pre-Pandemic"))){
  lines(c(4.25,5.75),c(VNT[i,"VNT50_WT"],VNT[i,"VNT50_SA"]),col="gray30")
}
axis(1,at=c(1,2.5,4.25,5.75),labels=rep(c("wt","RBD"),2))
mtext(c("Vaccinated","Infected/\nPre-Pandemic"),1,line=3,at=c(1.75,5))
axis(2,at=c(40,100,400,1000,4000),labels=c("40","100","400","1,000","4,000"))
#dev.off()

#Note that the lines connecting corresponding dots have to be connected corretly in Inkscape by moving ends !horizontally! only




#b) ACE2 competition - RBD WT vs ZAF

#svg(paste("Output/Fig4/Fig4_b.svg",sep=""),5,5)
par(mfrow=c(1,1),mar=c(4,4,3,1),pty="s")
plot(1-A$Anteo.RBD.WT,1-A$Anteo.RBD.ZAF,cex=0,ylim=c(-0.15,1),xlim=c(-0.15,1),
     main="IgG ACE2",ylab="RBD ZAF",xlab="RBD WT")
points(1-A[A$COVID.Status=="Infected","Anteo.RBD.WT"],
       1-A[A$COVID.Status=="Infected","Anteo.RBD.ZAF"],col=mycol[1],pch=1,lwd=1.5,cex=1)
points(1-A[A$COVID.Status=="Vaccinated"&A$dT2<=0,"Anteo.RBD.WT"],
       1-A[A$COVID.Status=="Vaccinated"&A$dT2<=0,"Anteo.RBD.ZAF"],col=mycol[2],pch=1,lwd=1.5,cex=1)
points(1-A[A$COVID.Status=="Vaccinated"&A$dT2>0,"Anteo.RBD.WT"],
       1-A[A$COVID.Status=="Vaccinated"&A$dT2>0,"Anteo.RBD.ZAF"],col=mycol[3],pch=1,lwd=1.5,cex=1)
abline(a=0,b=1,lty=2,col="gray30")
i <- cbind(1-A$Anteo.RBD.WT,1-A$Anteo.RBD.ZAF)
i <- lm(i[,2]~i[,1])
abline(i)
paste("y = ",round(i$coefficients[1],3)," + ",round(i$coefficients[2],3),"x",sep="") #This is put in the figure legend
legend("bottomright",
       paste("R2 = ",round(summary(i)$r.squared,3),sep=""),
       bty="n")
legend("topleft",c("Infected","Pre 2nd vaccination","Post 2nd vaccination"),fill=mycol)
#dev.off()


#c) ACE2 competition - RBD WT vs dT1

#svg(paste("Output/Fig4/Fig4_c.svg",sep=""),7.5,5)
par(mfrow=c(1,1),mar=c(4,4,3,1),pty="m")
plot(A$dT1,1-A$Anteo.RBD.WT,cex=0,ylim=c(-0.15,1),
     main="IgG ACE2",ylab="RBD WT",xlab="Time after first vaccination")
points(A[A$COVID.Status=="Vaccinated"&A$dT2<=0,"dT1"],
       1-A[A$COVID.Status=="Vaccinated"&A$dT2<=0,"Anteo.RBD.WT"],col=mycol[2],pch=1,lwd=1.5,cex=1)
points(A[A$COVID.Status=="Vaccinated"&A$dT2>0,"dT1"],
       1-A[A$COVID.Status=="Vaccinated"&A$dT2>0,"Anteo.RBD.WT"],col=mycol[3],pch=1,lwd=1.5,cex=1)
for (i in unique(tmp$DonorID)){
  lines(A[A$DonorID==i,"dT1"],1-A[A$DonorID==i,"Anteo.RBD.WT"],col="gray50")
}
legend("topleft",c("Pre 2nd vaccination","Post 2nd vaccination"),fill=mycol[2:3])
#dev.off()



#d) ACE2 competition - RBD ZAF vs dT1

#svg(paste("Output/Fig4/Fig4_d.svg",sep=""),7.5,5)
par(mfrow=c(1,1),mar=c(4,4,3,1),pty="m")
plot(A$dT1,1-A$Anteo.RBD.ZAF,cex=0,ylim=c(-0.15,1),
     main="IgG ACE2",ylab="RBD ZAF",xlab="Time after first vaccination")
points(A[A$COVID.Status=="Vaccinated"&A$dT2<=0,"dT1"],
       1-A[A$COVID.Status=="Vaccinated"&A$dT2<=0,"Anteo.RBD.ZAF"],col=mycol[2],pch=1,lwd=1.5,cex=1)
points(A[A$COVID.Status=="Vaccinated"&A$dT2>0,"dT1"],
       1-A[A$COVID.Status=="Vaccinated"&A$dT2>0,"Anteo.RBD.ZAF"],col=mycol[3],pch=1,lwd=1.5,cex=1)
for (i in unique(tmp$DonorID)){
  lines(A[A$DonorID==i,"dT1"],1-A[A$DonorID==i,"Anteo.RBD.ZAF"],col="gray50")
}
legend("topleft",c("Pre 2nd vaccination","Post 2nd vaccination"),fill=mycol[2:3])
#dev.off()






#####
#####Supplementary Figure 1 - Vaccinated sera results in a plateau effect in MULTICOV-AB######

#Read in of raw data
DS <- read.csv("input/DilSeries.csv",h=T,sep=",",stringsAsFactors=F)

mycol <- c(brewer.pal(10,"Paired"),"gray20")

#Plot
#svg(paste("Output/SupplFig1/SupplFig1_raw.svg",sep=""),7.5,5)
par(mfrow=c(1,1),mar=c(4,4,3,1),pty="m")
plot(1,1,cex=0,log="xy",xlim=c(350,32000),ylim=c(30,45000),xaxt="n",yaxt="n",
     xlab="Sample dilution factor",ylab="SARS-CoV-2 RBD (MFI)")
axis(1,at=c(400,800,1600,3200,6400,12800),labels=c("1:400","1:800","1:1,600","1:3,200","1:6,400","1:12,800"))
axis(2,at=c(40,400,4000,40000),labels=c("40","400","4,000","40,000"))
for (i in 1:nrow(DS)){
  lines(c(400,800,1600,3200,6400,12800),t(DS[i,4:9]),col=mycol[i])
}
for (i in 1:nrow(DS)){
  points(c(400,800,1600,3200,6400,12800),t(DS[i,4:9]),col=mycol[i],pch=1,lwd=1.5,cex=1)
}
legend(14000,40000,paste(DS$Donor.ID),fill=mycol,bty="n")
legend(23500,40000,paste("d",DS[1:10,"dT..1st.shot."],sep=""),bty="n")
#dev.off()



#####
#####Supplementary Figure 2 - Saliva IgG results are verified by a second ELISA######

mycol <- c("#ff4000","#00bfff","#004d67","#7F7F7F")
#Read-in of data from Saliva IgG ELISA
S2 <- read.csv("Input/ELISA_Saliva_Results.csv",h=T)

#Change all values at 0 to 0.1 to enable log scale display
S2$Saliva_IgG_adj <- S2$Saliva_IgG
S2[S2$Saliva_IgG_adj==0,"Saliva_IgG_adj"] <- 0.1

#Plot
#svg(paste("Output/SupplFig2/SupplFig2.svg",sep=""),8,5)
par(mfrow=c(1,2),mar=c(4,4,3,1))
#highlight sample type outliers
pwpch <- as.character(S2$Sample_Type)
pwpch[pwpch%in%c("Infected","Negative","Vaccinated")] <- 21
pwpch[pwpch%in%c("Inf_and_vac","Vaccinated_other")] <- 24
pwpch <- as.numeric(pwpch)
#a)
boxplot(log10(S2[S2$Sample_Type%in%c("Infected","Inf_and_vac"),"Saliva_IgG_adj"]),
        log10(S2[S2$Sample_Type=="Negative","Saliva_IgG_adj"]),
        log10(S2[S2$Sample_Type%in%c("Vaccinated","Vaccinated_other"),"Saliva_IgG_adj"]),
        outline=F,lwd=0.75,col=paste(mycol[c(1,4,2)],"30",sep=""),
        xlim=c(0.5,3.5),xaxt="n",yaxt="n",border="#00000000",
        main="Saliva IgG ELISA",ylab="Saliva IgG Signal")
axis(2,at=log10(c(0.1,1,10,100,1000)),labels=c("0.1","1","10","100","1,000"))
beeswarm(list(log10(S2[S2$Sample_Type%in%c("Infected","Inf_and_vac"),"Saliva_IgG_adj"]),
              log10(S2[S2$Sample_Type=="Negative","Saliva_IgG_adj"]),
              log10(S2[S2$Sample_Type%in%c("Vaccinated","Vaccinated_other"),"Saliva_IgG_adj"])),
         col="black",bg=mycol[c(1,4,2)],pwpch=pwpch,xaxt="n",add=T,cex=1)
axis(side=1,at=c(1,2,3),labels=c("Infected","Negative","Vaccinated"))
boxplot(log10(S2[S2$Sample_Type%in%c("Infected","Inf_and_vac"),"Saliva_IgG_adj"]),
        log10(S2[S2$Sample_Type=="Negative","Saliva_IgG_adj"]),
        log10(S2[S2$Sample_Type%in%c("Vaccinated","Vaccinated_other"),"Saliva_IgG_adj"]),
        outline=F,lwd=0.75,col="#00000000",add=T,yaxt="n",xaxt="n")
#dev.off()

#Calculate MWU test
MWU_ELISA <- data.frame(matrix(0,nrow=1,ncol=3))
names(MWU_ELISA) <- c("InfVsNeg","VacvsNeg","InfVsVac")
rownames(MWU_ELISA) <- c("IgG")
MWU_ELISA["IgG","InfVsNeg"] <- wilcox.test(S2[S2$Sample_Type%in%c("Infected","Inf_and_vac"),"Saliva_IgG"],
                                     S2[S2$Sample_Type=="Negative","Saliva_IgG"],
                                     alternative="two.sided")$p.value
MWU_ELISA["IgG","VacvsNeg"] <- wilcox.test(S2[S2$Sample_Type%in%c("Vaccinated","Vaccinated_other"),"Saliva_IgG"],
                                     S2[S2$Sample_Type=="Negative","Saliva_IgG"],
                                     alternative="two.sided")$p.value
MWU_ELISA["IgG","InfVsVac"] <- wilcox.test(S2[S2$Sample_Type%in%c("Vaccinated","Vaccinated_other"),"Saliva_IgG"],
                                     S2[S2$Sample_Type%in%c("Infected","Inf_and_vac"),"Saliva_IgG"],
                                     alternative="two.sided")$p.value
#write.table(MWU_ELISA,"Output/SupplFig2/MWU_ELISA.csv",sep=",",row.names=T,col.names=NA)



#####
#####Supplementary Figure 3 - Pfizer BNT-162b does not offer any cross-protection against endemic coronaviruses######

mycol <- c("#ff4000","#00bfff","#004d67","#7F7F7F")

#a) 229E
#svg(paste("Output/SupplFig3/SupplFig3_raw.svg",sep=""),9,6)
par(mfrow=c(2,2),mar=c(4,4,3,1),pty="m")
tmp <- M[M$Assay=="IgG",]
plot(tmp[tmp$COVID.Status=="Vaccinated","dT1"],tmp[tmp$COVID.Status=="Vaccinated","X229E.S1"],cex=0,log="y",
     main="hCoV 229E",ylab="hCoV 229E S1 (Normalized MFI)",xlab="Time post first vaccination (days)")
for (i in unique(tmp$DonorID)){
  lines(tmp[tmp$DonorID==i,"dT1"],tmp[tmp$DonorID==i,"X229E.S1"],col=mycol[4])
}
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"dT1"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"X229E.S1"],col=mycol[2],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"dT1"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"X229E.S1"],col=mycol[3],pch=1,lwd=1.5,cex=1)
legend("topleft",c("Pre 2nd vaccination","Post 2nd vaccination"),fill=mycol[2:3])

#b) NL63
tmp <- M[M$Assay=="IgG",]
plot(tmp[tmp$COVID.Status=="Vaccinated","dT1"],tmp[tmp$COVID.Status=="Vaccinated","NL63.S1"],cex=0,log="y",
     main="hCoV NL63",ylab="hCoV NL63 S1 (Normalized MFI)",xlab="Time post first vaccination (days)")
for (i in unique(tmp$DonorID)){
  lines(tmp[tmp$DonorID==i,"dT1"],tmp[tmp$DonorID==i,"NL63.S1"],col=mycol[4])
}
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"dT1"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"NL63.S1"],col=mycol[2],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"dT1"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"NL63.S1"],col=mycol[3],pch=1,lwd=1.5,cex=1)
legend("topleft",c("Pre 2nd vaccination","Post 2nd vaccination"),fill=mycol[2:3])

#c) OC43
tmp <- M[M$Assay=="IgG",]
plot(tmp[tmp$COVID.Status=="Vaccinated","dT1"],tmp[tmp$COVID.Status=="Vaccinated","OC43.S1"],cex=0,log="y",
     main="hCoV OC43",ylab="hCoV OC43 S1 (Normalized MFI)",xlab="Time post first vaccination (days)")
for (i in unique(tmp$DonorID)){
  lines(tmp[tmp$DonorID==i,"dT1"],tmp[tmp$DonorID==i,"OC43.S1"],col=mycol[4])
}
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"dT1"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"OC43.S1"],col=mycol[2],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"dT1"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"OC43.S1"],col=mycol[3],pch=1,lwd=1.5,cex=1)
legend("topleft",c("Pre 2nd vaccination","Post 2nd vaccination"),fill=mycol[2:3])

#d) HKU1
tmp <- M[M$Assay=="IgG",]
plot(tmp[tmp$COVID.Status=="Vaccinated","dT1"],tmp[tmp$COVID.Status=="Vaccinated","HKU1.S1"],cex=0,log="y",
     main="hCoV HKU1",ylab="hCoV HKU1 S1 (Normalized MFI)",xlab="Time post first vaccination (days)")
for (i in unique(tmp$DonorID)){
  lines(tmp[tmp$DonorID==i,"dT1"],tmp[tmp$DonorID==i,"HKU1.S1"],col=mycol[4])
}
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"dT1"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"HKU1.S1"],col=mycol[2],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"dT1"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"HKU1.S1"],col=mycol[3],pch=1,lwd=1.5,cex=1)
legend("topleft",c("Pre 2nd vaccination","Post 2nd vaccination"),fill=mycol[2:3])
#dev.off()


#####
#####Supplementary Figure 4 - RBD WT vs Mink and LA Mutants, IgG Response#####

mycol <- c("#ff4000","#00bfff","#004d67")

#a) IgG Signal RBD WT vs RBD Mink
#svg(paste("Output/SupplFig4/SupplFig4_a.svg",sep=""),5,5)
par(mfrow=c(1,1),mar=c(4,4,3,1),pty="s")
tmp <- MFI[MFI$Assay=="IgG"&M$COVID.Status%in%c("Vaccinated","Infected"),]
plot(tmp$Anteo.RBD.WT,tmp$Anteo.RBD.Mink,cex=0,xlim=range(tmp[,"Anteo.RBD.Mink"]),ylim=range(tmp[,"Anteo.RBD.Mink"]),
     main="IgG WT vs Mink",xlab="SARS-CoV-2 RBD WT (MFI)",ylab="SARS-CoV-2 RBD Mink (MFI)")
points(tmp[tmp$COVID.Status=="Infected","Anteo.RBD.WT"],
       tmp[tmp$COVID.Status=="Infected","Anteo.RBD.Mink"],col=mycol[1],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"Anteo.RBD.WT"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"Anteo.RBD.Mink"],col=mycol[2],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"Anteo.RBD.WT"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"Anteo.RBD.Mink"],col=mycol[3],pch=1,lwd=1.5,cex=1)
legend("topleft",c("Infected","Pre 2nd vaccination","Post 2nd vaccination"),fill=mycol)
abline(a=0,b=1,lty=2,col="gray30")
legend("bottomright",paste("Kendall's tau\n",round(cor(tmp[,"Anteo.RBD.WT"],tmp[,"Anteo.RBD.Mink"],method="kendall"),3)),bty="n")
#dev.off()

#b) IgG Signal RBD WT vs RBD LA
#svg(paste("Output/SupplFig4/SupplFig4_b.svg",sep=""),5,5)
par(mfrow=c(1,1),mar=c(4,4,3,1),pty="s")
tmp <- MFI[MFI$Assay=="IgG"&M$COVID.Status%in%c("Vaccinated","Infected"),]
plot(tmp$Anteo.RBD.WT,tmp$Anteo.RBD.LA,cex=0,xlim=range(tmp[,"Anteo.RBD.WT"]),ylim=range(tmp[,"Anteo.RBD.WT"]),
     main="IgG WT vs LA",xlab="SARS-CoV-2 RBD WT (MFI)",ylab="SARS-CoV-2 RBD LA (MFI)")
points(tmp[tmp$COVID.Status=="Infected","Anteo.RBD.WT"],
       tmp[tmp$COVID.Status=="Infected","Anteo.RBD.LA"],col=mycol[1],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"Anteo.RBD.WT"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2<=0,"Anteo.RBD.LA"],col=mycol[2],pch=1,lwd=1.5,cex=1)
points(tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"Anteo.RBD.WT"],
       tmp[tmp$COVID.Status=="Vaccinated"&tmp$dT2>0,"Anteo.RBD.LA"],col=mycol[3],pch=1,lwd=1.5,cex=1)
legend("topleft",c("Infected","Pre 2nd vaccination","Post 2nd vaccination"),fill=mycol)
abline(a=0,b=1,lty=2,col="gray30")
legend("bottomright",paste("Kendall's tau\n",round(cor(tmp[,"Anteo.RBD.WT"],tmp[,"Anteo.RBD.LA"],method="kendall"),3)),bty="n")
#dev.off()



#####
#####Supplementary Figure 5 - ACE2 and NBP####

mycol <- c("#ff4000","#00bfff","#004d67")

#NBP can only measure positive samples,so for ACE2 retain only samples in tmpN
tmpA <- A[A$SampleID%in%unique(N$SampleID),]


#a) WT reactivity
#svg(paste("Output/SupplFig5/SupplFig5_raw.svg",sep=""),5,5)
par(mfrow=c(1,1),mar=c(4,4,3,1),pty="s")
plot(1-tmpA$Anteo.RBD.WT,N$Anteo.RBD.WT,cex=0,xlim=c(-0.1,1),ylim=c(-0.1,1),
     main="IgG RBD WT",ylab="NeutroBodyPlex",xlab="ACE2 binding assay")
points(1-tmpA[tmpA$COVID.Status=="Infected","Anteo.RBD.WT"],
       N[N$COVID.Status=="Infected","Anteo.RBD.WT"],col=mycol[1],pch=1,lwd=1.5)
points(1-tmpA[tmpA$COVID.Status=="Vaccinated"&tmpA$dT2<=0,"Anteo.RBD.WT"],
       N[N$COVID.Status=="Vaccinated"&N$dT2<=0,"Anteo.RBD.WT"],col=mycol[2],pch=1,lwd=1.5)
points(1-tmpA[tmpA$COVID.Status=="Vaccinated"&tmpA$dT2>0,"Anteo.RBD.WT"],
       N[N$COVID.Status=="Vaccinated"&N$dT2>0,"Anteo.RBD.WT"],col=mycol[3],pch=1,lwd=1.5)
i <- cbind(1-tmpA$Anteo.RBD.WT,N$Anteo.RBD.WT)
i <- lm(i[,2]~i[,1])
abline(i)
paste("y = ",round(i$coefficients[1],3)," + ",round(i$coefficients[2],3),sep="")
legend("bottomright",
       paste("R2 = ",round(summary(i)$r.squared,3),sep=""),
       bty="n")
legend("topleft",c("Infected","Pre 2nd vaccination","Post 2nd vaccination"),fill=mycol)
#dev.off()

#####
#####Supplementary Figure 6 - ACE2 COMPETITION ACROSS ALL MUTANTS####

#a) ACE2 competition - RBD WT vs UK

#svg(paste("Output/SupplFig6/SupplFig6_a.svg",sep=""),5,5)
par(mfrow=c(1,1),mar=c(4,4,3,1),pty="s")
plot(1-A$Anteo.RBD.WT,1-A$Anteo.RBD.UK,cex=0,ylim=c(-0.15,1),xlim=c(-0.15,1),
     main="IgG ACE2",ylab="RBD UK",xlab="RBD WT")
points(1-A[A$COVID.Status=="Infected","Anteo.RBD.WT"],
       1-A[A$COVID.Status=="Infected","Anteo.RBD.UK"],col=mycol[1],pch=1,lwd=1.5,cex=1)
points(1-A[A$COVID.Status=="Vaccinated"&A$dT2<=0,"Anteo.RBD.WT"],
       1-A[A$COVID.Status=="Vaccinated"&A$dT2<=0,"Anteo.RBD.UK"],col=mycol[2],pch=1,lwd=1.5,cex=1)
points(1-A[A$COVID.Status=="Vaccinated"&A$dT2>0,"Anteo.RBD.WT"],
       1-A[A$COVID.Status=="Vaccinated"&A$dT2>0,"Anteo.RBD.UK"],col=mycol[3],pch=1,lwd=1.5,cex=1)
abline(a=0,b=1,lty=2,col="gray30")
i <- cbind(1-A$Anteo.RBD.WT,1-A$Anteo.RBD.UK)
i <- lm(i[,2]~i[,1])
abline(i)
paste("y = ",round(i$coefficients[1],3)," + ",round(i$coefficients[2],3),sep="")
legend("bottomright",
       paste("R2 = ",round(summary(i)$r.squared,3),sep=""),
       bty="n")
legend("topleft",c("Infected","Pre 2nd vaccination","Post 2nd vaccination"),fill=mycol)
#dev.off()



#b) ACE2 competition - RBD WT vs Mink

#svg(paste("Output/SupplFig6/SupplFig6_b.svg",sep=""),5,5)
par(mfrow=c(1,1),mar=c(4,4,3,1),pty="s")
plot(1-A$Anteo.RBD.WT,1-A$Anteo.RBD.Mink,cex=0,ylim=c(-0.15,1),xlim=c(-0.15,1),
     main="IgG ACE2",ylab="RBD Mink",xlab="RBD WT")
points(1-A[A$COVID.Status=="Infected","Anteo.RBD.WT"],
       1-A[A$COVID.Status=="Infected","Anteo.RBD.Mink"],col=mycol[1],pch=1,lwd=1.5,cex=1)
points(1-A[A$COVID.Status=="Vaccinated"&A$dT2<=0,"Anteo.RBD.WT"],
       1-A[A$COVID.Status=="Vaccinated"&A$dT2<=0,"Anteo.RBD.Mink"],col=mycol[2],pch=1,lwd=1.5,cex=1)
points(1-A[A$COVID.Status=="Vaccinated"&A$dT2>0,"Anteo.RBD.WT"],
       1-A[A$COVID.Status=="Vaccinated"&A$dT2>0,"Anteo.RBD.Mink"],col=mycol[3],pch=1,lwd=1.5,cex=1)
abline(a=0,b=1,lty=2,col="gray30")
i <- cbind(1-A$Anteo.RBD.WT,1-A$Anteo.RBD.Mink)
i <- lm(i[,2]~i[,1])
abline(i)
paste("y = ",round(i$coefficients[1],3)," + ",round(i$coefficients[2],3),sep="")
legend("bottomright",
       paste("R2 = ",round(summary(i)$r.squared,3),sep=""),
       bty="n")
legend("topleft",c("Infected","Pre 2nd vaccination","Post 2nd vaccination"),fill=mycol)
#dev.off()


#c) ACE2 competition - RBD WT vs LA

#svg(paste("Output/SupplFig6/SupplFig6_c.svg",sep=""),5,5)
par(mfrow=c(1,1),mar=c(4,4,3,1),pty="s")
plot(1-A$Anteo.RBD.WT,1-A$Anteo.RBD.LA,cex=0,ylim=c(-0.15,1),xlim=c(-0.15,1),
     main="IgG ACE2",ylab="RBD LA",xlab="RBD WT")
points(1-A[A$COVID.Status=="Infected","Anteo.RBD.WT"],
       1-A[A$COVID.Status=="Infected","Anteo.RBD.LA"],col=mycol[1],pch=1,lwd=1.5,cex=1)
points(1-A[A$COVID.Status=="Vaccinated"&A$dT2<=0,"Anteo.RBD.WT"],
       1-A[A$COVID.Status=="Vaccinated"&A$dT2<=0,"Anteo.RBD.LA"],col=mycol[2],pch=1,lwd=1.5,cex=1)
points(1-A[A$COVID.Status=="Vaccinated"&A$dT2>0,"Anteo.RBD.WT"],
       1-A[A$COVID.Status=="Vaccinated"&A$dT2>0,"Anteo.RBD.LA"],col=mycol[3],pch=1,lwd=1.5,cex=1)
abline(a=0,b=1,lty=2,col="gray30")
i <- cbind(1-A$Anteo.RBD.WT,1-A$Anteo.RBD.LA)
i <- lm(i[,2]~i[,1])
abline(i)
paste("y = ",round(i$coefficients[1],3)," + ",round(i$coefficients[2],3),sep="")
legend("bottomright",
       paste("R2 = ",round(summary(i)$r.squared,3),sep=""),
       bty="n")
legend("topleft",c("Infected","Pre 2nd vaccination","Post 2nd vaccination"),fill=mycol)
#dev.off()

#####