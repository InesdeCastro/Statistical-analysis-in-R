#Auditory_System_StatisticalAnalysis

library(readr)
Auditory_System <- read_delim("C:/Users/inesc/Bioinformatics/Genome comparison/Auditory_System.csv",  ";", escape_double = FALSE, trim_ws = TRUE)
View(Auditory_System)

dim(Auditory_System)
is.na(Auditory_System)
audi= na.exclude(Auditory_System)
View(audi)
dim(audi)

#******NUMBER OF FRAGMENTS: TESTING DIFFERENCES BETWEEN METHODS (PACBIO VS SOLEXA/ILLUMINA/SANGER)******


nr_frag_PB=audi$`Number of fragments_PB`
length(nr_frag_PB)

nr_frag_SIS=audi$`Number of fragments_SIS`
length(nr_frag_SIS)


##NORMALITY TEST------------------------------------------------------------------------------------
#PacBio values:

ks.test(nr_frag_PB, 'pnorm', mean(nr_frag_PB), sd(nr_frag_PB))
#p-value = 3.008e-09
shapiro.test(nr_frag_PB)
#p-value = 3.394e-12


#Solexa/Illumina/Sanger values:

ks.test(nr_frag_SIS,'pnorm', mean(nr_frag_SIS), sd(nr_frag_SIS))
#p-value = 0.001656

shapiro.test(nr_frag_SIS)
#p-value = 1.676e-09


#In both cases, data has not a normal distribution
#It is not possible to normalize the data, so non parametric tests will be performed



hist(nr_frag_PB) # --> Just to check the ditribution of data
hist(nr_frag_SIS) # --> Just to check the ditribution of data



## INDEPENDENT 2-GROUP TEST-------------------------------------------------------------------------
#H0: There are no significative differences in the number of fragments, between PacBio and Solexa/Illumina/Sanger
#H1: There are significative differences in the number of fragments, between PacBio and Solexa/Illumina/Sanger

wilcox.test(nr_frag_PB, nr_frag_SIS)
#p-value =  0.000504 -->  Reject H0  --> Significative differences. 

ks.test(nr_frag_PB, nr_frag_SIS)
#p-value = 0.03268 --> Reject H0  --> Significative differences.


#******INTEGRITY: TESTING DIFFERENCES BETWEEN METHODS (PACBIO VS SOLEXA/ILLUMINA/SANGER)******

int_PB= audi$`Integrity of gene_PB`
length(int_PB)


int_SIS=audi$`Integrity of gene_SIS`
length(int_SIS)

outlier(int_PB)
outlier(int_SIS)

##NORMALITY TEST------------------------------------------------------------------------------------

#PacBio values:
shapiro.test(int_PB)
#p-value =  1.229e-11

#Solexa/Illumina/Sanger values:
shapiro.test(int_SIS)
# p-value = 1.619e-06

#In both cases, data has not a normal distribution
#It is not possible to normalize the data, so non parametric tests will be performed


## INDEPENDENT 2-GROUP TEST-------------------------------------------------------------------------
#H0: There are no significative differences in the integrity of genes, between PacBio and Solexa/Illumina/Sanger
#H1: There are significative differences in the integrity of genes, between PacBio and Solexa/Illumina/Sanger


wilcox.test(int_PB, int_SIS)
#p-value = 0.007545

ks.test(int_PB, int_SIS)
#p-value = 0.03268

#There ARE significative statistical differences between the two methods!!!!!!!!

par(mfrow=c(1,2))
hist(int_PB, col=2, main=" Gene Integrity_PacBio")
hist(int_SIS, col=4,  main=" Gene Integrity_SIS")


#ARTEFACTS: TESTING DIFFERENCES BETWEEN METHODS (PACBIO VS SOLEXA/ILLUMINA/SANGER)----------------------------------

art_PB=audi$Artefacts_PB
art_SIS=audi$Artefacts_SIS


#outlier(art_PB) #nada de significativo
#outlier(art_SIS) #nada de significativo

##NORMALITY TEST------------------------------------------------------------------------------------

#PacBio values:
shapiro.test(art_PB)
#p-value = 7.688e-11

#Solexa/Illumina/Sanger values:
shapiro.test(art_SIS)
#p-value =  1.376e-10


#In both cases, data has not a normal distribution
#It is not possible to normalize the data, so non parametric tests will be performed


## INDEPENDENT 2-GROUP TEST-------------------------------------------------------------------------
#H0: There are no significative differences in the number of fragments founded, between PacBio and Solexa/Illumina/Sanger
#H1: There are significative differences in the number of fragments founded, between PacBio and Solexa/Illumina/Sanger


wilcox.test(art_PB, art_SIS)
#p-value = p-value =  0.9461  --> No significative differences.


ks.test(art_PB, art_SIS)
#p-value = 1   --> No significative differences.

########################################################TESTES DE GRAFICOS################

#par(mfrow=c(1,2))
#hist(int_PB,col=7, main="A. System Gene's Integrity", xlab= "Integrity % in PacBio")
#hist(int_SIS, col=10, main="A. System Gene's Integrity", xlab="Integrity % in Solexa/Illumina/Sanger")


#all_int= data.frame(int_PB, int_SIS)
#all_int

#par(mfrow=c(1,1))

#barplot(t(as.matrix(all_int)),beside = TRUE, main = "Gene's Integrity",legend.text = TRUE)
#hist(t(as.matrix(all_int)),beside = TRUE,plot = TRUE, main = "Gene's Integrity",legend.text = TRUE)


#pop é um dataframe com duas colunas
#names.arg=c()


######################################## BOXPLOTS SUGERIDOS ####################################################
#par(mfrow=c(1,1))

col2rgb("lightblue")
col2rgb(c("lightblue", "lightgreen", "pink", "orange", "red"))



c1_blue = rgb(173,216,230,max = 255, alpha = 100, names = "lt.blue")
c2_pink = rgb(255,192,203, max = 255, alpha = 110, names = "lt.pink")
c3_orange = rgb(255,165,0, max = 255, alpha = 80, names = "lt.orange")
c4_red = rgb(255,0,0, max = 255, alpha = 80, names = "lt.red")

par(mfrow=c(1,3))


#FRAGMENTS
boxplot(nr_frag_PB, nr_frag_SIS,
        main = "Number of Fragments",
        at = c(1,2),
        names = c("PacBio", "S./I./S."),
        ylab='Nr of Fragments',
        las = 2,
        col = c("lightblue4", "lightcoral"),
        border=c("lightblue4","lightcoral"),
        horizontal = FALSE,
        notch = FALSE
)


#INTEGRITY
boxplot(int_PB, int_SIS,
        main = "Gene Integrity",
        at = c(1,2),
        names = c("PacBio", "S./I./S."),
        ylab='Integrity %',
        las = 2,
        col = c("lightblue4", "lightcoral"),
        border=c("lightblue4","lightcoral"),
        horizontal = FALSE,
        notch = FALSE
)

mtext("Auditory System", side=3, line=3, adj= 1,padj=0, outer = F, cex = 1.0)

#ARTEFACTS
boxplot(art_PB, art_SIS,
        main = "Number of Artefacts",
        at = c(1,2),
        names = c("PacBio", "S./I./S."),
        ylab='Nr of Artefacts',
        las = 2,
        col = c("lightblue4", "lightcoral"),
        border=c("lightblue4","lightcoral"),
        horizontal = FALSE,
        notch = FALSE
)



####################HISTOGRAMS#######################################3

#FOR DIFFERENT BREAKPOINTSn
hgPB=hist(int_PB, plot=FALSE)
hgSIS=hist(int_SIS, plot=FALSE)
range(c(hgPB$breaks, hgSIS$breaks)) # Get range for x-axis
#0.55 1.05
max(c(hgPB$count, hgSIS$count))  # Get range for y-axis
#28
plot(hgPB, col = c1_blue, xlim = c(0.50, 1.15), ylim = c(0,30))
plot(hgSIS, add = TRUE, col = c2_pink)




plot(hgPB, col=c1_blue)

plot(hgSIS, col=c3_orange, add=TRUE, legend.text=T)
plot(hgSIS, col=c4_red, add=TRUE, legend.text=T)

plot(hgSIS, col=c2_pink, add=TRUE, legend.text=T)



######################################################################################################################
##Number of genes founded between Pacbio and SIS


auditory_pacbio=c(9,9,9,9)
# p-value =0.4533
auditory_sis=c(9,9,9,8)
# p-value =0.4533
shapiro.test(auditory_pacbio)
shapiro.test(auditory_sis)

wilcox.test(auditory_pacbio, auditory_sis)
#p-value = 0.4533
#No significative differences


