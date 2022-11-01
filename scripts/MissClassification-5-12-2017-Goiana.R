# require(magrittr)
# require(tidyr)
# require(dplyr)
# require(ggplot2)
# require(sp)
# require(leaflet)
# # require(simPop)
# require(plyr)
# require(survey)
# require(sampling)
# require(reshape)
# require(reshape2)
# require(grid)
# require(gridExtra)
# require(ncf)
# # require(pgirmess)
# require(spdep)
# require(data.table)



# Facilita la instalación de librerías
if(!require(install.load)){install.packages("install.load", dep=TRUE); require("install.load")}

# Librerías requeridas
install_load("magrittr", "tidyr", "dplyr", "ggplot2", "leaflet",
             "reshape", "reshape2", "grid", "gridExtra", "data.table",
             "latex2exp", "gamlss", "plot3D", "compiler", "moments",
             "tables", "fitdistrplus", "RColorBrewer", "lmomco", "lattice", 
             "seleniumPipes", "RSelenium", "readxl", "mailR", "mongolite", 
             "sp", "rgeos", "stringr", "latticeExtra", "plyr",  "reshapeGUI", 
             "gridBase", "gridExtra", "entropy", "EntropyEstimation", 
             "infotheo", "distrEx", "gamlss", "betareg", "readxl", "haven", 
             "survey", "sampling", "ncf", "spdep", "data.table")






#-------------------------------------------------
# Population parameter
# x.1  # Highly Cropland - Stratum I
# y.2  # Cropland - Stratum II
# z.3  # Non-Cropland - Stratum III

# Stratum  size
N.1 = 314 # Highly Cropland - Stratum I
N.2 = 232 # Cropland - Stratum II
N.3 = 376 # Non-Cropland - Stratum III

# N total
NT= N.1+N.2+N.3

# Here, The population by Stratum is Normally Distributed 
# y_k ~ N(mu.k, srt(v.k)), k=1,2,3


# Information based on Estimator 2 - Joao Eudes Miqueias - 11/10/2017
# Mean of Population by Stratum 
mu.1 = 324150850.90 / N.1   # Highly Cropland - Stratum I
mu.2 = 155530715.36 / N.2    # Cropland - Stratum II
mu.3 = 39838562 / N.3     # Non-Cropland - Stratum III

# Variance of Population by Stratum 
v.1 = 34449029.89^2 / N.1^2 
v.2 = 38739291.47^2 / N.2^2 
v.3 = 39838562^2 / N.3^2 # Estou colocando uma variância bem pequena próximo de zero pafa não dar problema, 
#-------------------------------------------------

#-------------------------------------------------
# Transition Proportions
# Let pij be the sample proportion of segments classified in stratum i 
# when using satellite imagery photo interpretation, and classified as 
# stratum j in the field. Based Table 2 of Report

# Transfering I for II and III
p.1.2 = 0.119
p.1.3 = 0.333 

# Transfering II for I and III
p.2.1 = 0.266
p.2.3 = 0.001 # (Truth 0.0)

# Tranfering III for I and II  
p.3.1 = 0.400 
p.3.2 = 0.667
#-------------------------------------------------

# Seed to Reproducible Research
set.seed(121)

#-------------------------------------------------
# We assumed an order relation between the means of populations.
# Here mu 1> mu 2> mu 3 

#----------------------------------------
# Population Model to Stratum I
# yk= 50 - 20I12 - 30I13(1-I12)  + ek 
I.1.2 = rbinom(N.1, 1, p.1.2)
I.1.3 = rbinom(N.1, 1, p.1.3)
e.1   = rnorm(N.1, 0, sqrt(v.1))

# Response Stratum I
y.k   = mu.1 - mu.2*I.1.2 - mu.3*I.1.3*(1-I.1.2) + e.1
ST.1  = data.frame(stratum=rep("I", N.1), y.k, I.1.2, I.1.3)

# Poststratification 
# PostStratum is a modification by misclassification
ST.1$PostStratum[ST.1$I.1.2+2*ST.1$I.1.3 == 0] <- "I"
ST.1$PostStratum[ST.1$I.1.2+2*ST.1$I.1.3 == 1] <- "II"
ST.1$PostStratum[ST.1$I.1.2+2*ST.1$I.1.3 == 2] <- "III"
ST.1$PostStratum[ST.1$I.1.2+2*ST.1$I.1.3 == 3] <- "II"
#----------------------------------------

#----------------------------------------
# Population Model to Stratum II
# yk= 30 + 20I21 - 10I23(1-I21)  + ek 
I.2.1 = rbinom(N.2, 1, p.2.1)
I.2.3 = rbinom(N.2, 1, p.2.3)
e.2   = rnorm(N.2, 0, sqrt(v.2))

# Response Stratum II
y.k   = mu.3 - (mu.3-mu.1)*I.2.1 - (mu.3-mu.2)*I.2.3*(1-I.2.1)+ e.2
ST.2  = data.frame(stratum=rep("II", N.2), y.k, I.2.1, I.2.3)

# Poststratification 
# PostStratum is a modification by misclassification
ST.2$PostStratum[ST.2$I.2.1+2*ST.2$I.2.3 == 0] <- "II"
ST.2$PostStratum[ST.2$I.2.1+2*ST.2$I.2.3 == 1] <- "I"
ST.2$PostStratum[ST.2$I.2.1+2*ST.2$I.2.3 == 2] <- "III"
ST.2$PostStratum[ST.2$I.2.1+2*ST.2$I.2.3 == 3] <- "I"
#----------------------------------------

#----------------------------------------
# Population Model to Stratum III
# yk= 20 + 10I32 + 30I31(1-I32)  + ek
I.3.1 = rbinom(N.3, 1, p.1.2)
I.3.2 = rbinom(N.3, 1, p.1.3)
e.3   = rnorm(N.3, 0, sqrt(v.3))

# Response Stratum III
# y.k = mu.1 - (mu.2-mu.3)*I.3.2 - (mu.3-mu.1)*I.3.1*(1-I.3.2)+ e.3 # OLD
y.k   = mu.2 + (mu.3-mu.2)*I.3.2 + (mu.1-mu.2)*I.3.1*(1-I.3.2)+ e.3
ST.3  = data.frame(stratum=rep("III", N.3), y.k, I.3.1, I.3.2)

# Poststratification 
# PostStratum is a modification by misclassification  
ST.3$PostStratum[ST.3$I.3.1+2*ST.3$I.3.2 == 0] <- "III" #"II"   OLD
ST.3$PostStratum[ST.3$I.3.1+2*ST.3$I.3.2 == 1] <- "I" #"II"  #"I"    OLD
ST.3$PostStratum[ST.3$I.3.1+2*ST.3$I.3.2 == 2] <- "II"   #"III"  OLD
ST.3$PostStratum[ST.3$I.3.1+2*ST.3$I.3.2 == 3] <- "II"  #"III"  OLD
#----------------------------------------

#-------------------------------------------------
# Artificial Population

# Join the stratum data sets
population = join_all(list(ST.1, ST.2, ST.3), by = 'stratum', type = 'full')

population$PostStratum = as.factor(population$PostStratum)

# Population - Totals by Strata - 
Total.pop = population %>%dplyr::group_by(stratum) %>% dplyr::summarize(Tot.pop = sum(y.k))

print(Total.pop)

# Total populational 
Total.Y = Total.pop %>% summarise_each (funs(sum) , -stratum)

print(Total.Y)

# Sample size of tranfering by misclassification 
# n.tmp <-  table(sample$stratum, sample$PostStratum)

# Total populational PostStratum
N.PsT=table(population$PostStratum)

# Stratum  total size PostStratum
N.PsT.1 = N.PsT[1]  # Highly Cropland - Stratum I
N.PsT.2 = N.PsT[2] # Cropland - Stratum II
N.PsT.3 = N.PsT[3] # Non-Cropland - Stratum III

#-------------------------------------------------


##------------ Parameter of sampling -------------------------------------
# Cases
# Sample size (Goiana dataset) by Stratum 

# # # # Size 60
# n.1 =  42 # Highly Cropland - Stratum I
# n.2 =  15 # Cropland - Stratum II
# n.3 =   3 # Non-Cropland - Stratum III


# # # # Size 90
# n.1 = 31 # Highly Cropland - Stratum I
# n.2 = 23 # Cropland - Stratum II
# n.3 = 36 # Non-Cropland - Stratum III

# # # # Size 120
n.1 = 42 # Highly Cropland - Stratum I
n.2 = 30 # Cropland - Stratum II
n.3 = 48 # Non-Cropland - Stratum III


# # # Size 180
# n.1 = 62 # Highly Cropland - Stratum I
# n.2 = 45 # Cropland - Stratum II
# n.3 = 73 # Non-Cropland - Stratum III


# sample size
size <- c(n.1, n.2, n.3)
#-------------------------------------------------


#-------------------------------------------------
# Monte Carlo Simulation

# Evaluate Estimators

# Stratum and Total
data.estimator= NULL


# Monte Carlo replics
R = 10000

for(replic in 1:R)
{

  # Extract the sample (Simple random sampling of units by stratum )
  st <- strata(population, stratanames = "stratum", size = size, method = "srswor")
  
  # Sample data of interest
  sample = getdata(population, st)
  
  #-------------Corrected sample sizes -----------------
  # Poststratification 
  
  # # Sample size of tranfering by misclassification 
  # n.1.2 = sum(sample$I.1.2, na.rm=T)
  # n.1.3 = sum(sample$I.1.3, na.rm=T)
  # n.2.1 = sum(sample$I.2.1, na.rm=T)
  # n.2.3 = sum(sample$I.2.3, na.rm=T)
  # n.3.1 = sum(sample$I.3.1, na.rm=T)
  # n.3.2 = sum(sample$I.3.2, na.rm=T)
  # 
  
  
  # Sample size of tranfering by misclassification 
  n.tmp <-  table(sample$stratum, sample$PostStratum)
  
  # sample %>%dplyr::group_by(stratum, PostStratum) %>%
  # dplyr::summarize(ntran = n()) %>% complete(stratum, PostStratum)
  
  n.1.1 = n.tmp[1,1]
  n.1.2 = n.tmp[1,2]
  n.1.3 = n.tmp[1,3]
  
  n.2.1 = n.tmp[2,1]
  n.2.2 = n.tmp[2,2]
  n.2.3 = n.tmp[2,3]
  
  n.3.1 = n.tmp[3,1]
  n.3.2 = n.tmp[3,2]
  n.3.3 = n.tmp[3,3]
  
  
  # Corrected population size (unweight)
  N.1.s = N.1 - n.1.2 - n.1.3 + n.2.1 + n.3.1
  N.2.s = N.2 - n.2.1 - n.2.3 + n.1.2 + n.3.2
  N.3.s = N.3 - n.3.1 - n.3.2 + n.1.3 + n.2.3
  
  # # Corrected sample size  (Marginais coluna corretas)
  # n.1.s = (n.1-n.1.2-n.1.3) + n.2.1 + n.3.1 # n.1 + n.2.1 + n.3.1  
  # n.2.s = (n.2-n.2.1-n.2.3) + n.1.2 + n.3.2 # n.2 + n.1.2 + n.3.2  
  # n.3.s = (n.3-n.3.1-n.3.2) + n.1.3 + n.2.3 # n.3 + n.1.3 + n.2.3  
  
  # Corrected sample size  (Marginais coluna corretas)
  n.1.s = n.1.1 + n.2.1 + n.3.1 # n.1 + n.2.1 + n.3.1  
  n.2.s = n.2.2 + n.1.2 + n.3.2 # n.2 + n.1.2 + n.3.2  
  n.3.s = n.3.3 + n.1.3 + n.2.3 # n.3 + n.1.3 + n.2.3  
  
  
  # Corrected population size (weight)
  N.1.w = N.1 - (n.1.2/n.1)*N.1 - (n.1.3/n.1)*N.1 + (n.2.1/n.2)*N.2 + (n.3.1/n.3)*N.3
  N.2.w = N.2 - (n.2.1/n.2)*N.2 - (n.2.3/n.2)*N.2 + (n.1.2/n.1)*N.1 + (n.3.2/n.3)*N.3
  N.3.w = N.3 - (n.3.1/n.3)*N.3 - (n.3.2/n.3)*N.3 + (n.1.3/n.1)*N.1 + (n.2.3/n.2)*N.2
  
  # Corrected population size (expansion)
  
  # Total sample size
  ns.t = n.1+n.2+n.3
  
  N.ex.1 = (NT/ ns.t)*(n.1.1+n.2.1+n.3.1)
  N.ex.2 = (NT/ ns.t)*(n.1.2+n.2.2+n.3.2)
  N.ex.3 = (NT/ ns.t)*(n.1.3+n.2.3+n.3.3)
    
  
  #-------------------------------------------
  
  #------------------------------------------
  # Weights
  
  # Naive
  w.n.1 = N.1/n.1
  w.n.2 = N.2/n.2
  w.n.3 = N.3/n.3
  w.naive = c(w.n.1, w.n.2, w.n.3)
  
  # Unweighted
  w.u.1 = N.1.s/n.1.s
  w.u.2 = N.2.s/n.2.s
  w.u.3 = N.3.s/n.3.s
  w.unw = c(w.u.1, w.u.2, w.u.3)
  
  # Weighted
  w.w.1 = N.1.w/n.1.s
  w.w.2 = N.2.w/n.2.s
  w.w.3 = N.3.w/n.3.s
  w.w = c(w.w.1, w.w.2, w.w.3)
  
  # Weigthed Post Stratum
  w.p.1 = N.PsT.1/n.1.s
  w.p.2 = N.PsT.2/n.2.s
  w.p.3 = N.PsT.3/n.3.s
  w.p = c(w.p.1, w.p.2, w.p.3)
  
  # Expansion
  w.e.1 = N.ex.1/n.1
  w.e.2 = N.ex.2/n.2
  w.e.3 = N.ex.3/n.3
  w.e = c(w.e.1, w.e.2, w.e.3)
  #-----------------------------------------------------
  
  # ------------ Total by Stratum ---------------------------------
  # Total obtained by correct clasification  
  Total.correct <-   sample %>%dplyr::group_by(stratum) %>%
    dplyr::summarize(Tk = sum(y.k))
  Total.correct$type = "After"
  
  
  # Poststratification 
  # Total obtained by misclasification  
  Total.mis <-   sample %>%dplyr::group_by(PostStratum) %>% 
    dplyr::summarize(Tk = sum(y.k))
  colnames(Total.mis)[1] <-"stratum"
  Total.mis$stratum = as.factor(Total.mis$stratum)
  Total.mis$type = "Post"
  #-------------------------------------------------------
  
  # --------- Estimators ---------------------------------
  # The Horvitz-Thompson estimator for the total 
  
  # -------------- By Stratum---------------------------------
  # By omission of missclasification
  
  # Naive
  Total.correct$Basic = Total.correct$Tk*w.naive
  
  # Unweight
  Total.correct$Unweigthed = Total.correct$Tk*w.unw
  
  # Weight
  Total.correct$Weigthed = Total.correct$Tk*w.w
  
  # Weight Post
  Total.correct$PostStrat = Total.correct$Tk*w.p
  
  # Weight Expansion
  Total.correct$Expansion = Total.correct$Tk*w.e
  
  
  
  # By consider missclasification - Poststratification 
  
  # Naive
  Total.mis$Basic = Total.mis$Tk*w.naive
  
  # Unweigth
  Total.mis$Unweigthed = Total.mis$Tk*w.unw
  
  # Weight
  Total.mis$Weigthed = Total.mis$Tk*w.w
  
  # Weight Post
  Total.mis$PostStrat = Total.mis$Tk*w.p
  
  # Weight Expansion
  Total.mis$Expansion = Total.correct$Tk*w.e
  #-------------------------------------------------------
  
  # -------------- Total ---------------------------------
  # Estimator of Total Population
  # After Postratification
  
  # Reorder columns
  Total.correct = Total.correct[, c(3,1,2,4:8)]
  Tot.pop.cor = Total.correct %>% summarise_each (funs(sum) , -stratum, -type)
  
  #-----------------------------------------------------------
  # I dont like this solution to add and reorder the data set 
  Tot.pop.cor$stratum = "Total"
  Tot.pop.cor$stratum = as.factor(Tot.pop.cor$stratum)
  Tot.pop.cor$type = "After"
  Tot.pop.cor = Tot.pop.cor[, c(8,7, 1:6)]
  #-----------------------------------------------------------
  
  # Postratification
  # Reorder columns    
  Total.mis = Total.mis[, c(3,1,2,4:8)]    
  Tot.pop.mis = Total.mis %>% summarise_each (funs(sum) , -stratum, -type)
  
  #-----------------------------------------------------------
  # I dont like this solution to add and reorder the data set
  Tot.pop.mis$stratum = "Total"
  Tot.pop.mis$stratum = as.factor(Tot.pop.mis$stratum)
  Tot.pop.mis$type = "Post"
  Tot.pop.mis = Tot.pop.mis[, c(8,7, 1:6)]
  
  #-----------------------------------------------------------

    
    #-------------------------------------------------------

# added to data set to Long format
  
    a <- rbind(Total.correct,Tot.pop.cor)
    a$Nws = c(N.1.w,N.2.w,N.3.w,N.1.w+N.2.w+N.3.w)
    
    b <- rbind(Total.mis, Tot.pop.mis) 
    b$Nws = c(N.1.w,N.2.w,N.3.w,N.1.w+N.2.w+N.3.w)
    
    temp.1 <- gather(a, Estimator, Area, Tk:Expansion, factor_key = TRUE)
    temp.2 <- gather(b, Estimator, Area, Tk:Expansion, factor_key = TRUE)
    
    data.estimator =rbind(data.estimator, cbind(rbind(temp.1,temp.2),replic))

}

# Data sets by evaluate EQM and Relative Bias

# By Stratum
data.estimator


#---------------------------------------------------------------------
#-------- BIAS - VARIANCE - EQM  -------------------------------------
#--------------- For Tables  -----------------------------------------
#---------------------------------------------------------------------

# Total by Strata
data.estimator$PopulPar=unlist(c(Total.pop$Tot.pop, Total.Y))


# MC.temp <-   data.estimator %>%dplyr::group_by(type, stratum, Estimator) %>% 
#   dplyr::summarize( Estimador = mean(Area) ,
#                     Bias = mean((Area-PopulPar)),
#                     Rel.bias= mean((Area-PopulPar)/PopulPar), 
#                     EQM = mean((Area-PopulPar)^2), 
#                     SD = sd(Area), 
#                     Size = mean(Nws))
# 
# MC.output <- subset(MC.temp, (Estimator !="Tk") )
# 
# 
# # wide format
# data_wide <- dcast(setDT(MC.output), stratum + Estimator ~ type, 
#                    value.var=list("Estimador", "Bias", "Rel.bias", "EQM", "SD") )
# 
# # Export the data
# write.table(data_wide,file="SaidaMC.csv",quote=T,append=F,sep=",",eol = "\n", 
#             na = "NA", dec = ".", row.names = T,col.names = T)
# #---------------------------------------------------------------------
# #---------------------------------------------------------------------




#---------------------------------------------------------------------
#--------------------- Graphics --------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
# Grafic by Total
# Select the total
# g.t=subset(data.estimator, (stratum == "Total" & Estimator !="Tk") )

g.t.1=subset(data.estimator, 
             (type=="After" & Estimator=="Basic" & stratum == "Total") )

g.t.2=subset(data.estimator, 
             (type=="Post" & Estimator != "Basic" & Estimator != "Tk"  & stratum == "Total"))

g.t = rbind(g.t.1, g.t.2)


MC.temp <-   g.t %>%dplyr::group_by(type, stratum, Estimator) %>% 
  dplyr::summarize( Estimador = mean(Area) ,
                    Bias = mean((Area-PopulPar)),
                    Rel.bias= mean((Area-PopulPar)/PopulPar), 
                    EQM = mean((Area-PopulPar)^2), 
                    SD = sd(Area), 
                    Size = mean(Nws))

MC.output <- subset(MC.temp, (Estimator !="Tk") )

# 
# # wide format
# data_wide <- dcast(setDT(MC.output), stratum + Estimator ~ type, 
#                    value.var=list("Estimador", "Bias", "Rel.bias", "EQM", "SD") )

# Export the data
# write.table(MC.output,file="SaidaTotal.csv") 
# 
# ,quote=T,append=F,sep=",",eol = "\n", 
#             na = "NA", dec = ".", row.names = T,col.names = T)


p <-  ggplot(g.t, aes(x=Estimator,  y=Area, fill=Estimator))
p <-  p+geom_boxplot(aes(color=Estimator),  notch = FALSE, alpha=0.3, width=0.5, 
                     outlier.colour = "darkblue", outlier.shape = 18, 
                     outlier.size = 5) #+facet_wrap(~stratum)
p <- p <- p+labs(x = " ", y="Area", title="Total")
p <- p +  theme(
  legend.key = element_rect(colour = "white", 
                            fill = "white"),
  legend.key.size = unit(1.1, "cm"),
  legend.text = element_text(face = "bold", 
                             size=15),
  legend.title = element_text(size=15, face="bold"),
  panel.grid.major = element_line(colour = "gray", 
                                  linetype = "dotted"),
  # panel.grid.minor = element_line(colour = "red", linetype = "dotted"),
  panel.background = element_rect(fill = "white", 
                                  colour="black"),
  strip.text.x = element_text(size=12, 
                              hjust=0.5, 
                              vjust=0.5,
                              face="bold", lineheight = 0.5),
  strip.background = element_rect(colour="black", fill="gray98"),
  axis.text=element_text(size=15, face="bold", colour="gray24"),
  
  axis.text.x  = element_text(angle=0, vjust=0.5, size=12),
  axis.title=element_text(size=15,face="bold"),
  plot.title = element_text(size = 13, colour = "black", face="bold", hjust=0.5),
  # legend.position = c(0.85, 0.85)
  legend.position ="right", #c(0.1, 0.85)
  legend.direction="vertical"
    
) 
pdf("Totales-n-60.pdf", width=10, height = 10)
p <-  p+geom_hline(aes(yintercept=Total.Y), color="darkblue") # Reference
# p <-  p+annotate("text", x = 0.5, y = as.numeric(paste(Total.Y))+1000, 
#                  label="Truth", parse=TRUE, col="darkred", size=6)
p
dev.off()

#---------------------------------------------------------------------




#---------------------------------------------------------------------
#---------------------------------------------------------------------
# Grafic by Stratum
# Select the Stratum
# g.t=subset(data.estimator, (stratum !="Total" & Estimator !="Tk") )

g.s.1=subset(data.estimator, 
             (type=="After" & Estimator=="Basic" & stratum != "Total") )

g.s.2=subset(data.estimator, 
             (type=="Post" & Estimator != "Basic" & Estimator != "Tk"  
              & stratum != "Total"))

g.s = rbind(g.s.1, g.s.2)



MC.temp2 <-   g.s %>%dplyr::group_by(type, stratum, Estimator) %>% 
  dplyr::summarize( Estimador = mean(Area) ,
                    Bias = mean((Area-PopulPar)),
                    Rel.bias= mean((Area-PopulPar)/PopulPar), 
                    EQM = mean((Area-PopulPar)^2), 
                    SD = sd(Area), 
                    Size = mean(Nws))

MC.output2 <- subset(MC.temp2, (Estimator !="Tk") )


# # wide format
# data_wide <- dcast(setDT(MC.output), stratum + Estimator ~ type, 
#                    value.var=list("Estimador", "Bias", "Rel.bias", "EQM", "SD") )

# Export the data
# write.table(data_wide,file="SaidaMC.csv",quote=T,append=F,sep=",",eol = "\n", 
#             na = "NA", dec = ".", row.names = T,col.names = T)


# h.1=data.estimator[(data.estimator$type=="After" & data.estimator$Estimator=="Basic" & 
#                      data.estimator$stratum != "Total"),]
# 
# h.2=data.estimator[(data.estimator$type=="Post" & 
#                       data.estimator$Estimator == "Unweigthed" & 
#                       data.estimator$Estimator == "Weigthed"  & data.estimator$stratum != "Total"),]
# 

p <-  ggplot(g.s, aes(x=Estimator,  y=Area, fill=Estimator))
p <-  p+geom_boxplot(aes(color=Estimator),  notch = FALSE, alpha=0.3, width=0.5, 
                     outlier.colour = "darkblue", outlier.shape = 18, 
                     outlier.size = 2)+facet_wrap(~stratum)
p<-p+geom_hline(aes(yintercept = PopulPar),color="darkblue")

p <- p <- p+labs(x = " ", y="Area", title="Strata")
p <- p +  theme(
  legend.key = element_rect(colour = "white", 
                            fill = "white"),
  legend.key.size = unit(1.1, "cm"),
  legend.text = element_text(face = "bold", 
                             size=15),
  legend.title = element_text(size=15, face="bold"),
  panel.grid.major = element_line(colour = "gray", 
                                  linetype = "dotted"),
  # panel.grid.minor = element_line(colour = "red", linetype = "dotted"),
  panel.background = element_rect(fill = "white", 
                                  colour="black"),
  strip.text.x = element_text(size=12, 
                              hjust=0.5, 
                              vjust=0.5,
                              face="bold", lineheight = 0.5),
  strip.background = element_rect(colour="black", fill="gray98"),
  axis.text=element_text(size=15, face="bold", colour="gray24"),
  
  axis.text.x  = element_text(angle=0, vjust=0.5, size=12),
  axis.title=element_text(size=15,face="bold"),
  plot.title = element_text(size = 13, colour = "black", face="bold", hjust=0.5),
  # legend.position = c(0.11, 0.15)
  legend.position ="right", #c(0.1, 0.85)
  legend.direction="vertical"
) 
pdf("Stratos-n-60.pdf", width=18, height = 8)
p
dev.off()


#---------------------------------------------------------------------



#---------------------------------------------------------------------
#-------- Sample size behaviour  -------------------------------------
#------------- PLOT ---------------------------------------------------

# Out MC

data.estimator$Nh=unlist(c(N.1, N.2, N.3, N.1+N.2+N.3 ))
data.estimator$Nhname=unlist(c("n.1", "n.2", "n.3", "n" ))
data.estimator$NhPost=unlist(c(N.PsT.1, N.PsT.2, N.PsT.3, N.PsT.1+N.PsT.2+N.PsT.3 ))

n.sample <- subset(data.estimator, stratum !="Total" )

p <-  ggplot(n.sample, aes(x=stratum,  y=Nws, fill=stratum))
p <-  p+geom_boxplot(aes(color=stratum),  notch = FALSE, alpha=0.3, width=0.5, 
                     outlier.colour = "darkblue", outlier.shape = 18, 
                     outlier.size = 5)+facet_wrap(~Nhname)+guides(color=FALSE)
p <-  p+geom_hline(aes(yintercept=Nh, color="darkblue")) # Reference
p <- p <- p+labs(x = " ", y="sample size", title=" ", fill="Strata")
p <- p +  theme(
  legend.key = element_rect(colour = "white", 
                            fill = "white"),
  legend.key.size = unit(1.1, "cm"),
  legend.text = element_text(face = "bold", 
                             size=15),
  legend.title = element_text(size=15, face="bold"),
  panel.grid.major = element_line(colour = "gray", 
                                  linetype = "dotted"),
  # panel.grid.minor = element_line(colour = "red", linetype = "dotted"),
  panel.background = element_rect(fill = "white", 
                                  colour="black"),
  strip.text.x = element_text(size=12, 
                              hjust=0.5, 
                              vjust=0.5,
                              face="bold", lineheight = 0.5),
  strip.background = element_rect(colour="black", fill="gray98"),
  axis.text=element_text(size=15, face="bold", colour="gray24"),
  
  axis.text.x  = element_text(angle=0, vjust=0.5, size=12),
  axis.title=element_text(size=15,face="bold"),
  plot.title = element_text(size = 13, colour = "black", face="bold", hjust=0.5),
  legend.position ="right", #c(0.1, 0.85)
  legend.direction="vertical"
  
  
) 
pdf("n-size-60.pdf", width=11, height = 8)
p
dev.off()


# teste=join_all(list(MC.output2, MC.output), by = 'type')
teste=rbind(MC.output2, MC.output)

write.table(teste,file="SaidaTotal-60.csv") 
