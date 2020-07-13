rm(list=ls())

#### Load packages 
library(dplyr)
library(survival)
library(data.table)
library(ggplot2)
library(questionr)  # oddsratio
library(networkD3)  # graphic multistate 
library(htmlwidgets)# save widgets

#### Read data 
setwd("C:/Users/jordi/Dropbox/Retreat GRBIO - Modeling/bbdd")
load("datos.RData")

####################################################################################################
# names(d_patients)
# d_patients <- d_patients %>% mutate_at(c("f_ingreso","f_entrada_uci","f_salida_uci","f_alta","f_ingreso_urg"), as.Date, format="%d/%m/%Y")

###########################################
# Noves variables
###########################################
#### Ingreso
d_patients$time_fu <- as.numeric(d_patients$f_alta-d_patients$f_ingreso)
d_patients$status_fu <-ifelse(d_patients$motivo_alta=="Fallecimiento",1,0)
d_patients$os <- Surv(d_patients$time_fu,d_patients$status_fu)
plot(d_patients$os)

# Com tractar pacients donats d'alta?
d_patients$time_fu[d_patients$motivo_alta=="Domicilio"] <- 30
d_patients$status_fu[d_patients$status_fu==1 &d_patients$time_fu>30 ] <-0
d_patients$time_fu[d_patients$time_fu>30 ] <-30
d_patients$os <- Surv(d_patients$time_fu,d_patients$status_fu)
plot(d_patients$os)
summary(survfit(d_patients$os ~1))                 # 83.6% OS at 30 days

#### UCI
d_patients$time_to_UCI <- as.numeric(d_patients$f_entrada_uci-d_patients$f_ingreso)
summary(d_patients$time_to_UCI)                    # 1 dia de mediana des de diagnòstic fins entrada UIC
sum(!is.na(d_patients$uci_dias))/nrow(d_patients)  # 8.6% van entrar a UCI
summary(d_patients$uci_dias)                       # 4 dies de mediana a UCI

###########################################
# Models de Cox - Anàlisi supervivencia
###########################################
d_patients$edad10 <- d_patients$edad/10
cox <- coxph(os ~ edad10 ,data=d_patients);summary(cox)
cox <- coxph(os ~ sexo ,data=d_patients);summary(cox)
cox <- coxph(os ~ h_admision ,data=d_patients);summary(cox)
cox <- coxph(os ~ h_const_primera ,data=d_patients);summary(cox)
cox <- coxph(os ~ temp_primera ,data=d_patients);summary(cox)
cox <- coxph(os ~ fc_primera ,data=d_patients);summary(cox)
cox <- coxph(os ~ sat_02_primera ,data=d_patients);summary(cox)   # Si                                                   
cox <- coxph(os ~ ta_max_primera ,data=d_patients);summary(cox)
cox <- coxph(os ~ ta_min_primera ,data=d_patients);summary(cox)   # Si

###########################################
# Models Logístics - Anàlisi UCI
###########################################
d_patients$UCI <- ifelse(is.na(d_patients$f_entrada_uci), 0,1)
mod <- glm(UCI ~ edad10 , family = "binomial",data=d_patients);odds.ratio(mod) #
mod <- glm(UCI ~ sexo , family = "binomial",data=d_patients);odds.ratio(mod)   #
mod <- glm(UCI ~ h_admision , family = "binomial",data=d_patients);odds.ratio(mod)
mod <- glm(UCI ~ h_const_primera , family = "binomial",data=d_patients);odds.ratio(mod)  
mod <- glm(UCI ~ temp_primera , family = "binomial",data=d_patients);odds.ratio(mod)    
mod <- glm(UCI ~ fc_primera  , family = "binomial",data=d_patients);odds.ratio(mod)
mod <- glm(UCI ~ sat_02_primera, family = "binomial",data=d_patients);odds.ratio(mod) 
mod <- glm(UCI ~ ta_max_primera  , family = "binomial",data=d_patients);odds.ratio(mod)
mod <- glm(UCI ~ ta_min_primera, family = "binomial",data=d_patients);odds.ratio(mod)


#### Interacció amb edat
mod <- glm(UCI ~ sexo*edad10 , family = "binomial",data=d_patients);odds.ratio(mod)  
mod <- glm(UCI ~ h_admision*edad10 , family = "binomial",data=d_patients);odds.ratio(mod)
mod <- glm(UCI ~ h_const_primera*edad10 , family = "binomial",data=d_patients);odds.ratio(mod)  
mod <- glm(UCI ~ temp_primera*edad10 , family = "binomial",data=d_patients);odds.ratio(mod)    
mod <- glm(UCI ~ fc_primera*edad10  , family = "binomial",data=d_patients);odds.ratio(mod)
mod <- glm(UCI ~ sat_02_primera*edad10, family = "binomial",data=d_patients);odds.ratio(mod)  
mod <- glm(UCI ~ ta_max_primera*edad10  , family = "binomial",data=d_patients);odds.ratio(mod)
mod <- glm(UCI ~ ta_min_primera*edad10, family = "binomial",data=d_patients);odds.ratio(mod)


###########################################
# Merge patients and medications
###########################################
# Labels
ATC7 <- c("L04AC07","H02AB09","J05AE20","P01BA02","L03AB08","V03AN01","V03AN04")
VAR <- c('Tocilizumab','Dexametasona','Antivirals','Hidroxicloriquina','Interferons','Gases')

# Medications
dmed0 <- data.table::dcast(d_medications %>% filter(id_atc7 %in% ATC7),patient_id~id_atc7)
# d_medications %>% filter(patient_id  %in% 9)
dmed0$V03AN <- dmed0$V03AN01 + dmed0$V03AN04
dmed0$V03AN01 <- NULL 
dmed0$V03AN04 <- NULL
toYes <- function(x) ifelse(x,'Yes','No')
dmed <- as.data.frame(cbind(dmed0[,1],apply(dmed0[,2:7],2,toYes)))
names(dmed) <- c('patient_id',VAR[c(2,3,5,1,4,6)])
dpat2 <- merge(d_patients,dmed,by='patient_id',all.x=TRUE)

#####################################################################################################################

### Plots:

##############################
### 1) Edat, mortalitat i sexe
##############################

d_patients$edad
d_patients$edad_group <- NA
d_patients$edad_group[d_patients$edad<50] <- "1.Menys 50"
d_patients$edad_group[d_patients$edad>=50 & d_patients$edad< 60 ] <- "2.50-60"
d_patients$edad_group[d_patients$edad>=60 & d_patients$edad< 70 ] <- "3.60-70"
d_patients$edad_group[d_patients$edad>=70 & d_patients$edad< 80 ] <- "4.70-80"
d_patients$edad_group[d_patients$edad>=80] <- "5.>80"
table(d_patients$edad_group)

table(d_patients$sexo)
women <-prop.table(table(d_patients$edad_group[d_patients$sexo=="FEMALE"],d_patients$status_fu[d_patients$sexo=="FEMALE"]),1)*100
men <- prop.table(table(d_patients$edad_group[d_patients$sexo=="MALE"],d_patients$status_fu[d_patients$sexo=="MALE"]),1)*100
women[,2]

windows(18,10)
barplot(women[,2],horiz = T, col = "#613D2D",xlim=c(-60,60))
barplot(-men[,2],add = TRUE,horiz = T, axes = F, col = "#D9CCB9")

##############################
### 2) Edat, mortalitat i UCI
##############################

mortality <-prop.table(table(d_patients$edad_group,d_patients$status_fu),1)*100
d_patients$uci <- ifelse(!is.na(d_patients$uci_dias),1,0)
table(d_patients$uci)
uci_perc <-prop.table(table(d_patients$edad_group,d_patients$uci),1)*100

cols1 <- RColorBrewer::brewer.pal(3, "Paired")
windows(18,10)
barplot(mortality[,2],horiz = T, col = cols1[2],xlim=c(-60,60))
barplot(-uci_perc[,2],add = TRUE,horiz = T, axes = F, col = cols1[1])

##############################
### 3) Forestplot OR
##############################
library(metafor)
#install.packages("meta")
library(meta)
library(grid)

setwd("~/Copia13.01.20/Retreat")
dat<-read.csv2("Foresplot.csv", sep=",",dec=".")

# Effect size (HR)
loghr <- log(dat$OR)
seloghr <- (log(dat$CI95..sup)-log(dat$CI95..inf))/3.92   #3.92 =1.96*2
dat <- cbind(dat, loghr, seloghr)

## Meta-anàlisi 1
meta1 <- metagen(loghr, seloghr, data=dat, studlab=paste(Variable), sm="OR")
print(meta1)
windows(20,10) 
forest(meta1, col.diamond="black", fontsize=13,  squaresize=1, plotwidth=unit(10, "cm"),leftcols=c("studlab"), rightcols=c("effect", "ci"),
       bysort=T,  comb.fixed=FALSE, ,spacing = 1.4,comb.random = F,text.random="Overall", TE.random=FALSE, seTE.random=FALSE,xlim=c(0.2,5),at=c(0.2,0.5,1,2,5)) 


###########################################
# 4. Graphics Florence Nithingale
###########################################
common_theme <- theme(axis.ticks = element_blank(),
                      axis.text.y = element_blank(),
                      axis.text.x = element_text(size=7,face='bold'))

##-- Treat x Sex
dpat3 <- expand.grid(Treat= VAR,Sex=c('MALE','FEMALE'))
dpat3$N <- dpat3$TOT <- NA
for(i in 1:nrow(dpat3)){
  dpat3$N[i] <- sum(dpat2$sexo==dpat3[i,2] & dpat2[,as.character(dpat3[i,1])]=='Yes',na.rm=TRUE)
  dpat3$TOT[i] <- sum(dpat2[,as.character(dpat3[i,1])]=='Yes',na.rm=TRUE)
}
dpat3$lab <- paste0(dpat3$N,' (',formatC(100*dpat3$N/dpat3$TOT,digits=0,format='f'),'%)')
ggplot(dpat3[dpat3$Treat!='Gases',],aes(x=Treat,y=N,fill=Sex,label=lab)) +
  geom_bar(width = 1, position="stack",stat="identity", color="black") +
  geom_text(size = 1, position = position_stack(vjust = 0.5)) +
  scale_y_sqrt() + coord_polar(start=-pi/2.5) + 
  ggtitle("Sexo según fármaco") + 
  xlab("") + ylab("") + common_theme

ggsave(filename = '../Plots/FN_Treat_x_Sex.png',width = 20,units = 'cm')

##-- Treat x Death at 30 days
dpat3 <- expand.grid(Treat=VAR,Death=0:1)
dpat3$N <- dpat3$TOT <- NA
for(i in 1:nrow(dpat3)){
  dpat3$N[i] <- sum(dpat2$status_fu==dpat3[i,2] & dpat2[,as.character(dpat3[i,1])]=='Yes',na.rm=TRUE)
  dpat3$TOT[i] <- sum(dpat2[,as.character(dpat3[i,1])]=='Yes',na.rm=TRUE)
  
}
dpat3$Death <- as.factor(ifelse(as.numeric(as.character(dpat3$Death)),'Yes','No'))
dpat3$lab <- paste0(dpat3$N,' (',formatC(100*dpat3$N/dpat3$TOT,digits=0,format='f'),'%)')
ggplot(dpat3[dpat3$Treat!='Gases',],aes(x=Treat,y=N,fill=Death,label=lab)) +
  geom_bar(width = 1, position="stack",stat="identity", color="black") +
  geom_text(size = 1, position = position_stack(vjust = 0.5)) +
  scale_y_sqrt() + coord_polar(start=-pi/2.5) + 
  ggtitle("Mortalidad según fármaco") + 
  xlab("") + ylab("") + common_theme

ggsave(filename = '../Plots/FN_Treat_x_Death.png',width = 20,units = 'cm')


##-- Treat x O2 x Death
dpat2$O2 <- cut(dpat2$sat_02_primera,breaks = c(0,90,95,100),include.lowest = TRUE)
dpat3 <- expand.grid(VAR,0:1,c('[0,90]', '(90,95]', '(95,100]'))
names(dpat3) <- c('Treat','Death','SAT')
dpat3$N <- dpat3$TOT <- NA
for(i in 1:nrow(dpat3)){
  dpat3$N[i] <- sum(dpat2$status_fu==dpat3[i,2] & dpat2$O2==dpat3[i,3] & dpat2[,as.character(dpat3[i,1])]=='Yes',na.rm=TRUE)
  dpat3$TOT[i] <- sum(dpat2$O2==dpat3[i,3] & dpat2[,as.character(dpat3[i,1])]=='Yes',na.rm=TRUE)
  
}
dpat3$Death <- as.factor(ifelse(as.numeric(as.character(dpat3$Death)),'Yes','No'))

dpat3$lab <- paste0(dpat3$N,' (',formatC(100*dpat3$N/dpat3$TOT,digits=0,format='f'),'%)')
ggplot(dpat3[dpat3$Treat!='Gases' & dpat3$Treat!='Dexametasona' & dpat3$Treat!='Interferons',],
              aes(x=Treat,y=N,fill=Death,label=lab)) +
  geom_bar(width = 1, position="stack", 
           stat="identity", color="black") +
  geom_text(size = 1, position = position_stack(vjust = 0.5)) +
  facet_wrap(.~SAT) +
  scale_y_sqrt() + coord_polar(start=-pi/2.5) + 
  ggtitle("Muerte según fármaco y saturación basal") + 
  xlab("") + ylab("") +
  theme(strip.text.x = element_text(size = 15, face='bold'),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size=7,face='bold'))

ggsave(filename = '../Plots/FN_Treat_x_Death_x_O2.png',width = 40,units = 'cm')


##-- Treat x Age
dpat2$age_cut <- cut(dpat2$edad,breaks = c(0,50,60,70,80,110),include.lowest = TRUE)
dpat3 <- expand.grid(Treat=VAR,Age=c('[0,50]','(50,60]','(60,70]','(70,80]','(80,110]'))
dpat3$N <- dpat3$TOT <- NA
for(i in 1:nrow(dpat3)){
  dpat3$N[i] <- sum(dpat2$age_cut==dpat3[i,2] & dpat2[,as.character(dpat3[i,1])]=='Yes',na.rm=TRUE)
  dpat3$TOT[i] <- sum(dpat2[,as.character(dpat3[i,1])]=='Yes',na.rm=TRUE)
}
dpat3$lab <- paste0(dpat3$N,' (',formatC(100*dpat3$N/dpat3$TOT,digits=0,format='f'),'%)')
ggplot(dpat3[dpat3$Treat!='Gases',],aes(x=Treat,y=N,fill=Age,label=lab)) +
  geom_bar(width = 1, position="stack",stat="identity", color="black") +
  geom_text(size = 1, position = position_stack(vjust = 0.5)) +
  scale_y_sqrt() + coord_polar(start=-pi/2.5) + 
  ggtitle("Edad según fármaco") + 
  xlab("") + ylab("") + common_theme

ggsave(filename = '../Plots/FN_Treat_x_Age.png',width = 20,units = 'cm')



###############################################################################################################################

############################
#### Naive prognostic model for 30-days mortality
############################

library("OptimalCutpoints")
library(mice)
set.seed(18)
imp=mice(d_patients[,c(1:2,13:25)],m=1,method = "pmm")
d_patients[,c(1:2,13:25)] <- complete(imp) 
d_patients$edad10 <- d_patients$edad/10
dd <- subset(d_patients, !is.na(status_fu))

optimal.cutpoint.Youden <- optimal.cutpoints(X = "edad10", status = "status_fu", tag.healthy = 0,
                                             methods = "Youden", data = dd, pop.prev = NULL, 
                                             control = control.cutpoints(), ci.fit = FALSE, conf.level = 0.95, trace = FALSE)

library(survminer)
## Punt de tall per edat
res.cut <- surv_cutpoint(dd, time = "time_fu", event = "status_fu",variables = c("edad10"),minprop=0.25)
cutoff<-summary(res.cut)   # Punto Ã³ptimo 2.5 (el definido actualmente es 3)
dat$age_cut <- ifelse(dat$edad10 > 7.5,1,0) ## 75

## Punt de tall per sat
res.cut <- surv_cutpoint(dd, time = "time_fu", event = "status_fu",variables = c("sat_02_primera"),minprop=0.25)
cutoff<-summary(res.cut)   # Punto Ã³ptimo 2.5 (el definido actualmente es 3)
dd$sat_cut <- ifelse(dd$sat_02_primera > 90,0,1) ## 90

## Calculem el score 0,1,2
d_patients$age_cut <- ifelse(d_patients$edad10 > 7.5,1,0) ## 75
d_patients$sat_cut <- ifelse(d_patients$sat_02_primera > 90,0,1) ## 75
d_patients$group_risk <-d_patients$age_cut +d_patients$sat_cut 
x <- prop.table(table(d_patients$group_risk,d_patients$status_fu ),1);x

barplot(x[,2],ylim=c(0,1),col=cols)

surv <- survfit(d_patients$os ~ d_patients$group_risk);surv
library(RColorBrewer)
cols <-brewer.pal(n = 8, name = "Reds")[c(3,5,8)]

library(survminer)
win.graph(16,10)
plot(surv1)
ggsurvplot(surv, data = d_patients,
           title = "",
           pval = F, pval.method = F,    # Add p-value &  method name
           palette = cols,                   # Use JCO journal color palette o AAAS palette, D3 etc.
           risk.table = T,                  # Add No at risk table
           cumcensor=F,
           tables.height = 0.15,               # Specify tables height
           tables.theme = theme_cleantable(),  # Clean theme for tables
           tables.y.text = F,
           conf.int = F, # Hide tables y axis text
           xlab= "Time (months)",
           ylab="Probability of event free (OS)",
           pval.size=4.5,
           risk.table.title="N. at risk",
           risk.table.fontsize=4.5,
           font.y=c(18),
           font.tickslab=18,
           size=1.9,
           font.x=c(18),
           linetype=c(1,1,1),
           legend=c(0.9,1.8),
           legend.title = "",           # Change legend titles
           legend.labs = c(" ", "  ","    "),  # Change legend labels
           font.legend=c(14,"bold"),
           break.time.by=3,
           xlim=c(0,30)) #surv.median.line="hv")


###############################################################################################################################

###########################################
# Transiciones multistate model
###########################################
#### Ns
# Patients to ICU
(a <- sum(!is.na(d_patients$f_entrada_uci)))
a/nrow(d_patients)

# Patients to Discharge without ICU
(a <- sum(is.na(d_patients$f_entrada_uci) & d_patients$motivo_alta=='Domicilio',na.rm = TRUE))
a/nrow(d_patients)

# Patients death without ICU
(a <- sum(is.na(d_patients$f_entrada_uci) & d_patients$motivo_alta=='Fallecimiento',na.rm = TRUE))
a/nrow(d_patients)

# Patients to Discharge with ICU
(a <- sum(!is.na(d_patients$f_entrada_uci) & d_patients$motivo_alta=='Domicilio',na.rm = TRUE))
a/sum(!is.na(d_patients$f_entrada_uci))

# Patients death with ICU
(a <- sum(!is.na(d_patients$f_entrada_uci) & d_patients$motivo_alta=='Fallecimiento',na.rm = TRUE))
a/sum(!is.na(d_patients$f_entrada_uci))


#### Plot multistate
# https://www.r-graph-gallery.com/321-introduction-to-interactive-sankey-diagram-2.html
links <- data.frame(
  source=c("Admission","Admission", "Admission", "ICU", "ICU"), 
  target=c("ICU","Death", "Discharge", "Death", "Discharge"), 
  value=c(217, 275, 1538, 58, 78)
)

nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1

p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE,fontSize=14,nodeWidth=30)
p

saveWidget(p, file=paste0( getwd(), "/HtmlWidget/sankeyBasic1.html"))
