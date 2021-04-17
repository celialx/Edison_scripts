```{r}
library(ggplot2)
library(ggpubr)
library(car)
library(haven)
library(data.table)
library(visreg)
library(dplyr)
library(tidyr)
library(reshape2)
library(rstatix)
library(JWileymisc)
library(lme4)
library(lmerTest)
library(multilevelTools)
library(tidyverse)
library(Rmisc)


SWS_data <- read.table("../Tables/CTET_SWdetection_thr90_byE_P2P_avDens_behav_vec_full_v3.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)

names(SWS_data)[names(SWS_data) == "Drug"] <- "Treatment"

SWS_data$Treatment<-as.factor(SWS_data$Treatment)

SWS_data$Treatment<-relevel(SWS_data$Treatment,"PLA")

SWS_data2<-aggregate(SWdens ~ SubID+Treatment+BlockN, data=SWS_data, FUN=mean)


```

```