# Load required libraries
library(metafor)
library(dplyr)
library(reshape2)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(readxl)
library(tidyr)
library(haven)
library(stringr)

# load WHO's data
dat_3d_WHO <- read_dta("WHO_data/ipvonly_studies_3doses.dta") %>%
  rename(N_1 = n1, N_2 = n2, N_3 = n3, p_1 = sc1p, p_2 = sc2p, p_3 = sc3p) %>%
  mutate(type = replace(type, type=="3IPV", "IPV"),
         type = replace(type, type=="3sIPV", "sIPV"),
         type = replace(type, type=="3fIPV", "fIPV"),
         schedule = replace(schedule, schedule=="8, 16, 24 weeks", "2, 4, 6 months"),
         year = as.numeric(str_sub(Journal_year, -4, -1))) %>%
  # clean up/shorten study name
  mutate(author = replace(author, author == "Luis Rivera et al. (Dominican Rep)", "Rivera L et al. (Dominican Rep)"),
         author = replace(author, author == "Josefina Cadorna-Carlos et al. (Philippines)", "Cadorna-Carlos J et al. (Philippines)"),
         author = replace(author, author == "Sonia Resik et al. (Cuba)", "Resik S et al. (Cuba)"),
         author = replace(author, author == "Gustavo Dayan et al. (Puerto Rico)", "Dayan G et al. (Puerto Rico)"),
         author = replace(author, author == "Jingjun Qui et al. (China)", "Qui J et al. (China)"),
         author = replace(author, author == "Guoyang Liao et al (China)", "Liao G et al (China)"),
         author = replace(author, author == "Guoyang Liao et al. (China)", "Liao G et al. (China)"),
         author = replace(author, author == "Yuemei Hu et al. (China)", "Hu Y et al. (China)"),
         author = replace(author, author == "Ali Saleem at el. (Pakistan)", "Saleem A et al. (Pakistan)"),
         author = replace(author, author == "Ananda et al. (Panama and Dominican Rep)", "Bandyopadhyay et al. (Panama & Dominican Rep)"),
         author = replace(author, author == "Chu kai et al. (China)", "Chu K et al. (China)"),
         author = replace(author, author == "Ali Jafer Mohammed et al. (Oman)", "Mohammed AJ et al. (Oman)"),
         author = replace(author, author == "Miguel O'Ryan et al. (Chile)", "O'Ryan M et al. (Chile)")
  ) %>%
  rename(author_country = author) %>%
  mutate(country = sub(".*\\(([^)]+)\\).*", "\\1", author_country),
         author = sub("[\\.\\(].*", "", author_country)) %>%
  mutate(country = replace(country, country == "Havana", "Cuba")) %>%
  select(-sc1,-sc2,-sc3,-slno,-schr,-ipvonly,-frac_full,-salk_sabin,-Journal_year,-author_country)


# load new data
new_dat0 <- read_xlsx("hexa_review.xlsx", 
                      na="", sheet = "extracted_data") %>% 
  filter(!is.na(study), !is.na(p)) %>% 
  rename(author = study) %>%
  rename(country = location) %>%
  mutate(author = paste(author," et al", sep="")) %>%
  mutate(author = replace(author, author == "Sanofi et al", "Sanofi study")) %>%
  mutate(type = ifelse(type=="Hexa", paste(type, " (", type_notes, ")", sep=""), type)) %>%
  select(author, year, country, study_id, schedule, type, num_doses, type_notes, strain, serotype, p, N, Notes) %>%
  filter(author != "LG Chem et al") %>% # confidential
  filter(!(author == "Capeding MR et al" & strain == "wild type"))

# format new data from long to wide
new_dat <- new_dat0 %>%
  pivot_wider(names_from = serotype, values_from = c(p, N)) 

# load OPV history data
vac_hist <- read_xlsx("C:/Users/qifangbi/OneDrive - Bill & Melinda Gates Foundation/Hexa/hexa_review.xlsx", na="", sheet = "OPV_history") %>%
  select(study_id, author, Journal_year, birthdose_OPV) %>%
  mutate(birthdose_OPV = replace(birthdose_OPV, birthdose_OPV != "N", "Y/NA")) %>%
  unique()
  
# merge OPV history with other data
dat_34d <- dat_3d_WHO %>%
  bind_rows(new_dat) %>% 
  # filter out WHO Thai to prevent duplicates
  left_join(vac_hist %>% filter(author!="WHO Collaborative study (Thailand)") %>% select(study_id,birthdose_OPV), 
            by = c("study_id")) %>%
  mutate(log_p_1 = log(p_1/100), se_log_p_1 = sqrt((1 - p_1/100) / (p_1/100 * N_1)),
         log_p_2 = log(p_2/100), se_log_p_2 = sqrt((1 - p_2/100) / (p_2/100 * N_2)), 
         log_p_3 = log(p_3/100), se_log_p_3 = sqrt((1 - p_3/100) / (p_3/100 * N_3)) 
         ) %>%
  mutate(author_year = paste(sub("\\s+$", "", author), year, sep=", ")) %>%
  arrange(num_doses, schedule, desc(type))

# correct a few mistakes in WHO data and filter out ineligible studies
dat_34d <- dat_34d %>%
  # remove two arms from Grassly
  filter(study_id != 108 | (study_id==108 & Notes == "infants born to mothers given no DTaP/IPV in pregnancy")) %>%
  # update WHO data - Luis Rivera et al. (Dominican Rep)
  mutate(p_1 = ifelse(study_id == 14, 100, p_1),
         N_1 = ifelse(study_id == 14, 206, N_1),
         N_2 = ifelse(study_id == 14, 203, N_2)) 


dat_34d %>% filter(num_doses == 4) %>% with(table(schedule, type))
dat_34d %>% filter(num_doses == 3) %>% with(table(schedule, type))

#write.csv(dat_34d, 'dat_34d.csv')

### rma method
tau_estim_method = "EB"


# 3 dose, IPV or IPV(hexa) only
dat_34d %>% 
  filter(num_doses == 3, type %in% c("IPV","Hexa")) %>% 
  with(table(schedule, type))

dat_3d_allIPV <- dat_34d %>% 
  filter(num_doses == 3, type %in% c("IPV")|str_detect(type, "Hexa")) 
  

dat_3d_allIPV %>% group_by(schedule) %>% summarise(count=n()) -> all_schedule_3d_allIPV

### plotting functions
# starting_row = Starting number
# gap_rows = Offset between sequences
gen_rma_rows <- function(c_vec, starting_row=1, gap_rows){
  result <- c()
  for (i in seq_along(c_vec)) {
    end <- starting_row + c_vec[i] - 1
    result <- c(result, seq(starting_row, end))
    if (i < length(c_vec)) {
      starting_row <- end + gap_rows + 1
    }
  }
  return(result)
}

gen_rma_schedule_rows <- function(c_vec, starting_row=1, gap_rows){
  result <- c()
  for (i in seq_along(c_vec)) {
    if (i == 1){
      result1 <- starting_row + c_vec[1]
      result <- c(result, result1)
    } 
    else{
      result1 <- result1 + gap_rows + c_vec[i] 
      result <- c(result, result1)
    }
  }
  return(result)
}

# make forest plots
make_figure_3d_type1 <- function(){
  log_3d <- rma(yi=log_p_1,
                sei=se_log_p_1,
                method=tau_estim_method,
                slab=author,
                dat=dat_3d_allIPV)
  
  forest(log_3d,
         slab=author_year,
         transf=function(x) exp(x),
         refline=NA,
         xlab="3-dose type 1 seroconversion",
         ilab=data.frame(country,type,ifelse(dat_3d_allIPV$study_id>100,"*","")),
         ilab.xpos=c(-0.25,0.25,-0.82), 
         ilab.lab=c("Country", "Vaccine", ""),
         ilab.pos=4,
         cex=.6,
         addcred=F,
         ylim=c(0,47),
         xlim=c(-0.8,1.5),
         rows=gen_rma_rows(unlist(all_schedule_3d_allIPV$count), 
                           starting_row=1, gap_rows=3), 
         alim=c(0.5, 1), 
         order=1:nrow(log_3d$data)
  )
  
  ## add text for IPV schedule
  text(-0.8, gen_rma_schedule_rows(c_vec=unlist(all_schedule_3d_allIPV$count),gap_rows=3) +0.2, 
       unlist(all_schedule_3d_allIPV$schedule), cex=.6, pos=4, font=2)

  
  dat_3d_allIPV_rma3 <- rma(yi=log_p_1,
                            sei=se_log_p_1,
                            method=tau_estim_method,
                            dat=dat_3d_allIPV,
                            subset=schedule=='6, 10, 14 weeks')
  
  dat_3d_allIPV_rma4 <- rma(yi=log_p_1,
                            sei=se_log_p_1,
                            method=tau_estim_method,
                            dat=dat_3d_allIPV,
                            subset=schedule=='2, 3, 4 months')
  
  addpoly(dat_3d_allIPV_rma4, row=8, cex=.6, transf=function(x) exp(x),
          mlab="Average seroconversion")
  addpoly(dat_3d_allIPV_rma3, row=25, cex=.6, transf=function(x) exp(x),
          mlab="Average seroconversion")

  
}

make_figure_3d_type2 <- function(){
  log_3d <- rma(yi=log_p_2,
                sei=se_log_p_2,
                method=tau_estim_method,
                slab=author,
                dat=dat_3d_allIPV)
  
  forest(log_3d,
         slab=author,
         transf=function(x) exp(x),
         refline=NA,
         xlab="3-dose type 2 seroconversion",
         ilab=data.frame(country,type,ifelse(dat_3d_allIPV$study_id>100,"*","")),
         ilab.xpos=c(-0.25,0.25,-0.82), 
         ilab.lab=c("Country", "Vaccine",""),
         ilab.pos=4,
         cex=.6,
         addcred=F,
         ylim=c(0,47),
         xlim=c(-0.8,1.5),
         rows=gen_rma_rows(unlist(all_schedule_3d_allIPV$count), 
                           starting_row=1, gap_rows=3), 
         alim=c(0.5, 1), 
         order=1:nrow(log_3d$data)
  )
  
  ## add text for IPV schedule
  text(-0.8, gen_rma_schedule_rows(c_vec=unlist(all_schedule_3d_allIPV$count),gap_rows=3) +0.2, unlist(all_schedule_3d_allIPV$schedule), cex=.6, pos=4, font=2)
  
  dat_3d_allIPV_rma3 <- rma(yi=log_p_2,
                            sei=se_log_p_2,
                            method=tau_estim_method,
                            dat=dat_3d_allIPV,
                            subset=schedule=='6, 10, 14 weeks')
  
  dat_3d_allIPV_rma4 <- rma(yi=log_p_2,
                            sei=se_log_p_2,
                            method=tau_estim_method,
                            dat=dat_3d_allIPV,
                            subset=schedule=='2, 3, 4 months')
  
  addpoly(dat_3d_allIPV_rma4, row=8, cex=.6, transf=function(x) exp(x),
          mlab="Average seroconversion")
  addpoly(dat_3d_allIPV_rma3, row=25, cex=.6, transf=function(x) exp(x),
          mlab="Average seroconversion")
  
}

make_figure_3d_type3 <- function(){
  log_3d <- rma(yi=log_p_3,
                sei=se_log_p_3,
                method=tau_estim_method,
                slab=author,
                dat=dat_3d_allIPV)
  
  forest(log_3d,
         slab=author,
         transf=function(x) exp(x),
         refline=NA,
         xlab="3-dose type 3 seroconversion",
         ilab=data.frame(country,type,ifelse(dat_3d_allIPV$study_id>100,"*","")),
         ilab.xpos=c(-0.25,0.25,-0.82), 
         ilab.lab=c("Country", "Vaccine",""),
         ilab.pos=4,
         cex=.6,
         addcred=F,
         ylim=c(0,47),
         xlim=c(-0.8,1.5),
         rows=gen_rma_rows(unlist(all_schedule_3d_allIPV$count), 
                           starting_row=1, gap_rows=3), 
         alim=c(0.5, 1), 
         order=1:nrow(log_3d$data)
  )
  
  ## add text for IPV schedule
  text(-0.8, gen_rma_schedule_rows(c_vec=unlist(all_schedule_3d_allIPV$count),gap_rows=3) +0.2, 
       unlist(all_schedule_3d_allIPV$schedule), cex=.6, pos=4, font=2)
  
  dat_3d_allIPV_rma3 <- rma(yi=log_p_3,
                            sei=se_log_p_3,
                            method=tau_estim_method,
                            dat=dat_3d_allIPV,
                            subset=schedule=='6, 10, 14 weeks')
  
  dat_3d_allIPV_rma4 <- rma(yi=log_p_3,
                            sei=se_log_p_3,
                            method=tau_estim_method,
                            dat=dat_3d_allIPV,
                            subset=schedule=='2, 3, 4 months')
  
  addpoly(dat_3d_allIPV_rma4, row=8, cex=.6, transf=function(x) exp(x),
          mlab="Average seroconversion")
  addpoly(dat_3d_allIPV_rma3, row=25, cex=.6, transf=function(x) exp(x),
          mlab="Average seroconversion")
  
}



### all 4-dose, IPV or sIPV
dat_34d %>% 
  filter(num_doses == 4) %>% 
  with(table(schedule, type))

dat_4d_allIPV <- dat_34d %>% filter(num_doses == 4)

dat_4d_allIPV %>% group_by(schedule) %>% summarise(count=n()) -> all_schedule_4d_allIPV

make_figure_4d_type1 <- function(){
  log_4d <- rma(yi=log_p_1,
                sei=se_log_p_1,
                method=tau_estim_method,
                slab=author,
                dat=dat_4d_allIPV)
  
  forest(log_4d,
         slab=author_year,
         transf=function(x) exp(x),
         refline=NA,
         xlab="4-dose type 1 seroconversion",
         ilab=data.frame(country,type,ifelse(dat_4d_allIPV$study_id>100,"*","")),
         ilab.xpos=c(-0.25,0.25,-0.82), 
         ilab.lab=c("Country", "Vaccine",""),
         ilab.pos=4,
         cex=.6,
         addcred=F,
         ylim=c(0,10),
         xlim=c(-1,2),
         rows=gen_rma_rows(unlist(all_schedule_4d_allIPV$count), 
                           starting_row=2, gap_rows=2), 
         alim=c(0.5,1), 
         order=1:nrow(log_4d$data)
  )
  
  ## add text for IPV schedule
  text(-1, gen_rma_schedule_rows(c_vec=unlist(all_schedule_4d_allIPV$count),
                                 starting_row=2,gap_rows=2) +0.2, 
       unlist(all_schedule_4d_allIPV$schedule), cex=.6, pos=4, font=2)
  
  dat_4d_allIPV_rma1 <- rma(yi=log_p_1,
                            sei=se_log_p_1,
                            method=tau_estim_method,
                            dat=dat_4d_allIPV,
                            subset=schedule=='2, 3, 4, 18 months')
  
  addpoly(dat_4d_allIPV_rma1, row=1, cex=.6, transf=function(x) exp(x),
          mlab="Average seroconversion")
}

make_figure_4d_type2 <- function(){
  log_4d <- rma(yi=log_p_2,
                sei=se_log_p_2,
                method=tau_estim_method,
                slab=author,
                dat=dat_4d_allIPV)
  
  forest(log_4d,
         slab=author,
         transf=function(x) exp(x),
         refline=NA,
         xlab="4-dose type 2 seroconversion",
         ilab=data.frame(country,type,ifelse(dat_4d_allIPV$study_id>100,"*","")),
         ilab.xpos=c(-0.25,0.25,-0.82), 
         ilab.lab=c("Country", "Vaccine",""),
         ilab.pos=4,
         cex=.6,
         addcred=F,
         ylim=c(0,10),
         xlim=c(-1,2),
         rows=gen_rma_rows(unlist(all_schedule_4d_allIPV$count), 
                           starting_row=2, gap_rows=2), 
         alim=c(0.5,1), 
         order=1:nrow(log_4d$data)
  )
  
  ## add text for IPV schedule
  text(-1, gen_rma_schedule_rows(c_vec=unlist(all_schedule_4d_allIPV$count),
                                 starting_row=2,gap_rows=2) +0.2, 
       unlist(all_schedule_4d_allIPV$schedule), cex=.6, pos=4, font=2)
  
  dat_4d_allIPV_rma1 <- rma(yi=log_p_2,
                            sei=se_log_p_2,
                            method=tau_estim_method,
                            dat=dat_4d_allIPV,
                            subset=schedule=='2, 3, 4, 18 months')
  
  addpoly(dat_4d_allIPV_rma1, row=1, cex=.6, transf=function(x) exp(x),
          mlab="Average seroconversion")
}

make_figure_4d_type3 <- function(){
  log_4d <- rma(yi=log_p_3,
                sei=se_log_p_3,
                method=tau_estim_method,
                slab=author,
                dat=dat_4d_allIPV)
  
  forest(log_4d,
         slab=author,
         transf=function(x) exp(x),
         refline=NA,
         xlab="4-dose type 3 seroconversion",
         ilab=data.frame(country,type,ifelse(dat_4d_allIPV$study_id>100,"*","")),
         ilab.xpos=c(-0.25,0.25,-0.82), 
         ilab.lab=c("Country", "Vaccine",""),
         ilab.pos=4,
         cex=.6,
         addcred=F,
         ylim=c(0,10),
         xlim=c(-1,2),
         rows=gen_rma_rows(unlist(all_schedule_4d_allIPV$count), 
                           starting_row=2, gap_rows=2), 
         alim=c(0.5,1), 
         order=1:nrow(log_4d$data)
  )
  
  ## add text for IPV schedule
  text(-1, gen_rma_schedule_rows(c_vec=unlist(all_schedule_4d_allIPV$count),
                                 starting_row=2,gap_rows=2) +0.2, 
       unlist(all_schedule_4d_allIPV$schedule), cex=.6, pos=4, font=2)
  
  dat_4d_allIPV_rma1 <- rma(yi=log_p_3,
                            sei=se_log_p_3,
                            method=tau_estim_method,
                            dat=dat_4d_allIPV,
                            subset=schedule=='2, 3, 4, 18 months')
  
  addpoly(dat_4d_allIPV_rma1, row=1, cex=.6, transf=function(x) exp(x),
          mlab="Average seroconversion")
}


### other analysis that explores impact of OPV vaccination history

rma(yi=log_p_1,
    sei=se_log_p_1,
    method=tau_estim_method,
    dat=dat_3d_allIPV,
    subset=schedule=='6, 10, 14 weeks',
    mods= ~ birthdose_OPV)


rma(yi=log_p_1,
    sei=se_log_p_1,
    method=tau_estim_method,
    dat=dat_3d_allIPV,
    subset=schedule=='2, 3, 4 months',
    mods= ~ birthdose_OPV)


png(file="figure/dose3serotype1.png", width = 2000, height = 2000, res=300)
make_figure_3d_type1()
dev.off()


png(file="figure/dose3serotype2.png", width = 2000, height = 2000, res=300)
make_figure_3d_type2()
dev.off()


png(file="figure/dose3serotype3.png", width = 2000, height = 2000, res=300)
make_figure_3d_type3()
dev.off()



png(file="figure/dose4serotype1.png", width = 1500, height = 1000,  res=300)
make_figure_4d_type1()
dev.off()


png(file="figure/dose4serotype2.png", width = 1500, height = 1000,  res=300)
make_figure_4d_type2()
dev.off()


png(file="figure/dose4serotype3.png", width = 1500, height = 1000,  res=300)
make_figure_4d_type3()
dev.off()