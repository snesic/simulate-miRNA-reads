library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)

## Read results from scoring_script:  
data <- read_tsv(args[1]) # PATH TO FILE
colnames(data) <- c('name','aligned','tool')

data = data %>% separate(tool, c("tool", "simulation"), ";")  #%>%
data$simulation = gsub("test_","",data$simulation, ignore.case = T)

data = data %>% mutate(aligned=case_when(is.na(.$aligned) ~ "NA",
                                         .$aligned=="yes" ~ "Yes",
                                         .$aligned=="no"  ~ "No",
                                         .$aligned=="multi-no"  ~ "multi-no",
                                         .$aligned=="multi-yes"  ~ "multi-yes"))

## Subsets from data

dt_gainRand_3p = subset(data,data['simulation']=='gainRand_3p')
dt_gainRand_5p = subset(data,data['simulation']=='gainRand_5p')

dt_lossRand_3p = subset(data,data['simulation']=='lossRand_3p')
dt_lossRand_5p = subset(data,data['simulation']=='lossRand_5p')

dt_gainTemp_3p = subset(data,data['simulation']=='gainTemp_3p')
dt_gainTemp_5p = subset(data,data['simulation']=='gainTemp_5p')

## Plot Graph

tri_graph <- function(data){
    
    str_name = gsub("dt_","",deparse(substitute(data)), ignore.case = T)
    if (str_name == 'data'){
        str_name = "all"
    }
    
    dt = data %>% group_by(tool,aligned) %>% summarise(total = n()) %>% as_data_frame()
    #
    dt_perc = group_by(dt,tool) %>% mutate(percent = round(100*(total/sum(total)),3))
    dt_perc$aligned <- sub(pattern="No", "2. Aligned Incorectly", dt_perc$aligned)
    dt_perc$aligned <- sub(pattern="NA", '1. Not Mapped', dt_perc$aligned)
    dt_perc$aligned <- sub(pattern="Yes", "5. Aligned Correctly", dt_perc$aligned)
    dt_perc$aligned <- sub(pattern="multi-yes", "3. Multi-aligned true", dt_perc$aligned)
    dt_perc$aligned <- sub(pattern="multi-no", "4. Multi-aligned false", dt_perc$aligned)
    #
    ggplot(dt_perc, aes(x = tool, y = percent, fill = aligned, label = percent)) +
      geom_bar(stat = "identity") +
      geom_text(size = 4, position = position_stack(vjust = 0.5),color='navyblue',check_overlap = TRUE) +
      scale_fill_manual(values=c("5. Aligned Correctly" = "darkgreen", "1. Not Mapped" = "firebrick4", "2. Aligned Incorectly" = "darkorange1","3. Multi-aligned true"='yellow', '4. Multi-aligned false'='peru') ) + ggtitle("Mapping statistics : ", str_name)
}

#
tri_graph(data)
#
tri_graph(dt_gainRand_3p)

tri_graph(dt_gainTemp_3p)
tri_graph(dt_gainTemp_5p)

tri_graph(dt_lossRand_3p)
tri_graph(dt_lossRand_5p)
