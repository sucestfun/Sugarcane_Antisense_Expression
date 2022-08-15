# compare AS expression level to non AS expression ------------------------

# High confidence expression per project ----------------------------------------------

library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

files = list.files("./Sig_express_030621/")

df_plot = data.frame()

df_summary = data.frame()

HC_SAS = vector(mode = "list", length = length(files))
names(HC_SAS) = files

# Warning! The loop may take a couple of minutes to wrap up
for (i in files) {
  #Load expression data
  exp_data = read.delim(paste("./Sig_express_030621/", i, sep = ""))
  
  #Get probe strand orientation
  strand = as.data.frame(str_split_fixed(exp_data$probe_name, "_", 4))
  strand = as.character(strand$V2)
  
  #Get transcript name and functional category
  transcript_name = as.character(exp_data$transcript_name)
  category = as.character(exp_data$category)
  probe_name = as.character(exp_data$probe_name)
  
  #Get clean exp matrix
  exp_data_1 = exp_data[,-c(1:6)]
  row.names(exp_data_1) = exp_data$probe_name
  exp_data_1 = cbind(transcript_name, probe_name, strand, category, exp_data_1)
  
  exp_data_2 = as.data.frame(exp_data_1[,2])
  names(exp_data_2) = "probe_name"
  
  #Keep genes with expression evidence on both replicates
  for (j in seq(5, ncol(exp_data_1), by = 2)) {
    x = exp_data_1[,j:(j+1)]
    thresh = !is.na(x)
    keep = rowSums(thresh) == 2
    counts.keep = x[keep,]
    probe_name = row.names(counts.keep)
    counts.keep = cbind(probe_name, counts.keep)
    exp_data_2 = full_join(exp_data_2, counts.keep, by = "probe_name", all=TRUE)
  }
  
  #Filter genes with no expression evidence
  thresh <- !is.na(exp_data_2)
  keep <- rowSums(thresh) >= 3
  exp_data_3 <- exp_data_2[keep,]
  row.names(exp_data_3) = exp_data_3$probe_name
  
  #Get annotation back
  annot = exp_data_1[,c(1:4)]
  annot2 = annot[ (annot$probe_name) %in% as.character(exp_data_3$probe_name), ]
  exp_data_4 = cbind(annot2, exp_data_3)
  row.names(exp_data_4) = exp_data_4$probe_name
  
  AS = subset(exp_data_4, strand  == "AS")
  
  SS = subset(exp_data_4, strand  == "SS")
  row.names(SS) = SS$transcript_name
  SS_2 = SS[as.character(AS$transcript_name),]
  SS_2 = na.omit(SS_2)
  
  exp_data_5 = rbind(AS, SS_2)
  row.names(exp_data_5) = exp_data_5$probe_name
  
  #Check if every SAS has both AS and SS
  x = as.data.frame(table(exp_data_5$transcript_name))
  y = subset(x, Freq == 2)
  nrow(x) == nrow(y) #if FALSE, double check
  
  #Keep only SAS with both AS and SS
  exp_data_6 = exp_data_5[ (exp_data_5$transcript_name) %in% as.character(y$Var1), ]
  
  #Double check if every SAS has both AS and SS
  x = as.data.frame(table(exp_data_6$transcript_name))
  y = subset(x, Freq == 2)
  nrow(x) == nrow(y) #if TRUE, OK
  
  #Keep only SAS with expression evidence for both AS and SS
  exp_data_7 = exp_data_6[order(exp_data_6$transcript_name),] 
  exp_data_7 = exp_data_7[,-c(1:5)]
  exp_data_7 = as.data.frame(t(exp_data_7))
  sample_name = as.character(row.names(exp_data_7))
  exp_data_8 = cbind(sample_name, exp_data_7)
  
  exp_data_9 = as.data.frame(sample_name)
  row.names(exp_data_9) = sample_name
  
  for (k in seq(2, nrow(exp_data_8), by = 2)) {
    x = exp_data_8[,k:(k+1)]
    thresh = !is.na(x)
    keep = rowSums(thresh) == 2
    counts.keep = x[keep,]
    sample_name = row.names(counts.keep)
    counts.keep = cbind(sample_name, counts.keep)
    exp_data_9 = full_join(exp_data_9, counts.keep, by = "sample_name", all=TRUE)
  }
  
  exp_data_10 = exp_data_9[,-1]
  row.names(exp_data_9) = exp_data_9[,1]
  exp_data_10 = exp_data_9[,-1]
  exp_data_10 = as.data.frame(t(exp_data_10))
  
  HC_probe = as.character(row.names(exp_data_10))
  
  exp_data_10 = exp_data_6[HC_probe,]
  exp_data_11 = exp_data_10[,-5]
  
  #Summary number of SAS
  df_summary_1 = as.data.frame(length(unique(exp_data_11$transcript_name)))
  names(df_summary_1) = "N_SAS"
  df_summary_1$N_SEG = nrow(exp_data)
  df_summary_1$project = str_sub(i, start=1, end=-9)
  df_summary_1$samples = ncol(exp_data_1[,-c(1:4)])
  
  df_summary = rbind(df_summary, df_summary_1)
  
  HC_SAS[[i]] = as.character(unique(exp_data_11$transcript_name))
  
  #Prepare data frame for plots
  
  exp_data_12 = gather(exp_data_11,
                       sample, exp,
                       5:ncol(exp_data_11),
                       factor_key=TRUE)
  
  exp_data_12$project = rep(str_sub(i, start=1, end=-9), nrow(exp_data_12))
  
  df_plot = rbind(df_plot, exp_data_12)
  
  #Rename SAS strand based on AS/SS or just SS expression
  
  df_plot_2 = data.frame()
  
  SAS = as.character(unique(df_plot$transcript_name))
  
  for (l in SAS) {
    a = subset(df_plot, transcript_name == l)
    sample = as.character(a$sample)
    
    for (m in sample) {
      x = subset(a, sample == m)
      if (is.na(x[1,6])) {
        x[2,3] = "SS_only"
      }
      df_plot_2 = rbind(df_plot_2, x)
    }
  }
}

save(df_plot_2, file = "df_plot_2.RData")
save(df_summary, file = "df_summary.RData")
write.table(df_summary, file = "df_summary_general.txt")


# Plot codes --------------------------------------------------------------

library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

load("df_plot_2.RData")

df_plot_2[df_plot_2 == "ancestral"] = "Ancestral"
df_plot_2[df_plot_2 == "circadian_RB855453"] = "Circadian I"
df_plot_2[df_plot_2 == "circadian_SP803280"] = "Circadian II"
df_plot_2[df_plot_2 == "drought_field"] = "Drought I"
df_plot_2[df_plot_2 == "drought_SP901638"] = "Drought II"
df_plot_2[df_plot_2 == "drought_SP803280"] = "Drought III"
df_plot_2[df_plot_2 == "ethylene"] = "Ethylene"
df_plot_2[df_plot_2 == "growth_maturation"] = "G&M"

df_plot_1 = df_plot_2
df_plot_1[df_plot_1 == "SS_only"] = "SS"


p1 = ggplot(df_plot_1, aes(x=strand, y=exp, fill=strand)) + 
  scale_fill_manual(values=c("#F8766D", "#4dac26")) +
  geom_violin(trim=FALSE, alpha = 0.3) +
  geom_boxplot(width=0.1, outlier.alpha = 0.2) +
  facet_wrap(~ project, scales = "free") +
  labs(x = "Probe stand",
       y = "log2(Intensity)") +
  guides(fill=FALSE) +
  theme_bw()


# T-test ------------------------------------------------------------------

df_mean_test = data.frame()

projects = as.character(unique(df_plot_2$project))

for (i in projects) {
  
  x = subset(df_plot_2, project == i)
  SAS_list = as.character(unique(x$transcript_name))
  
  for (j in SAS_list) {
    
    df = df_plot_2 %>%
      filter(project == i & transcript_name == j) %>%
      filter(strand != "AS") %>%
      select(transcript_name, strand, sample, exp, project)
    
    if (length(unique(factor(df$strand))) == 2) {
      
      res = t.test(exp ~ strand, data = df)
      pval = unlist(res$p.value)
      means = unlist(res$estimate)
      
      temp_df = df[1,c(1,5)]
      
      if (means[1] > means[2]) {
        temp_df$mean = "SS>SS_only" 
      }
      
      if (means[1] < means[2]) {
        temp_df$mean = "SS<SS_only" 
      }
      
      if (pval >= 0.05) {
        temp_df$pval = "ns" 
      }
      
      if (pval <= 0.05) {
        temp_df$pval = "*" 
      }
      
      if (pval <= 0.01) {
        temp_df$pval = "**" 
      }
      
      if (pval <= 0.001) {
        temp_df$pval = "***" 
      }
      
      df_mean_test = rbind(df_mean_test, temp_df)
      
    }
    
    else {
      
      temp_df = df[1,c(1,5)]
      temp_df$pval = NA
      temp_df$mean = NA
      df_mean_test = rbind(df_mean_test, temp_df)
      
    }
    
  }
  
}


# Summary table -----------------------------------------------------------

df_summary = data.frame()

for (i in projects) {
  df_summary_temp_1 = subset(df_mean_test, project == i)
  df_summary_temp_1 = na.omit(df_summary_temp_1)
  df_summary_temp_1 = subset(df_summary_temp_1, pval != "ns")
  
  if (nrow(df_summary_temp_1) != 0) {
    df_1_temp = as.data.frame(table(df_summary_temp_1$mean))
    df_1_temp$project = i
    
    df_summary = rbind(df_summary, df_1_temp)
  }
  
}

write.table(df_mean_test, file = "df_mean_test.txt")
write.table(df_summary, file = "df_summary.txt")

#summary_table = read.delim("summary_table.txt", header = T)
#names(summary_table) = c("Project", "N_SEG", "Both_rep", "N_SAS", "SS<SS_only", "SS>SS_only", "No_SS_only")

# Specific examples -------------------------------------------------------

# No difference
df1 = df_plot_2 %>%
  filter(project == "Circadian I" & transcript_name == "SCACCL6007G08.g")

# SS < SS_only
df2 = df_plot_2 %>%
  filter(project == "Ethylene" & transcript_name == "SCACAD1038A12.g")

# SS > SS_only
df3 = df_plot_2 %>%
  filter(project == "Drought I" & transcript_name == "SCACLR1036F06.g")

#Merge examples data frames

df4 = rbind(df1, df2, df3)

p2 = ggplot(df4, aes(x=strand, y=exp, fill=strand)) + 
  scale_fill_manual(values=c("#F8766D", "#ffc425", "#00BFC4")) +
  geom_violin(trim=FALSE, alpha = 0.3) +
  geom_boxplot(width=0.1, outlier.alpha = 0.2) +
  facet_wrap(~ project + transcript_name, scales = "free", ncol = 1) +
  labs(x = "Probe stand",
       y = "log2(Intensity)") +
  guides(fill=FALSE) +
  theme_bw()


# Final plot --------------------------------------------------------------

pdf(file = "Figure_6.pdf", width = 10, height = 7)
ggarrange(p1, p2,widths = c(2.8,1.5),
          labels = c("(a)", "(b)"))
dev.off()



