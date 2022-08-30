# compare AS expression level to non AS expression ------------------------

# High confidence expression per project ----------------------------------------------

# Load required packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(ggupset)


# Reading data files

files = list.files("./Sig_express_030621/")

df_plot = data.frame()

df_summary = data.frame()

# Warning! The loop may take a couple of hours to wrap up
for (i in files) {
  # Load expression data
  exp_data = read.delim(paste("./Sig_express_030621/", i, sep = ""))
  
  # Get probe strand orientation
  strand = as.data.frame(str_split_fixed(exp_data$probe_name, "_", 4))
  strand = as.character(strand$V2)
  
  # Get transcript name and functional category
  transcript_name = as.character(exp_data$transcript_name)
  category = as.character(exp_data$category)
  probe_name = as.character(exp_data$probe_name)
  
  # Get clean exp matrix
  exp_data_1 = exp_data[,-c(1:6)]
  row.names(exp_data_1) = exp_data$probe_name
  exp_data_1 = cbind(transcript_name, probe_name, strand, category, exp_data_1)
  
  # Filter genes with no expression evidence - Sanity check!
  thresh = !is.na(exp_data_1[,c(5:ncol(exp_data_1))])
  keep = rowSums(thresh) >= 1 # Expression evidence in at least one sample
  exp_data_2 = exp_data_1[keep,]
  row.names(exp_data_2) = exp_data_2$probe_name
  
  # Check if every SAS has both AS and SS
  check_SS_AS = as.data.frame(table(exp_data_2$transcript_name))
  check_SS_AS_ok = subset(check_SS_AS, Freq == 2)
  check_SS_AS_ok = as.character(check_SS_AS_ok$Var1)
  exp_data_3 = exp_data_2[ exp_data_2$transcript_name %in% check_SS_AS_ok , ]
  
  # Check if every SAS has both AS and SS
  x = as.data.frame(table(exp_data_3$transcript_name))
  y = subset(x, Freq == 2)
  nrow(x) == nrow(y) #if TRUE, it's OK!
  if(nrow(x) == nrow(y)) {
    print(paste("Experiment", i, "is okay!", sep = " "))
  } else {
    print(paste("ATENTION: Experiment", i, "is NOT okay! Please double-check.", sep = " "))
  }
  
  # Summary number of SAS
  df_summary_1 = as.data.frame(length(unique(exp_data_3$transcript_name)))
  names(df_summary_1) = "N_SAS"
  df_summary_1$N_SEG = nrow(exp_data)
  df_summary_1$project = str_sub(i, start=1, end=-9)
  df_summary_1$samples = ncol(exp_data_3[,-c(1:4)])
  
  df_summary = rbind(df_summary, df_summary_1)
  
  # Prepare data frame for plots
  
  exp_data_4 = gather(exp_data_3,
                      sample, exp,
                      5:ncol(exp_data_3),
                      factor_key=TRUE)
  
  exp_data_4$project = rep(str_sub(i, start=1, end=-9), nrow(exp_data_4))
  
  df_plot = rbind(df_plot, exp_data_4)
  
  # Rename SAS strand based on AS/SS or just SS expression
  
  df_plot_2 = data.frame()
  
  SAS = as.character(unique(df_plot$transcript_name))
  
  for (l in SAS) {
    a = subset(df_plot, transcript_name == l)
    sample = as.character(unique(a$sample))
    
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

# Select leaf samples -----------------------------------------------------

sample_annot = read.delim("sample_annot.txt", header = T)
sample_annot = sample_annot[,-2]

df_plot_3 = inner_join(df_plot_2, sample_annot, by = "sample")

df_plot_4 = subset(df_plot_3, tissue == "L1")

# Plot codes --------------------------------------------------------------

# Load required packages

library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

# Rename dataframe columns

df_plot_4[df_plot_4 == "ancestral"] = "Ancestral"
df_plot_4[df_plot_4 == "circadian_RB855453"] = "Circadian I"
df_plot_4[df_plot_4 == "circadian_SP803280"] = "Circadian II"
df_plot_4[df_plot_4 == "drought_field"] = "Drought I"
df_plot_4[df_plot_4 == "drought_SP901638"] = "Drought II"
df_plot_4[df_plot_4 == "drought_SP803280"] = "Drought III"
df_plot_4[df_plot_4 == "ethylene"] = "Ethylene"
df_plot_4[df_plot_4 == "growth_maturation"] = "G&M"

df_plot_1 = df_plot_4
df_plot_1[df_plot_1 == "SS_only"] = "SS"

# Plot 1 (Figure 5a)

p1 = ggplot(df_plot_1, aes(x=strand, y=exp, fill=strand)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  scale_fill_manual(values=c("#F8766D", "#4dac26")) +
  geom_violin(trim=FALSE, alpha = 0.3) +
  geom_boxplot(width=0.1, outlier.alpha = 0.2, outlier.size = 0.5) +
  facet_wrap(~ project) +
  labs(x = "Probe strand",
       y = "log2(Intensity)") +
  guides(fill="none") +
  theme_bw()
p1

# T-test ------------------------------------------------------------------

df_mean_test = data.frame()

projects = as.character(unique(df_plot_4$project))

# Warning! The loop may take a while to wrap up
for (i in projects) {
  
  x = subset(df_plot_4, project == i)
  SAS_list = as.character(unique(x$transcript_name))
  
  for (j in SAS_list) {
    
    df = df_plot_4 %>%
      filter(project == i & transcript_name == j) %>%
      filter(strand != "AS") %>%
      select(transcript_name, probe_name, strand, sample, exp, project, sample_2, condition, tissue)
    
    SS = subset(df, strand == "SS")
    SS = na.omit(SS)
    SS_only = subset(df, strand == "SS_only")
    SS_only = na.omit(SS_only)
    
    # Run test only if we have expression evidence of at least two samples per SAS (SS and SS_only)
    # Treatment
    if (length(unique(factor(df$strand))) == 2 &&
        nrow(SS) >= 2 &&
        nrow(SS_only) >= 2) {
      
      res = t.test(exp ~ strand, data = df)
      pval = unlist(res$p.value)
      means = unlist(res$estimate)
      
      temp_df = data.frame(transcript_name = j,
                           project = i,
                           mean_SS = means[1],
                           mean_SS_only = means[2],
                           pval = pval)
      
      df_mean_test = rbind(df_mean_test, temp_df)
      
    }
    
    else {
      
      temp_df = data.frame(transcript_name = j,
                           project = i,
                           mean_SS = NA,
                           mean_SS_only = NA,
                           pval = NA)
      
      df_mean_test = rbind(df_mean_test, temp_df)
      
    }
    
  }
  
}

df_mean_test = na.omit(df_mean_test)
df_mean_test


# FDR correction (alpha 0.05)-----------------------------------------------
Sig = sapply(p.adjust(df_mean_test$pval, "fdr"), function(x) x <= 0.05)
Sig[Sig == TRUE] = "*"
Sig[Sig == FALSE] = "ns"
df_mean_test_fdr = cbind(df_mean_test, Sig)
df_mean_test_fdr$transcript_name_project = paste(df_mean_test_fdr$transcript_name, df_mean_test_fdr$project, sep = "_")

df_mean_test_final = data.frame()

for (i in unique(df_mean_test_fdr$transcript_name_project)) {
  x = subset(df_mean_test_fdr, transcript_name_project == i )
  if(x$Sig == "ns") {
    x$mean = "Neutral"
  }
  if(x$Sig == "*" && x$mean_SS > x$mean_SS_only) {
    x$mean = "Concordant"
  }
  if(x$Sig == "*" && x$mean_SS < x$mean_SS_only) {
    x$mean = "Discordant"
  }
  
  df_mean_test_final = rbind(df_mean_test_final,x)
}


# Specific examples shown in Figure 5b  -------------------------------------------------------

df_plot_4$transcript_name_project = paste(df_plot_4$transcript_name, df_plot_4$project, sep = "_")

# Discordant expression (SS < SS_only)
df1 = df_plot_4 %>%
  filter(transcript_name_project == "SCBFRZ2017D06.g_Ethylene")

# Neutral expression (ns)
df2 = df_plot_4 %>%
  filter(transcript_name_project == "SCRFLB1054D01.g_Ancestral")

# Concordant expression (SS > SS_only)
df3 = df_plot_4 %>%
  filter(transcript_name_project == "SCSBAD1052H02.g_Circadian I")


# Merge examples data frames

df4 = rbind(df1, df2, df3)

# Plot 2 (Figure 5b)
p2 = ggplot(df4, aes(x=strand, y=exp, fill=strand)) + 
  scale_fill_manual(values=c("#F8766D", "#ffc425", "#00BFC4")) +
  geom_violin(trim=FALSE, alpha = 0.3) +
  geom_boxplot(width=0.1, outlier.alpha = 0.2, outlier.size = 0.5) +
  #geom_boxplot() +
  facet_wrap(~ project + transcript_name, scales = "free", ncol = 1) +
  labs(x = "Probe strand",
       y = "log2(Intensity)") +
  guides(fill="none") +
  theme_bw()
p2

# Final plots --------------------------------------------------------------

# Figure 5a and 5b
pdf(file = "Figure_5.pdf", width = 10, height = 7)
ggarrange(p1, p2,widths = c(2.8,1.5),
          labels = c("(a)", "(b)"))
dev.off()

# Figure 5c (Upset plot) --------------------------------------------------------------

# Load required package
library(UpSetR)

data = df_mean_test_final[,c(1,2,8)]
names(data) = c("SAS", "Project", "Expression")
row.names(data) = c(1:nrow(data))

# Preparing data
# Making a list with all transcripts grouped by sense/antisense expression pattern
y = list(
  Concordant = data[which(data[,3]=="Concordant"),1],
  Discordant = data[which(data[,3]=="Discordant"),1],
  Neutral = data[which(data[,3]=="Neutral"),1])
y

# Figure 5c in pdf format

pdf(file = "Figure_5c.pdf", width = 10, height = 7)

upset(fromList(y), nsets=3, keep.order="True",

      queries = list(list(query = intersects,
                          params = list ("Concordant"), color = "darkblue", active = F),
                     list(query = intersects,
                          params = list ("Discordant"), color = "darkred", active = F),
                     list(query = intersects,
                          params = list ("Neutral"), color = "khaki2", active = F),
                     list(query = intersects,
                          params = list ("Neutral", "Concordant"), color = "lightblue2", active = F),
                     list(query = intersects,
                          params = list ("Neutral", "Discordant"), color = "pink1", active = F)),
      main.bar.color = "gray40", mainbar.y.label = "Number of SAS", 
      sets.bar.color = c("khaki2", "darkblue", "darkred"), 
      sets.x.label = "Set Size", point.size = 4.8, line.size = 2.0,
      mb.ratio = c(0.70,0.30), show.numbers = "yes", shade.color = "white", 
      matrix.dot.alpha =0.8, scale.intersections = "identity",
      scale.sets = "identity", text.scale = 2.0, set_size.show = F)

dev.off()