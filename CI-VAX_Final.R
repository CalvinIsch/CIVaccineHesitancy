library(ggplot2); library(dplyr); library(ggpubr); library(igraph)
library(tidyverse); library(DescTools); library(lmtest); library(sandwich)
library(ggpattern); library(sf); library(viridis); library(maps); library(usmap)

################################## Loading data ###################################
df = read.csv('data/all_data_2023.09.27.csv')
df_confeds <- df %>% filter(confederate == 'yes')

df$err_1_abs <- abs(df$err_1)
df$err_2_abs <- abs(df$err_2)
df$err_3_abs <- abs(df$err_3)
df$change_abs_err <- df$err_3_abs- df$err_1_abs
df$mag_revision <- abs(df$value_3 - df$value_1)

# Remove confederates
df <- df %>% filter(confederate == 'no')


################################## Figure 1 ###################################
# Create Star Network
star_g <- make_star(40, center=1, mode="undirected")
vertex_colors <- rep("darkgray", vcount(star_g))
vertex_colors[1] <- "red"
vertex_colors_2 <- rep(c('darkgray','navy','darkgreen','brown','black'), vcount(star_g)/5)
vertex_colors_2[1] <- "red"

# Create Egalitarian Network (Regular Graph with each node having 4 edges)
egal <- read.csv('data/neighbors.csv')
egal_g <- make_empty_graph(n = nrow(egal), directed = TRUE)

for(i in 1:nrow(egal)){
  node <- egal$Node[i]
  neighbors <- c(egal$N1[i], egal$N2[i], egal$N3[i], egal$N4[i])
  for(neigh in neighbors){
    egal_g <- add_edges(egal_g, c(node, neigh))
  }
}
egal_g <- as.undirected(egal_g, mode="collapse")
N <- vcount(egal_g)
half_N <- N/2

# Define radii for inner and outer circles
r_inner <- 0.8
r_outer <- 1
angles <- seq(0, 2*pi, length.out = half_N + 1)
angles <- angles[-(half_N + 1)]  
x_inner <- r_inner * cos(angles)
y_inner <- r_inner * sin(angles)
x_outer <- r_outer * cos(angles)
y_outer <- r_outer * sin(angles)
x_coords <- c(x_inner, x_outer)
y_coords <- c(y_inner, y_outer)

label_size = 4
vertex_size = 8
png(filename="figures/Figure_1.png", width=1200, height=800)
par(mfrow=c(2,2),mar=c(0,5,5,0),mgp=c(-2.5, 1.5, 0))
plot(star_g, vertex.label=NA, vertex.size=vertex_size, vertex.color=vertex_colors)
title(main="Centralized", ylab="Homophilous", cex.lab=label_size, cex.main=label_size, font.main=2, font.lab=2, las=1)
plot(egal_g, layout=matrix(c(x_coords, y_coords), ncol=2), vertex.label=NA, vertex.size=vertex_size, vertex.color=vertex_colors)
title(main="Egalitarian", cex.main=label_size, font.main=2)
plot(star_g, vertex.label=NA, vertex.size=vertex_size, vertex.color=vertex_colors_2)
title(ylab="Diverse", cex.lab=label_size, font.lab=2, las=1)
plot(egal_g, layout=matrix(c(x_coords, y_coords), ncol=2), vertex.label=NA, vertex.size=vertex_size, vertex.color=vertex_colors_2)
dev.off()

############################## Figure 2 ##################################
df_table_networks <- df %>% 
  select(network, structure, change_abs_err, composition, value_1, value_3,err_1_abs) %>%
  group_by(structure,network, composition) %>%
  dplyr::summarise(
    mean_r1 = mean(value_1, na.rm = T),
    mean_r3 = mean(value_3, na.rm = T),
    mean_change = mean(change_abs_err, na.rm = T),
    n_people = length(value_1)
  )

# Calculate the crowd estimates in case we choose to look at those.
df_table_networks$Err_1_Net <- abs(df_table_networks$mean_r1 - 35)
df_table_networks$Err_3_Net <- abs(df_table_networks$mean_r3 - 35)
df_table_networks$Mean_Change_Net <- df_table_networks$Err_3_Net - df_table_networks$Err_1_Net

df_table_networks$structure <- factor(df_table_networks$structure, levels=c('c','d'))
levels(df_table_networks$structure) <- c('Centralized','Egalitarian')

df_summary = df_table_networks %>%
  group_by(structure) %>%
  dplyr::summarise(
    grand_mean_ind = mean(mean_change),
    se_ind = sd(mean_change)/sqrt(length(mean_change)),
    cilow_ind=t.test(mean_change)$conf.int[1],
    cihi_ind=t.test(mean_change)$conf.int[2]
  )


ggplot(df_summary, aes(x = structure, y = grand_mean_ind, fill = structure)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_hline(yintercept = 0, color ='black') + 
  geom_errorbar(aes(ymin = cilow_ind,
                    ymax = cihi_ind),
                width = 0.2, position = position_dodge(0.9)) +
  labs(y = "Change in absolute error",
       fill = "Category") +
  scale_fill_manual(values=c("#de2d26", "#3182bd")) + 
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12), 
        axis.title.x=element_blank(),
        panel.background = element_rect(fill='white'),
        plot.background = element_rect(fill='white', color=NA),
        legend.position="right", 
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(size=10), 
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_segment(aes(x = 2.75, xend = 2.75, y = -100, yend = -500),
               lineend = 'round', linejoin = 'round',
               linewidth = 1,
               arrow = arrow(type = "closed", length = unit(.1, "inches"))) +
  geom_text(aes(x = 2.75, y = -550, label='Improves'), size = 5) +
  geom_segment(aes(x = 0.25, xend = 0.25, y = 700, yend =1100),
               lineend = 'round', linejoin = 'round',
               linewidth = 1,
               arrow = arrow(type = "closed", length = unit(.1, "inches"))) +
  geom_text(aes(x = 0.25, y = 1180, label='Worsens'), size = 5)  +
  coord_cartesian(xlim=c(0.4,2.6))
ggsave("figures/Figure-2-Network-types.png",width = 6,height = 3.5,units = "in",dpi = 300)

################################### Differences by Network ######################
df_c <- df_table_networks %>% filter(structure == 'Centralized')
df_e <- df_table_networks %>% filter(structure == 'Egalitarian')

# Crowd Error - SI
print(paste('Centralized: average change in error was ',mean(df_c$Mean_Change_Net), 'which is a',
            100*mean(df_c$Mean_Change_Net)/mean(df_table_networks$Err_1_Net), '% increase.'))
wilcox.test(df_c$Mean_Change_Net)
print(paste('Egalitarian: average change in error was ',mean(df_e$Mean_Change_Net), 'which is a',
            100*mean(df_e$Mean_Change_Net)/mean(df_table_networks$Err_1_Net), '% decrease'))
wilcox.test(df_e$Mean_Change_Net)
wilcox.test(df_e$Mean_Change_Net, df_c$Mean_Change_Net)

# Individual Error - Main text
print(paste('Centralized: average change in error was ',mean(df_c$mean_change), 'which is a',
            100*mean(df_c$mean_change)/mean(df_table_networks$mean_r1), '% increase.'))
wilcox.test(df_c$mean_change)
print(paste('Egalitarian: average change in error was ',mean(df_e$mean_change), 'which is a',
            100*mean(df_e$mean_change)/mean(df_table_networks$mean_r1), '% increase.'))
wilcox.test(df_e$mean_change)

wilcox.test(df_e$mean_change, df_c$mean_change)


############################ Differences by Homophily ##########################
df_m <- df_table_networks %>% filter(composition == 'm')
df_h <- df_table_networks %>% filter(composition == 'aa')

# Crowd Error - SI
wilcox.test(df_m$Mean_Change_Net)
wilcox.test(df_h$Mean_Change_Net)


# Individual Error -Main text
print(paste(mean(df_h$mean_change),mean(df_m$mean_change)))
wilcox.test(df_m$mean_change)
wilcox.test(df_h$mean_change)
wilcox.test(df_m$mean_change, df_h$mean_change)


############################## Differences Both homophily and structure ###
df_ch <- df_table_networks %>% filter(composition == 'aa' & structure == 'Centralized')
df_cm <- df_table_networks %>% filter(composition == 'm' & structure == 'Centralized')
df_eh <- df_table_networks %>% filter(composition == 'aa' & structure == 'Egalitarian')
df_em <- df_table_networks %>% filter(composition == 'm' & structure == 'Egalitarian')

# Crowd Error
print(paste(mean(df_ch$Mean_Change_Net),mean(df_cm$Mean_Change_Net),
            mean(df_eh$Mean_Change_Net),mean(df_em$Mean_Change_Net)))
wilcox.test(df_ch$Mean_Change_Net)
wilcox.test(df_cm$Mean_Change_Net)
wilcox.test(df_eh$Mean_Change_Net)
wilcox.test(df_em$Mean_Change_Net)
wilcox.test(df_cm$Mean_Change_Net, df_em$Mean_Change_Net)
wilcox.test(df_ch$Mean_Change_Net, df_eh$Mean_Change_Net)

# Individual Error - Main text
print(paste('Centralized-Diverse: average change in error was ',mean(df_cm$mean_change), 'which is a',
            100*mean(df_cm$mean_change)/mean(df_table_networks$mean_r1), '% increase.'))
print(paste('Egalitarian-Diverse: average change in error was ',mean(df_em$mean_change), 'which is a',
            100*mean(df_em$mean_change)/mean(df_table_networks$mean_r1), '% increase.'))
wilcox.test(df_cm$mean_change, df_em$mean_change)
print(paste('Centralized-Homo: average change in error was ',mean(df_ch$mean_change), 'which is a',
            100*mean(df_ch$mean_change)/mean(df_table_networks$mean_r1), '% increase.'))
print(paste('Egalitarian-Homo: average change in error was ',mean(df_eh$mean_change), 'which is a',
            100*mean(df_eh$mean_change)/mean(df_table_networks$mean_r1), '% increase.'))
wilcox.test(df_ch$mean_change, df_eh$mean_change)
wilcox.test(df_ch$mean_change)
wilcox.test(df_cm$mean_change)
wilcox.test(df_eh$mean_change)
wilcox.test(df_em$mean_change)


##################################### Figure 3 #######################################
df_table_networks$composition <- factor(df_table_networks$composition, levels=c('aa','m'))
levels(df_table_networks$composition) <- c('Homophilous','Diverse')

df_table_networks$Conditions <- paste(df_table_networks$composition, df_table_networks$structure)
df_table_networks$Conditions <- factor(df_table_networks$Conditions,
                                       levels=c("Homophilous Centralized","Diverse Centralized","Homophilous Egalitarian",
                                                "Diverse Egalitarian"))


df_summary_2 = df_table_networks %>%
  group_by(Conditions) %>%
  dplyr::summarise(
    grand_mean_ind = mean(mean_change),
    se_ind = sd(mean_change)/sqrt(length(mean_change)),
    cilow=t.test(mean_change)$conf.int[1],
    cihi=t.test(mean_change)$conf.int[2]
  )

df_summary_2 <- df_summary_2%>%
  separate(Conditions, into = c("Composition", "Network"), sep = " ", remove=FALSE)
levels(df_summary_2$Conditions) <- c("Homo.\n Cent.","Diverse\n Cent.","Homo.\n Egal.",
                                     "Diverse\n Egal.")

# Your existing code with modifications for patterns
ggplot(df_summary_2, aes(x = Conditions, y = grand_mean_ind, fill = Network, pattern = Composition)) +
  geom_bar_pattern(aes(pattern = Composition), stat = "identity", position = "dodge", color='black', 
                   pattern_density = 0.1, 
                   pattern_spacing = 0.02) + 
  scale_pattern_manual(values = c("stripe", "circle")) +
  geom_hline(yintercept = 0, color = 'black') + 
  geom_errorbar(aes(ymin =  cilow,
                    ymax = cihi),
                width = 0.2, position = position_dodge(0.9)) +
  labs(y = "Change in absolute error",
       fill = "Network Type") +
  scale_fill_manual(values = c("Centralized" = "#de2d26", "Egalitarian" = "#3182bd")) +
  guides(pattern = guide_legend(override.aes = list(fill = "white")),
         fill = guide_legend(override.aes = list(pattern = "none")))+
  theme(axis.text=element_text(size=10), 
        axis.title=element_text(size=12), 
        axis.title.x=element_blank(),
        panel.background = element_rect(fill='white'),
        plot.background = element_rect(fill='white', color=NA),
        legend.position="right", 
        plot.title = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        axis.text.x=element_text(size=10), 
        panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  geom_text(aes(x = 0.8, y = -600, label='Decrease\nin Error'), size = 5) +
  geom_text(aes(x = 4.05, y = 1400, label='Increase\nin Error'), size = 5) +
  coord_cartesian(xlim = c(0.5,4.5))
ggsave("figures/Figure-3-Networks-Composition-Pattern.png",width = 6,height = 3.5,units = "in",dpi = 300)


################################## Figure 4 ###################################
df_network <- df %>% filter(structure != 'control') %>% drop_na(c('mag_revision'))
df_1 <- df_network %>% 
  select(instance_id, node_id, network, composition, structure, err_1_abs,mag_revision) %>% drop_na()
df_1 <- df_1 %>% mutate(quantile = ntile(err_1_abs, 10))
q <- quantile(df_1$err_1_abs, probs = seq(0, 1, 1/10),na.rm=TRUE)

labs = c("")
for (i in 1:10){
  labs <- c(labs,paste(min(df_1[(df_1$quantile == i ), ]$err_1_abs),'-\n',max(df_1[(df_1$quantile == i ), ]$err_1_abs)))
}

fig4 <- df_1 %>%
  group_by(quantile, structure) %>%
  summarise(mean.ch = mean(mag_revision, na.rm = TRUE),
            sd.error = sd(mag_revision, na.rm = TRUE),
            n.error = n()) %>%
  mutate(se.error = sd.error / sqrt(n.error),
         lower.ci.error = mean.ch - qt(1 - (0.05 / 2), n.error - 1) * se.error,
         upper.ci.error = mean.ch + qt(1 - (0.05 / 2), n.error - 1) * se.error)
fig4$Structure <- factor(fig4$structure, levels=c('c','d','control'))
levels(fig4$Structure) <- c('Centralized','Egalitarian','Control')
fig4 <- fig4 %>% filter(Structure == "Egalitarian")

ggplot(fig4,aes(x=quantile,y=mean.ch,ymin=lower.ci.error, ymax =upper.ci.error)) +
  geom_point(size=2,color='#1F77B4') + geom_line(size=0.5,color='#1F77B4') + 
  ylab('Magnitude of revision') + xlab('Initial Error (Quantiles)') + 
  theme(axis.text=element_text(size=15), 
        axis.title=element_text(size=15), 
        axis.title.x=element_text(size=9),
        axis.title.y=element_text(size=9),
        panel.background = element_rect(fill='white'),
        plot.background = element_rect(fill='white', color=NA),
        plot.title = element_text(hjust = 0),
        #axis.ticks.x = element_blank(),
        axis.text.y=element_text(size=8),
        axis.text.x=element_text(size=8),
        legend.position = "top",
        legend.box = "horizontal",
        legend.key = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) + 
  geom_errorbar(width = 0.2, size = 0.5, color='#1F77B4') +
  scale_x_continuous(breaks = c(1:10), labels = labs[2:11])
ggsave("figures/Figure-4-Revision-Coefficient.png",width = 5,height = 3.5,units = "in",dpi = 300)

JonckheereTerpstraTest(df_1$mag_revision,df_1$quantile)



################################## Figure 5 ###################################
df_group <- df %>% group_by(network) %>% 
  drop_na(c('err_1_abs','err_3_abs','mag_revision')) %>%
  summarise(mean_err_1 = mean(err_1_abs), 
            mean_err_3 = mean(err_3_abs),
            change_error = mean_err_3 - mean_err_1,
            r_c = cor(err_1_abs,mag_revision,method='spearman'),
            homo = first(composition),
            structure = first(structure),
            mean_1 = mean(value_1),
            sd_1 = sd(value_1),
            mean_3 = mean(value_3),
            change_mean = mean_3 - mean_1)

df_group$Structure <- factor(df_group$structure, levels=c('c','d'), labels=c('Centralized','Egalitarian'))


a <- ggplot(df_group, aes(x=r_c, color = Structure, fill = Structure)) +
  geom_histogram(binwidth = .1, alpha=0.5, position="identity") +
  scale_fill_manual(values=c("#de2d26", "#3182bd")) + 
  scale_color_manual(values=c("#de2d26", "#3182bd")) + 
  xlab('Revision Coefficient') + ylab('Frequency') + 
  theme_bw() + # Add a dashed line at y = 0
  geom_vline(xintercept = 0, linewidth=1, linetype="dashed", color = "black") +
  geom_hline(yintercept = 0, linewidth=1, color = "black") +
  coord_cartesian(xlim=c(-0.86, 0.86), ylim=c(0, 5.9)) +
  theme(legend.position = "none")

b <- ggplot(df_group, aes(x=change_error, color = Structure, fill = Structure)) +
  geom_histogram(binwidth = 300, alpha=0.5, position="identity") +
  scale_fill_manual(values=c("#de2d26", "#3182bd")) + 
  scale_color_manual(values=c("#de2d26", "#3182bd")) + 
  xlab('Change in Error') + ylab('') + 
  theme_bw() + # Add a dashed line at y = 0
  geom_vline(xintercept = 0, linewidth=1, linetype="dashed", color = "black") +
  geom_hline(yintercept = 0, linewidth=1, color = "black") +
  annotate("text", x=1200, y=5.5, label="Increase Error", color="black", size=3.5) + # Add "Less Accurate" annotation
  annotate("text", x=-1200, y=5.5, label="Reduce Error", color="black", size=3.5) + 
  coord_cartesian(xlim=c(-2200, 2200), ylim=c(0, 5.9)) 
ggarrange(a, b, ncol = 2, nrow = 1, widths = c(1, 1.45))
ggsave("figures/Figure_5-revision-coef-hists.png",width = 7,height = 3.5,units = "in",dpi = 300)

df_cent <- df_group %>% filter(structure =='c')
df_egal <- df_group %>% filter(structure =='d')
cor.test(df_cent$r_c,df_cent$change_error,method='spearman')
cor.test(df_egal$r_c,df_egal$change_error,method='spearman')
mean(df_cent$change_error)
wilcox.test(df_cent$change_error)
mean(df_egal$change_error)
wilcox.test(df_egal$change_error)

### Included still? - CALVIN
df_t <- df %>% filter(structure == 'c')
means_1 <- df_t %>%
  group_by(instance_id) %>%
  summarize(mean_value_1 = mean(value_1, na.rm = TRUE),
            mean_value_3 = mean(value_3, na.rm = TRUE))

df_t <- df %>% filter(structure == 'd')
means_2 <- df_t %>%
  group_by(instance_id) %>%
  summarize(mean_value_1 = mean(value_1, na.rm = TRUE),
            mean_value_3 = mean(value_3, na.rm = TRUE))

means_1$Structure <- 'Centralized'
means_2$Structure <- 'Egalitarian'
means <- rbind(means_1, means_2)

means_long <- means %>%
  pivot_longer(cols = c(mean_value_1, mean_value_3), 
               names_to = "round", 
               values_to = "mean_value") %>%
  mutate(round = ifelse(round == "mean_value_1", 1, 3))

ggplot(means) +
  geom_segment(aes(x = 1, y = mean_value_1, xend = 3, yend = mean_value_3, color = Structure), 
               size = 0.5, alpha = 0.15) + # Thin, slightly transparent lines
  geom_point(aes(x = 1, y = mean_value_1, color = Structure, alpha = 0.1)) +
  geom_point(aes(x = 3, y = mean_value_3, color = Structure, alpha = 0.1)) +
  geom_smooth(data = means_long, aes(x = round, y = mean_value, color = Structure), 
              method = "lm", se = FALSE, size = 2) + # Thicker line for lm function
  theme_bw() + xlab('Round') + ylab('Mean Estimate') + 
  scale_color_manual(values=c("#de2d26", "#3182bd")) +
  theme(legend.position = "none") + facet_wrap(~Structure) + 
  scale_x_continuous(breaks = c(1, 3), limits = c(.7, 3.3)) + 
  coord_cartesian(ylim = c(100, 4900))
ggsave("figures/Figure_5b-Network-changes.png",width = 8,height = 4,units = "in",dpi = 300)


################################## Chat Analysis ###################################
df_anger <- read_csv('data/Network Chat Anger (n 32 nets) 23.09.28.csv')
df_anger <- df_anger %>% drop_na()

# Get network types
df_ac <- df %>% select(composition, structure, network) %>% filter(structure == 'c' & composition == 'aa')
df_ad <- df %>% select(composition, structure, network) %>% filter(structure == 'd' & composition == 'aa')
df_mc <- df %>% select(composition, structure, network) %>% filter(structure == 'c' & composition == 'm')
df_md <- df %>% select(composition, structure, network) %>% filter(structure == 'd' & composition == 'm')
labels <- c('homophilous centralized','homophilous egalitarian','diverse centralized',
            'diverse egalitarian')
new <- c()
for (i in df_anger$instance_id) {
  if (i %in% unique(df_ac$network)) {
    new =c(new,labels[1])
  } else if (i %in% unique(df_ad$network)){
    new =c(new,labels[2])
  } else if (i %in% unique(df_mc$network)){
    new =c(new,labels[3])
  }  else if (i %in% unique(df_md$network)) {
    new = c(new,labels[4])
  }
}
df_anger$network <- factor(new, levels = labels)

df_network_improvements <- df_table_networks %>% select(network, mean_change)
df_network_improvements$instance_id <- df_network_improvements$network
df_anger <- merge(df_anger,df_network_improvements, by='instance_id')
df_anger$network <- df_anger$network.x

df_ch <- df_anger %>% filter(network == 'homophilous centralized') 
df_eh <- df_anger %>% filter(network == 'homophilous egalitarian') 
df_cm <- df_anger %>% filter(network == 'diverse centralized') 
df_em <- df_anger %>% filter(network == 'diverse egalitarian') 

# No difference in diverse
wilcox.test(df_cm$liwc.anger,df_em$liwc.anger)

# Difference in homophilous
print(paste(mean(df_ch$liwc.anger),mean(df_eh$liwc.anger)))
wilcox.test(df_ch$liwc.anger,df_eh$liwc.anger)

# Correlation
cor.test(df_anger$liwc.anger,df_anger$mean_change,method='spearman')





################################## SUPPLEMENTAL ###################################
######################### Fig S1. Map of participants ##########
state_freq <- as.data.frame(table(df$state))
colnames(state_freq) <- c("state", "freq")

convert_state_abbreviation_to_name <- function(abbreviations) {
  state_indices <- match(abbreviations, state.abb)
  full_names <- state.name[state_indices]
  lower_case_names <- tolower(full_names)
  lower_case_names <- ifelse(abbreviations == "DC", "district of columbia", lower_case_names)
  
  return(lower_case_names)
}
state_freq$ID <- convert_state_abbreviation_to_name(state_freq$state)

plot_usmap(regions = "states", data = state_freq, values = 'freq') + 
  scale_fill_gradient(name = "Number of\nParticipants",
                      low = "white", high = "#1F456E", trans='log10', 
                      limits = c(1, 700), na.value = "lightgray") +
  theme( axis.text = element_blank(),
         axis.title = element_blank(),
         axis.ticks = element_blank(),
         panel.background = element_rect(fill = "white", colour = "white"),
         plot.background = element_rect(fill = "white", colour = "white"),
         legend.position = "right",
         legend.text = element_text(size = 16),
         legend.title = element_text(size = 20))
ggsave("figures/Figure-S1-Map.png",width = 10,height = 7,units = "in",dpi = 300)


######################### Fig S8. Centralized networks impact ##########
df_group <- df %>% group_by(network) %>% 
  drop_na(c('err_1_abs','err_3_abs','mag_revision')) %>%
  summarise(mean_err_1 = mean(err_1_abs), 
            mean_err_3 = mean(err_3_abs),
            change_error = mean_err_3 - mean_err_1,
            r_c = cor(err_1_abs,mag_revision,method='spearman'),
            homo = first(composition),
            structure = first(structure),
            mean_1 = mean(value_1),
            sd_1 = sd(value_1),
            mean_3 = mean(value_3),
            change_mean = mean_3 - mean_1)

# To find p-value that the average initial guesses are above the truth, we simply 
# perform a wilcoxon test on mean_1 minus the truth (35)
wilcox.test(df_group$mean_1 - 35)

# Now see if significant for each individual condition
wilcox.test((df_group %>% filter(structure == 'c', homo == 'aa') %>% pull(mean_1)) - 35)
wilcox.test((df_group %>% filter(structure == 'd', homo == 'aa') %>% pull(mean_1)) - 35)
wilcox.test((df_group %>% filter(structure == 'c', homo == 'm') %>% pull(mean_1)) - 35)
wilcox.test((df_group %>% filter(structure == 'd', homo == 'm') %>% pull(mean_1)) - 35)

# Now let's make a nice plot for the centralized ones
df_cent <- df %>% filter(structure == 'c')

# Plot histogram of value 1
means_1 <- df_cent %>%
  group_by(instance_id) %>%
  summarize(mean_value = mean(value_1, na.rm = TRUE))

# Plot with vertical lines at the mean for each network
ggplot(df_cent, aes(x = value_1)) +
  geom_histogram(binwidth = 1000, alpha = 0.5) + 
  facet_wrap(~instance_id) +
  geom_vline(data = means_1, aes(xintercept = mean_value), color = "blue", size = 1.5) +
  geom_segment(aes(x = 35, xend = 35, y = 0, yend = 17.5),  color = "black", size = 1.5) +
  annotate("segment", x = 8000, xend = 8000, y = 17.5, yend = 5, color = "red", size = 1.5, 
           arrow = arrow(type = "closed", length = unit(0.1, "inches")))+ 
  annotate("text", x = 5700, y = 19, label = "Central Node", color = "red", hjust = 0, size = 3) +
  annotate("text", x = -600, y = 19, label = "Truth", color = "black", hjust = 0, size = 3) +
  geom_text(data = means_1, aes(x = mean_value, y = Inf, label = "Mean"), color = "blue", hjust = -0.2, vjust = 10, size = 3) +
  theme_bw() + xlab('Estimate 1') + ylab("") +
  theme(strip.text = element_blank(),  # Turns off facet titles
        plot.title = element_blank())

ggsave("figures/Figure-S8-central-figure-SI.png",width = 10,height = 8,units = "in",dpi = 300)



################################# Sensitivity Analysis #############################
print(paste('N =',length(df$X)))
df_n <- df  %>% drop_na(c('value_1','value_2','value_3'))
print(paste('N =',length(df_n$X)))
df_table_networks <- df_n %>% 
  select(network, structure, change_abs_err, composition, value_1, value_3,err_1_abs) %>%
  group_by(structure,network, composition) %>%
  dplyr::summarise(
    mean_r1 = mean(value_1, na.rm = T),
    mean_r3 = mean(value_3, na.rm = T),
    mean_change = mean(change_abs_err, na.rm = T),
    n_people = length(value_1)
  )
#Crowd estimates
df_table_networks$Err_1_Net <- abs(df_table_networks$mean_r1 - 35)
df_table_networks$Err_3_Net <- abs(df_table_networks$mean_r3 - 35)
df_table_networks$Mean_Change_Net <- df_table_networks$Err_3_Net - df_table_networks$Err_1_Net

df_table_networks$structure <- factor(df_table_networks$structure, levels=c('c','d'))
levels(df_table_networks$structure) <- c('Centralized','Egalitarian')


df_c <- df_table_networks %>% filter(structure == 'Centralized')
df_e <- df_table_networks %>% filter(structure == 'Egalitarian')

# Individual Error
print(paste('Centralized: average change in error was ',mean(df_c$mean_change), 'which is a',
            100*mean(df_c$mean_change)/mean(df_table_networks$mean_r1), '% increase.'))
wilcox.test(df_c$mean_change)
print(paste('Egalitarian: average change in error was ',mean(df_e$mean_change), 'which is a',
            100*mean(df_e$mean_change)/mean(df_table_networks$mean_r1), '% increase.'))
wilcox.test(df_e$mean_change)
wilcox.test(df_e$mean_change, df_c$mean_change)

#Differences by Homophily 
df_h <- df_table_networks %>% filter(composition == 'aa')
df_m <- df_table_networks %>% filter(composition == 'm')


# Individual Error 
print(paste('Homo: average change in error was ',mean(df_h$mean_change), 'which is a',
            100*mean(df_h$mean_change)/mean(df_table_networks$mean_r1), '% increase.'))
wilcox.test(df_h$mean_change)
print(paste('Egalitarian: average change in error was ',mean(df_m$mean_change), 'which is a',
            100*mean(df_m$mean_change)/mean(df_table_networks$mean_r1), '% increase.'))
wilcox.test(df_m$mean_change)


# Differences Both homophily and structure
df_ch <- df_table_networks %>% filter(composition == 'aa' & structure == 'Centralized')
df_cm <- df_table_networks %>% filter(composition == 'm' & structure == 'Centralized')
df_eh <- df_table_networks %>% filter(composition == 'aa' & structure == 'Egalitarian')
df_em <- df_table_networks %>% filter(composition == 'm' & structure == 'Egalitarian')

# Individual Error
print(paste('Egalitarian-Diverse: average change in error was ',mean(df_em$mean_change), 'which is a',
            100*mean(df_em$mean_change)/mean(df_table_networks$mean_r1), '% increase.'))
wilcox.test(df_em$mean_change)
print(paste('Centralized-Diverse: average change in error was ',mean(df_cm$mean_change), 'which is a',
            100*mean(df_cm$mean_change)/mean(df_table_networks$mean_r1), '% increase.'))
wilcox.test(df_cm$mean_change)
wilcox.test(df_cm$mean_change, df_em$mean_change)
print(paste('Egalitarian-Homo: average change in error was ',mean(df_eh$mean_change), 'which is a',
            100*mean(df_eh$mean_change)/mean(df_table_networks$mean_r1), '% increase.'))
wilcox.test(df_eh$mean_change)
print(paste('Centralized-Homo: average change in error was ',mean(df_ch$mean_change), 'which is a',
            100*mean(df_ch$mean_change)/mean(df_table_networks$mean_r1), '% increase.'))
wilcox.test(df_ch$mean_change)
wilcox.test(df_ch$mean_change, df_eh$mean_change)


########################## Results with clustered standard errors ###########
df = read.csv('data/all_data_2023.09.27.csv')
df_confeds <- df %>% filter(confederate == 'yes')
df$err_1_abs <- abs(df$err_1)
df$err_3_abs <- abs(df$err_3)
df$change_abs_err <- df$err_3_abs- df$err_1_abs
df <- df %>% filter(confederate == 'no')

df$structure <- factor(df$structure, levels=c('c','d'))
levels(df$structure) <- c('Centralized','Egalitarian')

lm_fit <- lm(change_abs_err ~ structure, data = df)
clustered_se <- coeftest(lm_fit, vcov = vcovCL(lm_fit, cluster = df$instance_id))

# Table S3
print(clustered_se)



# Homophilly
df$composition <- factor(df$composition, levels=c('aa','m'))
levels(df$composition) <- c('Homo','Diverse')

lm_fit <- lm(change_abs_err ~ composition, data = df)
clustered_se <- coeftest(lm_fit, vcov = vcovCL(lm_fit, cluster = df$instance_id))

# Table S4
print(clustered_se)


############################### Demographics, Table S1 ###########################
df = read.csv('data/all_data_2023.09.27.csv')
df_confeds <- df %>% filter(confederate == 'yes')
df$err_1_abs <- abs(df$err_1)
df$err_3_abs <- abs(df$err_3)
df$change_abs_err <- df$err_3_abs- df$err_1_abs
df <- df %>% filter(confederate == 'no')

df$ethnicity_recode <- case_when(
  df$ethnicity == "black" ~ "black",
  df$ethnicity %in% c("white", "white/caucasian") ~ "white/caucasian",
  TRUE ~ "other/mult-racial"
)

df <- df %>% drop_na('change_abs_err')

# Group and summarize
summary_table_1 <- df %>%
  group_by(network, composition, structure, ethnicity_recode) %>%
  tally() %>%
  spread(key = ethnicity_recode, value = n, fill = 0)


summary_table_2 <- df %>%
  group_by(network, composition, structure, gender) %>%
  tally() %>%
  spread(key = gender, value = n, fill = 0)

df$vaccine_recode <- case_when(
  df$health_vaccine_self == "yes" ~ 1,
  TRUE ~ 0
)

summary_table_3 <- df %>%
  group_by(network, composition, structure) %>%
  summarize(
    avg_age = mean(age, na.rm = TRUE),
    percent_vaccine = round(mean(vaccine_recode),3) * 100,
    count = n()
  )

summary_table <- merge(summary_table_1,summary_table_2, by=c('network', 'composition', 'structure'))
summary_table <- merge(summary_table,summary_table_3, by=c('network', 'composition', 'structure'))

summary_table$structure <- factor(summary_table$structure, levels=c('c','d'))
levels(summary_table$structure) <- c('Centralized','Egalitarian')
summary_table$composition <- factor(summary_table$composition, levels=c('aa','m'))
levels(summary_table$composition) <- c('Homophilous','Diverse')

write.csv(summary_table,"results/TableS1-Demographics.csv")

############################## Table S2 ########################################
df = read.csv('data/all_data_2023.09.27.csv')
df_confeds <- df %>% filter(confederate == 'yes')

df$err_1_abs <- abs(df$err_1)
df$err_2_abs <- abs(df$err_2)
df$err_3_abs <- abs(df$err_3)
df$change_abs_err <- df$err_3_abs- df$err_1_abs
df$mag_revision <- abs(df$value_3 - df$value_1)

# Remove confederates
df <- df %>% filter(confederate == 'no')

df_table_networks <- df %>% 
  select(network, structure, change_abs_err, composition, err_1_abs, err_3_abs, mag_revision) %>%
  group_by(network, structure, composition) %>%
  dplyr::summarise(
    err_1 = mean(err_1_abs, na.rm = T),
    err_3 = mean(err_3_abs, na.rm = T),
    mean_change = mean(change_abs_err, na.rm = T),
    mean_change_perc = mean(change_abs_err, na.rm = T) / mean(err_1_abs, na.rm = T),
    rc = cor.test(err_1_abs,mag_revision)$estimate
  )

df_table_networks$structure <- factor(df_table_networks$structure, levels=c('c','d','control'))
levels(df_table_networks$structure) <- c('Centralized','Egalitarian','Control')
df_table_networks <- df_table_networks %>% filter(structure != 'Control')

df_table_networks$composition <- factor(df_table_networks$composition, levels=c('aa','m'))
levels(df_table_networks$composition) <- c('Homophilous','Diverse')
write.csv(df_table_networks,"results/TableS2.csv")



############################## Table S5 ########################################
df_anger_2 <- read.csv("data/no_confeds_chats_2023.09.27.csv")

df_anger_2 <- df_anger_2 %>% select(instance_id, structure, condition, content)
df_anger_2$structure <- factor(df_anger_2$structure, levels=c('c','d'))
levels(df_anger_2$structure) <- c('Centralized','Egalitarian')
df_anger_2$condition <- factor(df_anger_2$condition, levels=c('aa.c','aa.d','m.c','m.d'))
levels(df_anger_2$condition) <- c('Homophilous','Homophilous','Diverse','Diverse')
df_chats <- df_anger_2 %>%
  group_by(instance_id, structure, condition) %>%
  summarise(concatenated_content = paste(content, collapse = " / ")) %>%
  ungroup()
write.csv(df_chats,"results/TableS5-Chat_Messages.csv")


################################# Figure S9 Chat ###############################

all_data <- "data/chat_code_final.xlsx"
sheets <- readxl::excel_sheets(all_data)
df <- readxl::read_excel(all_data, sheet = sheets[1])
df <- df %>% select(instance_id, Code_h)
# Go through the rest of the sheets and add them 1 by 1.
for (x in 2:length(sheets)) {
  df_2 <- readxl::read_excel(all_data, sheet = sheets[x])
  df_2 <- df_2 %>% select(instance_id, Code_h)
  df <- rbind(df,df_2)
}

df_separated <- df %>%
  separate_rows(Code_h, sep = ";\\s*")

df_separated <- df_separated %>%
  mutate(Code_h = case_when(
    Code_h == "vaccination intention" ~ "vaccine other",
    Code_h == "vaccine attitude" ~ "vaccine other",
    Code_h == "vaccine information" ~ "vaccine other",
    Code_h == "vaccine safety" ~ "vaccine other",
    Code_h == "vaccine efficacy" ~ "vaccine other",
    Code_h == "vaccine harm" ~ "vaccine other",
    Code_h == "vaccine mandate" ~ "vaccine other",
    Code_h == "vaccine confidence" ~ "vaccine other",
    Code_h == "vaccine trust" ~ "vaccine other",
    TRUE ~ Code_h  # Default case to keep other labels unchanged
  ))

code_counts <- df_separated %>%
  count(Code_h) %>%
  arrange(desc(n)) 
code_counts$Category <- factor(code_counts$Code_h, levels=code_counts$Code_h)

ggplot(code_counts, aes(x = Category, y = n, fill = Category)) +
  geom_bar(stat = "identity") + # Tell ggplot to use the counts ('n') as bar heights
  theme_minimal() + # A minimalistic theme
  labs(
    y = "Count"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x axis labels for better readability
  ) +
  scale_fill_grey(start = 0.2, end = 0.8)
ggsave("figures/Figure-S7_Chat-Annotations.png",width = 10,height = 7,units = "in",dpi = 300)


################################# Figure S9 Neighborhoods ###############################
df = read.csv('data/all_data_2023.09.27.csv')
df_t <- read.csv('data/Node-Location.csv')

df <- left_join(df, df_t, by = 'node_id') # Make sure confeds are in

# Only need egalitarian networks
df <- df %>% filter(structure =='d') %>% 
  filter(instance_id != 411) # remove one network of 40

# Get the "neighborhood" guess for each node based on neighbors
egal <- data.frame(
  Node = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30),
  N1 = c(8, 3, 2, 11, 2, 8, 20, 1, 16, 17, 2, 17, 3, 11, 1, 5, 10, 1, 12, 6, 9, 1, 7, 14, 4, 7, 6, 9, 2, 16),
  N2 = c(15, 5, 13, 15, 16, 20, 23, 6, 21, 21, 4, 19, 14, 13, 3, 9, 12, 5, 20, 7, 10, 4, 8, 22, 6, 10, 7, 10, 3, 18),
  N3 = c(18, 11, 15, 22, 18, 25, 26, 11, 27, 26, 8, 23, 17, 24, 4, 22, 13, 21, 21, 13, 18, 16, 12, 26, 12, 24, 9, 17, 5, 19),
  N4 = c(22, 29, 29, 25, 29, 27, 27, 23, 28, 28, 14, 25, 20, 25, 23, 30, 28, 30, 30, 19, 19, 24, 15, 30, 14, 27, 26, 29, 28, 24)
)

get_neighborhood <- function(n) {
  neighbors <- c(egal$N1[n], egal$N2[n], egal$N3[n], egal$N4[n])
  return(neighbors)
}

# For each network, get the neighborhood guesses
df$Neighbor_Mean_1 <- NA
df$Neighbor_Mean_2 <- NA
df$Neighbor_Mean_3 <- NA
for (i in 1:nrow(df)){
  df_t <- df %>% filter(instance_id == df$instance_id[i])
  neighbors <- c(get_neighborhood(df$number_in_pattern[i]), df$number_in_pattern[i])
  df$Neighbor_Mean_1[i] <- mean(df_t %>% filter(number_in_pattern %in% neighbors) %>% pull(err_1), na.rm = T)
  df$Neighbor_Mean_2[i] <- mean(df_t %>% filter(number_in_pattern %in% neighbors) %>% pull(err_2), na.rm = T)
  df$Neighbor_Mean_3[i] <- mean(df_t %>% filter(number_in_pattern %in% neighbors) %>% pull(err_3), na.rm = T)
}

df$Change_Neigh_1 <- df$Neighbor_Mean_2 - df$Neighbor_Mean_1
df$Change_Neigh_2 <- df$Neighbor_Mean_3 - df$Neighbor_Mean_1


# Because confederates were randomly assigned, we need to make a new ordering to show
# where people are compared to confed
df$confed_position <- NA
for (i in 1:nrow(df)){
  df_t <- df %>% filter(instance_id == df$instance_id[i])
  confed_val <- df_t %>% filter(confederate == 'yes') %>% pull(number_in_pattern)
  neighbors <- get_neighborhood(confed_val)
  if (df$confederate[i] == 'yes') {
    df$confed_position[i] <- 2
  } else if (df$number_in_pattern[i] %in% neighbors){
    df$confed_position[i] <- 1
  } else {
    df$confed_position[i] <- 0
  }
}

df <- df %>%
  mutate(confed_position = factor(confed_position, 
        labels = c("Other", "Confederate's neighbors", "Confederate")))

# Group by confed_position and calculate means
df_plot <- df %>% group_by(confed_position) %>%
  summarise(mean_change_1 = mean(Change_Neigh_1, na.rm = TRUE), 
            mean_change_2 = mean(Change_Neigh_2, na.rm = TRUE),
            se_change_1 = sd(Change_Neigh_1, na.rm = TRUE) / sqrt(n()),
            se_change_2 = sd(Change_Neigh_2, na.rm = TRUE) / sqrt(n()))

a <- ggplot(df_plot, aes(x = confed_position, y = mean_change_1, fill = confed_position)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, color = "black") + 
  geom_errorbar(aes(ymin = mean_change_1 - se_change_1, ymax = mean_change_1 + se_change_1), 
                position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = c("darkgreen", "orange", "red")) +
  annotate("text", label = "Increased error", x = 1.25, y = 450, color = "black") +
  annotate("text", label = "Reduced error", x = 3, y = -500, color = "black") +
  labs(x = "Neighborhood", y = "Change in Neighborhood Error", fill = "Neighborhood",
       title = "Round 1 to 2 change") +  theme_bw() + theme(legend.position = "none") +
  coord_cartesian(ylim = c(-600, 500))


b <- ggplot(df_plot, aes(x = confed_position, y = mean_change_2, fill = confed_position)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, color = "black") + 
  geom_errorbar(aes(ymin = mean_change_2 - se_change_2, ymax = mean_change_2 + se_change_2), 
                position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_manual(values = c("darkgreen", "orange", "red")) +
  
  annotate("text", label = "Increased error", x = 1.25, y = 450, color = "black") +
  annotate("text", label = "Reduced error", x = 3, y = -500, color = "black") +
  labs(x = "Neighborhood", y = "", fill = "Neighborhood",
       title = "Round 1 to 3 change") +
  theme_bw() + coord_cartesian(ylim = c(-600, 500))

ggarrange(a, b, ncol = 2, widths = c(1, 1.5))
ggsave("figures/Figure-S9_neighborhood_change_means.png", width = 10, height = 5)

# stats: Error neighborhood is lower than error individual
mean(df %>% filter(confed_position == "Other") %>% pull(Neighbor_Mean_1), na.rm=T)
mean(df %>% filter(confed_position == "Other") %>% pull(err_1), na.rm=T)
wilcox.test(df %>% filter(confed_position == "Other") %>% pull(Neighbor_Mean_1), 
            df %>% filter(confed_position == "Other") %>% pull(err_1), paired = T)

# stats: Error neighborhood improves over time for those not attached to confederate
mean(df %>% pull(Neighbor_Mean_1), na.rm=T) - mean(df %>% pull(Neighbor_Mean_3), na.rm=T)
wilcox.test(df %>% pull(Neighbor_Mean_1), df %>% pull(Neighbor_Mean_3), paired = T)

cor.test(df$Neighbor_Mean_1, df$err_change_total, na.rm=T, method='spearman')

mean(df %>% filter(confed_position == "Other") %>% pull(Change_Neigh_2), na.rm=T)
wilcox.test(df %>% filter(confed_position == "Other") %>% pull(Change_Neigh_2))

mean(df %>% filter(confed_position == "Confederate's neighbors") %>% pull(Change_Neigh_2), na.rm=T)
wilcox.test(df %>% filter(confed_position == "Confederate's neighbors") %>% pull(Change_Neigh_2))

mean(df %>% filter(confed_position == "Confederate") %>% pull(Change_Neigh_2), na.rm=T)
wilcox.test(df %>% filter(confed_position == "Confederate") %>% pull(Change_Neigh_2))


################################# Figure S10 Simulations ###############################
# Generate a normal distribution based on mean and sd of our data
df = read.csv('data/all_data_2023.09.27.csv')
mean_value <- mean(df %>% filter(confederate == "no") %>% pull(value_1), na.rm = TRUE)
sd_value <- sd(df %>% filter(confederate == "no") %>% pull(value_1), na.rm = TRUE)

get_value <- function(){
  value <- rnorm(1, mean = mean_value, sd = sd_value)
  if(value < 0){
    value <- 0
  }
  return(value)
}

# Simulate degroot model
df <- data.frame()
sim_degroot <- function(n) {
  # 30 initial guesses, confed is 8000
  guess_1 <- sapply(1:30, function(x) get_value())
  guess_1[1] <- 8000
  
  # Generate a correlation between self-weights and errors ~.6
  errors_1 <- abs(guess_1 - 35)
  self_weights = -errors_1*.5 + rnorm(30, 0, 1000) + rnorm(30, 0, 1000) 
  
  # Standardize self_weights between 0 and 1
  self_weights <- self_weights[-1]
  self_weights <- (self_weights - min(self_weights)) / (max(self_weights) - min(self_weights))
  self_weights <- c(1, self_weights) # Confederate has full self-weight.
  
  neighborhood_means <- c()
  for (i in 1:30){
    node <- egal$Node[i]
    neighbors <- c(egal$N1[i], egal$N2[i], egal$N3[i], egal$N4[i])
    neighborhood_means <- c(neighborhood_means, mean(guess_1[neighbors]))
  }
  guess_2 <- self_weights * guess_1 + (1 - self_weights) * neighborhood_means
  
  neighborhood_means_2 <- c()
  for (i in 1:30){
    node <- egal$Node[i]
    neighbors <- c(egal$N1[i], egal$N2[i], egal$N3[i], egal$N4[i])
    neighborhood_means_2 <- c(neighborhood_means_2, mean(guess_2[neighbors]))
  }
  guess_3 <- self_weights * guess_2 + (1 - self_weights) * neighborhood_means_2
  df_t <- data.frame(X = 1:30, Change_1 = (guess_2 - 35) - (guess_1 - 35), Change_2 = (guess_3 - 35) - (guess_1 - 35))
  df_t$Instance <- n
  return(df_t)
}
for (i in 1:1000){
  d_t <- sim_degroot(i)
  df <- rbind(df, d_t)
}

df <- df %>% group_by(X) %>% summarise(Change_1 =mean(Change_1), Change_2 = mean(Change_2))


### Make nice layout
egal_g <- make_empty_graph(n = nrow(egal), directed = TRUE)

for(i in 1:nrow(egal)){
  node <- egal$Node[i]
  neighbors <- c(egal$N1[i], egal$N2[i], egal$N3[i], egal$N4[i])
  for(neigh in neighbors){
    egal_g <- add_edges(egal_g, c(node, neigh))
  }
}
egal_g <- as.undirected(egal_g, mode="collapse")
N <- vcount(egal_g)
half_N <- N/2

# Define radii for inner and outer circles
r_inner <- 0.8
r_outer <- 1
angles <- seq(0, 2*pi, length.out = half_N + 1)
angles <- angles[-(half_N + 1)]  
x_inner <- r_inner * cos(angles)
y_inner <- r_inner * sin(angles)
x_outer <- r_outer * cos(angles)
y_outer <- r_outer * sin(angles)
x_coords <- c(x_inner, x_outer)
y_coords <- c(y_inner, y_outer)

label_size = 4
vertex_size = 8

df_network = df
# Get nice color gradient
color_gradient <- colorRampPalette(c("darkgreen",'green',"lightgreen","yellow", "orange","darkorange","red"))
i <- 500
values <- seq(-i, i, length.out = i*2)
colors <- color_gradient(length(values))
get_color <- function(value) {
  if (value < -i){
    value <- -i
  } else if (value > i-1){
    value <- i-1
  }
  scaled_value <- value + i+1
  colors[scaled_value]
}


png(filename="figures/Figure-S10_degroot-sims.png", width=2400, height=600)
par(mfrow=c(1,3), mar=c(5, 5, 5, 5), mgp=c(3, 1, 0))

# Round 1
#guesses <- df_network %>% select(value_1)
#guesses <- c(df_network %>% pull(value_1), rep(NA, 30 - length(df_network$X)))
example_colors <- rep("gray", 29)
example_colors <- c("red", example_colors)
plot(egal_g,  layout=matrix(c(x_coords, y_coords), ncol=2), vertex.label=NA, vertex.size=vertex_size, vertex.color=example_colors)
title(main="Round 1", cex.main=label_size, font.main=1.5)

# Round 2
#guesses <- df_network %>% select(value_2)
guesses <- c(df_network %>% pull(Change_1), rep(NA, 30 - length(df_network$X)))
example_colors <- sapply(guesses, get_color)
example_colors[1] <- "red"
plot(egal_g,  layout=matrix(c(x_coords, y_coords), ncol=2), vertex.label=NA, vertex.size=vertex_size, vertex.color=example_colors)
title(main="Round 2", cex.main=label_size, font.main=1.5)

# Round 3
#guesses <- df_network %>% select(value_3)
guesses <- c(df_network %>% pull(Change_2), rep(NA, 30 - length(df_network$X)))
example_colors <- sapply(guesses, get_color)
example_colors[1] <- "red"
plot(egal_g,  layout=matrix(c(x_coords, y_coords), ncol=2), vertex.label=NA, vertex.size=vertex_size, vertex.color=example_colors)
title(main="Round 3", cex.main=label_size, font.main=1.5)

dev.off()




##### Centralized ####
df <- data.frame()
sim_degroot_cent <- function(n) {
  # 30 initial guesses
  guess_1 <- sapply(1:30, function(x) get_value())
  guess_1[1] <- 8000
  
  errors_1 <- abs(guess_1 - 35)
  self_weights = -errors_1*.5 + rnorm(30, 0, 1000) + rnorm(30, 0, 1000) 
  
  # Standardize self_weights between 0 and 1
  self_weights <- self_weights[-1]
  self_weights <- (self_weights - min(self_weights)) / (max(self_weights) - min(self_weights))
  self_weights <- c(1, self_weights) # Confederate has full self-weight.
  
  # No neighborhoods in centralized
  guess_2 <- self_weights * guess_1 + (1 - self_weights) * 8000
  guess_3 <- self_weights * guess_2 + (1 - self_weights) * 8000
  df_t <- data.frame(X = 1:30, Change_1 = (guess_2 - 35) - (guess_1 - 35), Change_2 = (guess_3 - 35) - (guess_1 - 35))
  df_t$Instance <- n
  return(df_t)
}
for (i in 1:1000){
  d_t <- sim_degroot_cent(i)
  df <- rbind(df, d_t)
}

df <- df %>% group_by(X) %>% summarise(Change_1 =mean(Change_1), Change_2 = mean(Change_2))

star_g <- make_star(30, center=1, mode="undirected")


n <- vcount(star_g)
layout <- matrix(NA, nrow = n, ncol = 2)
layout[1, ] <- c(0, 0)
angle <- seq(0, 2 * pi, length.out = n)
layout[-1, ] <- cbind(cos(angle[-1]), sin(angle[-1]))

df_network = df

png(filename="degroot_cent-2.png", width=2400, height=600)
par(mfrow=c(1,3), mar=c(5, 5, 5, 5), mgp=c(3, 1, 0))



# Round 1
example_colors <- rep("gray", 29)
example_colors <- c("red", example_colors)
plot(star_g,  layout=layout, vertex.label=NA, vertex.size=vertex_size, vertex.color=example_colors)
title(main="Round 1", cex.main=label_size, font.main=1.5)

# Round 2
#guesses <- df_network %>% select(value_2)
guesses <- c(df_network %>% pull(Change_1), rep(NA, 30 - length(df_network$X)))[-1]
example_colors <- sapply(guesses, get_color)
example_colors <- c("red", example_colors)
plot(star_g, layout=layout,  vertex.label=NA, vertex.size=vertex_size, vertex.color=example_colors)
title(main="Round 2", cex.main=label_size, font.main=1.5)

# Round 3
#guesses <- df_network %>% select(value_3)
guesses <- c(df_network %>% pull(Change_2), rep(NA, 30 - length(df_network$X)))[-1]
example_colors <- sapply(guesses, get_color)
example_colors <- c("red", example_colors)
plot(star_g, layout=layout, vertex.label=NA, vertex.size=vertex_size, vertex.color=example_colors)
title(main="Round 3", cex.main=label_size, font.main=1.5)
dev.off()


