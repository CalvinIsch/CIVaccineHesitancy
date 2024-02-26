library(ggplot2); library(dplyr); library(ggpubr); library(igraph)
library(tidyverse); library(DescTools); library(lmtest); library(sandwich)
library(ggpattern); library(sf); library(viridis); library(maps); 
library(gridpattern)

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


############################## Figure 1 ##################################

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
       fill = "Network Type") +
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
  geom_text(aes(x = 0.4, y = -500, label='Decrease\nin Error'), size = 5) +
  geom_text(aes(x = 0.4, y = 1180, label='Increase\nin Error'), size = 5)  +
  coord_cartesian(xlim=c(0.5,2.4))
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


# Figure 3 ############################################
df_table_networks$composition <- factor(df_table_networks$composition, levels=c('aa','m'))
levels(df_table_networks$composition) <- c('Homophilous','Diverse')

df_table_networks$Conditions <- paste(df_table_networks$composition, df_table_networks$structure)
df_table_networks$Conditions <- factor(df_table_networks$Conditions,
                                       levels=c("Homophilous Centralized","Diverse Centralized","Homophilous Egalitarian",
                                                "Diverse Egalitarian"))



df_summary_1 = df_table_networks %>%
  group_by(Conditions) %>%
  dplyr::summarise(
    grand_mean_ind = mean(mean_change),
    se_ind = sd(mean_change)/sqrt(length(mean_change)),
    cilow=t.test(mean_change)$conf.int[1],
    cihi=t.test(mean_change)$conf.int[2]
  )
levels(df_summary_1$Conditions) <- c("Homo.\n Cent.","Diverse\n Cent.","Homo.\n Egal.",
                                          "Diverse\n Egal.")

ggplot(df_summary_1, aes(x = Conditions, y = grand_mean_ind, fill = Conditions)) +
  geom_bar(stat = "identity", position = "dodge", color='black') +
  geom_hline(yintercept = 0, color = 'black') + 
  geom_errorbar(aes(ymin =  cilow, #grand_mean - qnorm(0.975) * se,
                    ymax = cihi), #grand_mean + qnorm(0.975) * se),
                width = 0.2, position = position_dodge(0.9)) +
  labs(y = "Change in absolute error",
       fill = "Category") +
  scale_fill_manual(values=c("#de2d26", "#fc9272", "#3182bd", "#9ecae1", "#E4E4E4", "#808080")) + 
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
  geom_text(aes(x = 0.6, y = -670, label='Decrease\nin Error'), size = 5) +
  geom_text(aes(x = 4.15, y = 1400, label='Increase\nin Error'), size = 5) +
  coord_cartesian(xlim = c(0.4,4.6))
ggsave("figures/Figure-3-Networks-Composition.png",width = 6,height = 3.5,units = "in",dpi = 300)

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


################### Figure 4 ##########################
df_network <- df %>% filter(structure != 'control') %>% drop_na(c('mag_revision'))
df_1 <- df_network %>% 
  select(instance_id, node_id, network, composition, structure, err_1_abs,mag_revision) %>% drop_na()
df_1 <- df_1 %>% mutate(quantile = ntile(err_1_abs, 10))
q <- quantile(df_1$err_1_abs, probs = seq(0, 1, 1/10),na.rm=TRUE)

labs = c("")
for (i in 1:10){
  labs <- c(labs,paste(min(df_1[(df_1$quantile == i ), ]$err_1_abs),'-\n',max(df_1[(df_1$quantile == i ), ]$err_1_abs)))
}

labs = c("")
for (i in 1:10){
  labs <- c(labs,i)
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
  ylab('Magnitude of revision') + xlab('Initial Error (Deciles)') + 
  theme(axis.text=element_text(size=15), 
        axis.title=element_text(size=15), 
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
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



# Chat analysis ################################
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

df_a <- df_anger %>% select(network, liwc.anger)
df_c <- df_a %>% filter(grepl('centralized', network))
df_d <- df_a %>% filter(grepl('egalitarian', network))
print(paste(mean(df_c$liwc.anger),mean(df_d$liwc.anger)))
print(wilcox.test(df_c$liwc.anger,df_d$liwc.anger, paired = F))

df_network_improvements <- df_table_networks %>% select(network, mean_change)
df_network_improvements$instance_id <- df_network_improvements$network
df_anger <- merge(df_anger,df_network_improvements, by='instance_id')
df_anger$network <- df_anger$network.x

df_ch <- df_anger %>% filter(network == 'homophilous centralized') 
df_eh <- df_anger %>% filter(network == 'homophilous egalitarian') 
df_cm <- df_anger %>% filter(network == 'diverse centralized') 
df_em <- df_anger %>% filter(network == 'diverse egalitarian') 

wilcox.test(df_cm$liwc.anger,df_em$liwc.anger)
print(paste(mean(df_ch$liwc.anger),mean(df_eh$liwc.anger)))
wilcox.test(df_ch$liwc.anger,df_eh$liwc.anger)

cor.test(df_anger$liwc.anger,df_anger$mean_change,method='spearman')




# Individual level - SI
df_anger_2 <- read.csv("data/no_confeds_chats_2023.09.27.csv")
d_a <- df_anger_2 %>% select(node_id, liwc.anger)
d_b <- df %>% select(node_id,err_change_total)
d_a_2 <- merge(d_a, d_b) %>% drop_na()
aggregated_data <- d_a_2 %>%
  group_by(node_id) %>%
  summarise(
    avg_liwc.anger = mean(liwc.anger, na.rm = TRUE),
    avg_err_change_total = mean(err_change_total, na.rm = TRUE)
  )
ggplot(aggregated_data, aes(x = avg_liwc.anger, y = avg_err_change_total)) +
  geom_point(size = 1, alpha = 0.7) + # Color-coded points based on y-values, with some transparency
  geom_smooth(method = "lm", se = FALSE, color = "blue") + # Linear regression line without the standard error
  labs(
    title = 'Individual Scatterplot',
    x = "LIWC Anger",
    y = "Change in Error"
  ) +
  theme_minimal() + # Minimal theme
  theme(
    legend.position = "none" # No legend for color
  )

new <- c()
for (i in df_anger_2$instance_id) {
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
df_anger_2$network <- factor(new, levels = labels)

df_c <- df_anger_2 %>% filter(grepl('centralized', network))
df_d <- df_anger_2 %>% filter(grepl('egalitarian', network))
print(paste(mean(df_c$liwc.anger),mean(df_d$liwc.anger)))
print(wilcox.test(df_c$liwc.anger,df_d$liwc.anger, paired = F))

df_h <- df_anger_2 %>% filter(grepl('homophilous', network))
df_d <- df_anger_2 %>% filter(grepl('diverse', network))
print(paste(mean(df_h$liwc.anger),mean(df_d$liwc.anger)))
print(wilcox.test(df_h$liwc.anger,df_d$liwc.anger, paired = F))



# Create nice excel file
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


# Nice table for SI
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



## Sensitivity Analysis in SI
print(paste('N =',length(df$X)))
df <- df  %>% drop_na(c('value_1','value_2','value_3'))
print(paste('N =',length(df$X)))
df_table_networks <- df %>% 
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
print(clustered_se)
# Table S3


# Homophilly
df$composition <- factor(df$composition, levels=c('aa','m'))
levels(df$composition) <- c('Homo','Diverse')

lm_fit <- lm(change_abs_err ~ composition, data = df)
clustered_se <- coeftest(lm_fit, vcov = vcovCL(lm_fit, cluster = df$instance_id))
print(clustered_se)
### Creates Table S4





######################### Map of participants ##########
library(usmap)
 
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




# Table Demographics and n_complete by network ######################3 
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







################ Chat annotations ######
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

