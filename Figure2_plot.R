#####################
#Plots 
library(ggplot2)
library(reshape2)
library(ggpubr)

#INSTRUCTIONS:
#First run all code in Figure2a_code.R, Figure2b_code.R, and Figure2c_code.R
#This code assumes particular directories for each of the results, but can be easily changed. 


##### 
n <- 250
#Plot 1: Gaussian, restricting each learner to only predict on test points from same cluster
filename1 <- "gaussian_restrict_update"
means_restrict <- read.csv(paste0("~/Desktop/rf_theory/", filename1, "/cs_means1.csv"), skip = 1, header = F)
means_restrict[, 1] = means_restrict[, 1] - 1
colnames(means_restrict) <- c("ncoef_list", "Merged", "Ensemble")
df.melted <- melt(means_restrict, id = "ncoef_list")

sds_restrict <- read.csv(paste0("~/Desktop/rf_theory/", filename1, "/cs_sds1.csv"), skip = 1, header = F)
sds_restrict[, 1] = sds_restrict[, 1] - 1
colnames(sds_restrict) <- c("ncoef_list", "Merged", "Ensemble")

ncoef_list <- c(2, 4, 8, 16, 32)

se.ensemble <- data.frame(cbind(ncoef_list, means_restrict$Ensemble - 1.96*sds_restrict$Ensemble/(sqrt(n)),
                                means_restrict$Ensemble + 1.96*sds_restrict$Ensemble/(sqrt(n))))
se.merged <- data.frame(cbind(ncoef_list, means_restrict$Merged - 1.96*sds_restrict$Merged/(sqrt(n)),
                              means_restrict$Merged + 1.96*sds_restrict$Merged/(sqrt(n))))

q1 <- ggplot() + 
  geom_line(data = df.melted, aes(x = ncoef_list, y = value, color = variable)) +
  geom_ribbon(data=se.merged,aes(x=ncoef_list,ymin=V2,ymax=V3),fill="#2A9D8F",alpha=0.2, inherit.aes = FALSE) +
  geom_ribbon(data=se.ensemble,aes(x=ncoef_list,ymin=V2,ymax=V3),fill="red4",alpha=0.2, inherit.aes = FALSE) + 
  theme_classic() + xlab("Number of clusters") + 
  #ylab("") + 
  scale_y_continuous("", limits = c(2.8, 4.5), breaks  = scales::pretty_breaks(n = 4)) + 
  #scale_x_continuous(xlab, limits = xlim, breaks  = scales::pretty_breaks(n = 5)) +
  labs(title = "Gaussian clusters") +
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)),
        plot.title = element_text(size = rel(1.3)), legend.text=element_text(size=rel(1.3)), 
        legend.title=element_text(size=rel(1.5)), legend.position = "none") +
  scale_color_manual(values=c("#2A9D8F","red4"), labels = c("Merged", "Ensemble"), name="Method")
show(q1)

#######
#Plot 2: Uniform clusters

filename2 <- "uniform"
means_uniform <- read.csv(paste0("~/Desktop/rf_theory/", filename2, "/cs_means1.csv"), skip = 1, header = F)
means_uniform[, 1] = means_uniform[, 1] - 1
colnames(means_uniform) <- c("ncoef_list", "Merged", "Ensemble")
df.melted1 <- melt(means_uniform, id = "ncoef_list")

sds_uniform <- read.csv(paste0("~/Desktop/rf_theory/", filename2, "/cs_sds1.csv"), skip = 1, header = F)
sds_uniform[, 1] = sds_uniform[, 1] - 1
colnames(sds_uniform) <- c("ncoef_list", "Merged", "Ensemble")

ncoef_list <- c(2, 4, 8, 16, 32)

se.ensemble1 <- data.frame(cbind(ncoef_list, means_uniform$Ensemble - 1.96*sds_uniform$Ensemble/(sqrt(n)),
                                 means_uniform$Ensemble + 1.96*sds_uniform$Ensemble/(sqrt(n))))
se.merged1 <- data.frame(cbind(ncoef_list, means_uniform$Merged - 1.96*sds_uniform$Merged/(sqrt(n)),
                               means_uniform$Merged + 1.96*sds_uniform$Merged/(sqrt(n))))

q2 <- ggplot() + 
  geom_line(data = df.melted1, aes(x = ncoef_list, y = value, color = variable)) +
  geom_ribbon(data=se.merged1,aes(x=ncoef_list,ymin=V2,ymax=V3),fill="#2A9D8F",alpha=0.2, inherit.aes = FALSE) +
  geom_ribbon(data=se.ensemble1,aes(x=ncoef_list,ymin=V2,ymax=V3),fill="red4",alpha=0.2, inherit.aes = FALSE) + 
  theme_classic() + xlab("Number of clusters") + ylab("") + 
  labs(title = "Uniform clusters") +
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)),
        plot.title = element_text(size = rel(1.3)), legend.text=element_text(size=rel(1.3)), 
        legend.title=element_text(size=rel(1.5)), legend.position = "none") +
  scale_color_manual(values=c("#2A9D8F","red4"), labels = c("Merged", "Ensemble"), name="Method")
show(q2)


###################
#Plot 4: Restricted laplace

filename3 <- "laplace"
means_allpoints <- read.csv(paste0("~/Desktop/rf_theory/", filename3, "/cs_means1.csv"), header = T)
#colnames(means_allpoints) <- c("ncoef_list", "Merged", "Ensemble")
means_mat1 <- as.matrix(cbind(means_allpoints$Merged, means_allpoints$Ensemble))
data_means1 <- data.frame(cbind(ncoef_list, means_mat1))
df.melted1 <- melt(data_means1, id = "ncoef_list")

sds_allpoints <- read.csv(paste0("~/Desktop/rf_theory/", filename3, "/cs_sds1.csv"), header = T)
sds_allpoints[, 1] = sds_allpoints[, 1] - 1
#colnames(sds_allpoints) <- c("ncoef_list", "Merged", "Ensemble")

ncoef_list <- c(2, 4, 8, 16, 32)
n1 <- 250

se.ensemble1 <- data.frame(cbind(ncoef_list, means_allpoints$Ensemble - 1.96*sds_allpoints$Ensemble/(sqrt(n1)),
                                 means_allpoints$Ensemble + 1.96*sds_allpoints$Ensemble/(sqrt(n1))))
se.merged1 <- data.frame(cbind(ncoef_list, means_allpoints$Merged - 1.96*sds_allpoints$Merged/(sqrt(n1)),
                               means_allpoints$Merged + 1.96*sds_allpoints$Merged/(sqrt(n1))))

q4 <- ggplot() + 
  geom_line(data = df.melted1, aes(x = ncoef_list, y = value, color = variable)) +
  geom_ribbon(data=se.merged1,aes(x=ncoef_list,ymin=V2,ymax=V3),fill="#2A9D8F",alpha=0.2, inherit.aes = FALSE) +
  geom_ribbon(data=se.ensemble1,aes(x=ncoef_list,ymin=V2,ymax=V3),fill="red4",alpha=0.2, inherit.aes = FALSE) + 
  theme_classic() + xlab("Number of clusters") + ylab("") + 
  labs(title = "Laplace Clusters") +
  theme(axis.title=element_text(size=rel(1.3)), axis.text=element_text(size=rel(1.3)),
        plot.title = element_text(size = rel(1.3)), legend.text=element_text(size=rel(1.3)), 
        legend.title=element_text(size=rel(1.5)), legend.position = "none") +
  scale_color_manual(values=c("#2A9D8F","red4"), labels = c("Merged", "Ensemble"), name="Method")
show(q4)



#composite plot
rf_theory_composite <- ggarrange(q2, q1, q4, labels = c("A", "B", "C"),
                                 ncol = 3, nrow = 1, common.legend = TRUE, legend = "bottom")
show(rf_theory_composite)

annotate_figure(rf_theory_composite, left = text_grob("Average RMSE", rot = 90, size = 14),
                top = text_grob("Random Forest SCLs", size = 20, face = "bold"))

