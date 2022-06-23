############################# niche structure ##############################


##### correlation structure

scale_env_values = scale_spp_env[,2:5]

for (i in ncol(scale_env_values)){
  cor.test(scale_env_values[,1], scale_env_values[,4], method= "spearman")
}

###### PCA 
correlation = cor(scale_spp_env[,2:5])
pca = prcomp(scale_spp_env[,2:5], center = T)
stdev = pca$sdev / sum(pca$sdev)
pc_number =c("1", "2", "3", "4")
stdev = data.frame(pc_number, stdev)

### boots values
boot_stdev_df = data.frame(matrix(0, nrow=1, ncol=4))
boot_rotation_list = vector('list', 99)

for (i in 1:99){
  boot_num = sample(x=1:nrow(scale_spp_env[,2:5]), size = nrow(scale_spp_env[,2:5]), replace = T)
  boot_env_values = scale_spp_env[boot_num,2:5]
  boot_pca = prcomp(boot_env_values, center = T)
  boot_rotation_list[[i]] = boot_pca$rotation
  boot_stdev = boot_pca$sdev/ sum(boot_pca$sdev)
  boot_stdev_df = rbind(boot_stdev_df, boot_stdev)
}
boot_stdev_df = boot_stdev_df[-1,]

# generating stdev intervals
boot_stdev_means = apply(boot_stdev_df, 2,FUN=mean)
boot_stdev_quantiles = apply(boot_stdev_df, 2,FUN=function(x){quantile(x, probs = c(0.05,0.95))} )
boot_stdev_quantiles= t(boot_stdev_quantiles)
stdev = data.frame(stdev, boot_stdev_quantiles)

###### randomization procedure 
rand_correlation_list = vector('list', 99)
rand_stdev_df = data.frame(matrix(0, nrow=1, ncol=4))
rand_rotation_list = vector('list', 99)

for (i in 1:99){
  rand_env_values = scale_spp_env[,2:5]
  for (j in ncol(rand_env_values) ){
    rand_num = sample(x=1:nrow(rand_env_values), replace=F)
    rand_env_values[,j] = rand_env_values[rand_num, j]
  }
  rand_correlation_list[[i]] = cor(rand_env_values)
  rand_pca = prcomp(rand_env_values, center = T)
  rand_rotation_list[[i]] = rand_pca$rotation
  rand_stdev = rand_pca$sdev/ sum(rand_pca$sdev)
  rand_stdev_df = rbind(rand_stdev_df, rand_stdev)
}
rand_stdev_df = rand_stdev_df[-1,]

### collecting stdev statistics
rand_stdev_means = apply(rand_stdev_df, 2,FUN=mean)
rand_stdev_quantiles = apply(rand_stdev_df, 2,FUN=function(x){quantile(x, probs = c(0.05,0.95))} )
rand_stdev_stats= t(rbind(rand_stdev_means, rand_stdev_quantiles))
rand_stdev_stats = data.frame(pc_number, rand_stdev_stats)

### plotting
tiff("PC_variance.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data=rand_stdev_stats, aes(x=pc_number, y=rand_stdev_means*100) ) +
  geom_point(data=rand_stdev_stats, aes(y=rand_stdev_means*100), color="gray", size = 2) +
  geom_errorbar(aes(ymin=X5.*100, ymax=X95.*100), col="gray", size=1, width=0.25)+
  geom_point(data=stdev, aes(y=stdev*100), color="red", size = 1.5) +
  geom_errorbar(data=stdev, aes(ymin=X5.*100, ymax=X95.*100), col="red", size=1, width=0)+
  labs(x="PC axis", y="% explained variance")+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"),axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=12))
dev.off()

### collecting rotation statistics
rand_rotation_df = rand_rotation_list[[1]][,2]
for (i in 2:length(rand_rotation_list) ){
  rand_rotation_df = rbind(rand_rotation_df , rand_rotation_list[[i]][,2])
}
rand_rotation_means = apply(rand_rotation_df, 2,FUN=mean)
rand_rotation_quantiles = apply(rand_rotation_df, 2,FUN=function(x){quantile(x, probs = c(0.05,0.95))} )
rand_rotation_stats= t(rbind(rand_rotation_means, rand_rotation_quantiles))
variable = c("1","2","3","4")
rand_rotation_stats = data.frame(variable, rand_rotation_stats)

# observed rotation
boot_rotation_df = boot_rotation_list[[1]][,2]
for (i in 2:length(boot_rotation_list) ){
  boot_rotation_df = rbind(boot_rotation_df , boot_rotation_list[[i]][,2])
}

boot_rotation_quantiles = apply(boot_rotation_df, 2,FUN=function(x){quantile(x, probs = c(0.05,0.95))} )
boot_rotation_stats= t(rbind(pca$rotation[,2], boot_rotation_quantiles))
variable = c("1","2","3","4")
boot_rotation_stats = data.frame(variable, boot_rotation_stats)
colnames(boot_rotation_stats)[2] = "rotation"


### ploting
tiff("PC2_loadings.tiff", units="in", width=3.5, height=3, res=600)
ggplot(data=rand_rotation_stats, aes(x=variable, y=rand_rotation_means) ) +
  geom_point(data=rand_rotation_stats, aes(y=rand_rotation_means), color="gray", size = 2) +
  geom_errorbar(aes(ymin= X5., ymax= X95.), col="gray", size=1, width=0.25)+
  geom_point(data=boot_rotation_stats, aes(y=rotation), color="red", size = 1.5) +
  geom_errorbar(data=boot_rotation_stats, aes(ymin= X5., ymax= X95.), col="red", size=1, width=0)+
  labs(x="variable", y="PC2 loadings")+
  scale_x_discrete(labels=c("1" = "temp", "2" = "precip", "3" = "solar", "4"="soil"))+
  theme(panel.background=element_rect(fill="white"), panel.grid=element_line(colour=NULL),panel.border=element_rect(fill=NA,colour="black"), axis.title=element_text(size=14,face="bold"),axis.text.x=element_text(size=12))
dev.off()

rand_correlation_list[[1]][]