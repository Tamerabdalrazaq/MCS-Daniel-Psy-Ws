rm(list = ls()) 

library(dplyr)
library(ggplot2)

#Visualizations
## visualize the PoL effect ----
#### distance ----

Path = "~/models_181026/"
saving_path = "~/Figures/"
#####microwaves ---- 
#n-1
model_2 = readRDS(paste0(Path,"n1_m_PoL_Effect_distance_2_18102026.RDS"))
model_3 = readRDS(paste0(Path,"n1_m_PoL_Effect_distance_3_18102026.RDS"))   
model_4 = readRDS(paste0(Path,"n1_m_PoL_Effect_distance_4_18102026.RDS"))
model_5 = readRDS(paste0(Path,"n1_m_PoL_Effect_distance_5_18102026.RDS"))   
model_6 = readRDS(paste0(Path,"n1_m_PoL_Effect_distance_6_18102026.RDS"))   

### visualize angle models
samp_2_dis = insight::get_parameters(model_2)%>%mutate(Saccade = 1)
samp_3_dis = insight::get_parameters(model_3)%>%mutate(Saccade = 2)
samp_4_dis = insight::get_parameters(model_4)%>%mutate(Saccade = 3)
samp_5_dis = insight::get_parameters(model_5)%>%mutate(Saccade = 4)
samp_6_dis = insight::get_parameters(model_6)%>%mutate(Saccade = 5)

samp_all <- bind_rows(samp_2_dis, samp_3_dis, samp_4_dis, samp_5_dis,samp_6_dis)

#install.packages("paletteer")
library(paletteer)
my_p = paletteer_d("awtools::a_palette",n =5, direction = -1)

library(ggdist)
t1 = ggplot(samp_all, aes(x = b_Intercept,y =as.factor(Saccade), fill =as.factor(Saccade))) +
  stat_halfeye() +
  xlim(-25, 25) +
  #coord_flip() +
  #facet_wrap(~ Saccade, nrow = 1) +
  theme(legend.position = "none") + xlab("PoL effect (distance)") + ylab("Saccade order") + 
  scale_y_discrete(limits = rev,labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) + 
  #  scale_x_discrete(labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) +
  #scale_fill_viridis_d(option = "plasma") + 
  scale_fill_paletteer_d("awtools::a_palette")+ 
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  theme_classic() +
  theme(legend.position = "none", text = element_text(size=25), 
        axis.title.x = element_text(face="bold"), 
        axis.title.y = element_text(face="bold"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour=my_p),
        plot.title = element_text(hjust = 0.5))  

ggsave(filename = paste0(saving_path,"m_n1.png"),
       plot = t1,width = 1920, height = 1080, units = "px")

##n-2 
model_2 = readRDS(paste0(Path,"n2_m_PoL_Effect_distance_2_19102026.RDS"))
model_3 = readRDS(paste0(Path,"n2_m_PoL_Effect_distance_3_19102026.RDS"))
model_4 = readRDS(paste0(Path,"n2_m_PoL_Effect_distance_4_19102026.RDS"))
model_5 = readRDS(paste0(Path,"n2_m_PoL_Effect_distance_5_19102026.RDS"))
model_6 = readRDS(paste0(Path,"n2_m_PoL_Effect_distance_6_19102026.RDS"))

### visualize angle models
samp_2_dis = insight::get_parameters(model_2)%>%mutate(Saccade = 1)
samp_3_dis = insight::get_parameters(model_3)%>%mutate(Saccade = 2)
samp_4_dis = insight::get_parameters(model_4)%>%mutate(Saccade = 3)
samp_5_dis = insight::get_parameters(model_5)%>%mutate(Saccade = 4)
samp_6_dis = insight::get_parameters(model_6)%>%mutate(Saccade = 5)

samp_all <- bind_rows(samp_2_dis, samp_3_dis, samp_4_dis, samp_5_dis,samp_6_dis)

#install.packages("paletteer")
library(paletteer)


library(ggdist)
t3 = ggplot(samp_all, aes(x = b_Intercept,y =as.factor(Saccade), fill =as.factor(Saccade))) +
  stat_halfeye() +
  xlim(-25, 25) +
  #coord_flip() +
  #facet_wrap(~ Saccade, nrow = 1) +
  theme(legend.position = "none") + xlab("PoL effect (distance)") + ylab("Saccade order") + 
  scale_y_discrete(limits = rev,labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) + 
  #  scale_x_discrete(labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) +
  #scale_fill_viridis_d(option = "plasma") + 
  scale_fill_paletteer_d("awtools::a_palette")+ 
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  theme_classic() +
  theme(legend.position = "none", text = element_text(size=25), 
        axis.title.x = element_text(face="bold"), 
        axis.title.y = element_text(face="bold"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour=my_p),
        plot.title = element_text(hjust = 0.5)) 

ggsave(filename = paste0(saving_path,"m_n2.png"),
       plot = t3,width = 1920, height = 1080, units = "px")

##n-3
#microwaves
model_2 = readRDS(paste0(Path,"n3_m_PoL_Effect_distance_2_19102026.RDS"))
model_3 = readRDS(paste0(Path,"n3_m_PoL_Effect_distance_3_19102026.RDS"))
model_4 = readRDS(paste0(Path,"n3_m_PoL_Effect_distance_4_19102026.RDS"))
model_5 = readRDS(paste0(Path,"n3_m_PoL_Effect_distance_5_19102026.RDS"))
model_6 = readRDS(paste0(Path,"n3_m_PoL_Effect_distance_6_19102026.RDS"))

### visualize angle models
samp_2_dis = insight::get_parameters(model_2)%>%mutate(Saccade = 1)
samp_3_dis = insight::get_parameters(model_3)%>%mutate(Saccade = 2)
samp_4_dis = insight::get_parameters(model_4)%>%mutate(Saccade = 3)
samp_5_dis = insight::get_parameters(model_5)%>%mutate(Saccade = 4)
samp_6_dis = insight::get_parameters(model_6)%>%mutate(Saccade = 5)

samp_all <- bind_rows(samp_2_dis, samp_3_dis, samp_4_dis, samp_5_dis,samp_6_dis)

#install.packages("paletteer")
library(paletteer)


library(ggdist)
t5 = ggplot(samp_all, aes(x = b_Intercept,y =as.factor(Saccade), fill =as.factor(Saccade))) +
  stat_halfeye() +
  xlim(-25, 25) +
  #coord_flip() +
  #facet_wrap(~ Saccade, nrow = 1) +
  theme(legend.position = "none") + xlab("PoL effect (distance)") + ylab("Saccade order") + 
  scale_y_discrete(limits = rev,labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) + 
  #  scale_x_discrete(labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) +
  #scale_fill_viridis_d(option = "plasma") + 
  scale_fill_paletteer_d("awtools::a_palette")+ 
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  theme_classic() +
  theme(legend.position = "none", text = element_text(size=25), 
        axis.title.x = element_text(face="bold"), 
        axis.title.y = element_text(face="bold"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour=my_p),
        plot.title = element_text(hjust = 0.5)) 

ggsave(filename = paste0(saving_path,"m_n3.png"),
       plot = t5,width = 1920, height = 1080, units = "px")


#####Clocks ---- 
## n-1
model_2 = readRDS(paste0(Path,"n1_c_PoL_Effect_distance_2_19102026.RDS"))
model_3 = readRDS(paste0(Path,"n1_c_PoL_Effect_distance_3_19102026.RDS"))
model_4 = readRDS(paste0(Path,"n1_c_PoL_Effect_distance_4_19102026.RDS"))
model_5 = readRDS(paste0(Path,"n1_c_PoL_Effect_distance_5_19102026.RDS"))
model_6 = readRDS(paste0(Path,"n1_c_PoL_Effect_distance_6_19102026.RDS"))

### visualize angle models
samp_2_dis = insight::get_parameters(model_2)%>%mutate(Saccade = 1)
samp_3_dis = insight::get_parameters(model_3)%>%mutate(Saccade = 2)
samp_4_dis = insight::get_parameters(model_4)%>%mutate(Saccade = 3)
samp_5_dis = insight::get_parameters(model_5)%>%mutate(Saccade = 4)
samp_6_dis = insight::get_parameters(model_6)%>%mutate(Saccade = 5)

samp_all <- bind_rows(samp_2_dis, samp_3_dis, samp_4_dis, samp_5_dis,samp_6_dis)

#install.packages("paletteer")
library(paletteer)


library(ggdist)
t2 = ggplot(samp_all, aes(x = b_Intercept,y =as.factor(Saccade), fill =as.factor(Saccade))) +
  stat_halfeye() +
  xlim(-25, 25) +
  #coord_flip() +
  #facet_wrap(~ Saccade, nrow = 1) +
  theme(legend.position = "none") + xlab("PoL effect (distance)") + ylab("Saccade order") + 
  scale_y_discrete(limits = rev,labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) + 
  #  scale_x_discrete(labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) +
  #scale_fill_viridis_d(option = "plasma") + 
  scale_fill_paletteer_d("awtools::a_palette")+ 
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  theme_classic() +
  theme(legend.position = "none", text = element_text(size=25), 
        axis.title.x = element_text(face="bold"), 
        axis.title.y = element_text(face="bold"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour=my_p),
        plot.title = element_text(hjust = 0.5)) 

ggsave(filename = paste0(saving_path,"c_n1.png"),
       plot = t2,width = 1920, height = 1080, units = "px")

#clocks n-2
model_2 = readRDS(paste0(Path,"n2_c_PoL_Effect_distance_2_19102026.RDS"))
model_3 = readRDS(paste0(Path,"n2_c_PoL_Effect_distance_3_19102026.RDS"))
model_4 = readRDS(paste0(Path,"n2_c_PoL_Effect_distance_4_19102026.RDS"))
model_5 = readRDS(paste0(Path,"n2_c_PoL_Effect_distance_5_19102026.RDS"))
model_6 = readRDS(paste0(Path,"n2_c_PoL_Effect_distance_6_19102026.RDS"))

### visualize angle models
samp_2_dis = insight::get_parameters(model_2)%>%mutate(Saccade = 1)
samp_3_dis = insight::get_parameters(model_3)%>%mutate(Saccade = 2)
samp_4_dis = insight::get_parameters(model_4)%>%mutate(Saccade = 3)
samp_5_dis = insight::get_parameters(model_5)%>%mutate(Saccade = 4)
samp_6_dis = insight::get_parameters(model_6)%>%mutate(Saccade = 5)

samp_all <- bind_rows(samp_2_dis, samp_3_dis, samp_4_dis, samp_5_dis,samp_6_dis)

#install.packages("paletteer")
library(paletteer)


library(ggdist)
t4 = ggplot(samp_all, aes(x = b_Intercept,y =as.factor(Saccade), fill =as.factor(Saccade))) +
  stat_halfeye() +
  xlim(-25, 25) +
  #coord_flip() +
  #facet_wrap(~ Saccade, nrow = 1) +
  theme(legend.position = "none") + xlab("PoL effect (distance)") + ylab("Saccade order") + 
  scale_y_discrete(limits = rev,labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) + 
  #  scale_x_discrete(labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) +
  #scale_fill_viridis_d(option = "plasma") + 
  scale_fill_paletteer_d("awtools::a_palette")+ 
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  theme_classic() +
  theme(legend.position = "none", text = element_text(size=25), 
        axis.title.x = element_text(face="bold"), 
        axis.title.y = element_text(face="bold"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour=my_p),
        plot.title = element_text(hjust = 0.5)) 

ggsave(filename =  paste0(saving_path,"/c_n2.png"),
       plot = t4,width = 1920, height = 1080, units = "px")


##n-3
model_2 = readRDS(paste0(Path,"n3_c_PoL_Effect_distance_2_19102026.RDS"))
model_3 = readRDS(paste0(Path,"n3_c_PoL_Effect_distance_3_19102026.RDS"))
model_4 = readRDS(paste0(Path,"n3_c_PoL_Effect_distance_4_19102026.RDS"))
model_5 = readRDS(paste0(Path,"n3_c_PoL_Effect_distance_5_19102026.RDS"))
model_6 = readRDS(paste0(Path,"n3_c_PoL_Effect_distance_6_19102026.RDS"))

### visualize angle models
samp_2_dis = insight::get_parameters(model_2)%>%mutate(Saccade = 1)
samp_3_dis = insight::get_parameters(model_3)%>%mutate(Saccade = 2)
samp_4_dis = insight::get_parameters(model_4)%>%mutate(Saccade = 3)
samp_5_dis = insight::get_parameters(model_5)%>%mutate(Saccade = 4)
samp_6_dis = insight::get_parameters(model_6)%>%mutate(Saccade = 5)

samp_all <- bind_rows(samp_2_dis, samp_3_dis, samp_4_dis, samp_5_dis,samp_6_dis)

#install.packages("paletteer")
library(paletteer)

library(ggdist)
t6 = ggplot(samp_all, aes(x = b_Intercept,y =as.factor(Saccade), fill =as.factor(Saccade))) +
  stat_halfeye() +
  xlim(-25, 26) +
  #coord_flip() +
  #facet_wrap(~ Saccade, nrow = 1) +
  theme(legend.position = "none") + xlab("PoL effect (distance)") + ylab("Saccade order") + 
  scale_y_discrete(limits = rev,labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) + 
  #  scale_x_discrete(labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) +
  #scale_fill_viridis_d(option = "plasma") + 
  scale_fill_paletteer_d("awtools::a_palette")+ 
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  theme_classic() +
  theme(legend.position = "none", text = element_text(size=25), 
        axis.title.x = element_text(face="bold"), 
        axis.title.y = element_text(face="bold"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour=my_p),
        plot.title = element_text(hjust = 0.5)) 

ggsave(filename = paste0(saving_path,"c_n3.png"),
       plot = t6,width = 1920, height = 1080, units = "px")

####saccdic direction ---- 
###### Microwave ----
#n-1
model_2 = readRDS(paste0(Path,"n1_m_PoL_Effect_angles_2_18102026.RDS"))
model_3 = readRDS(paste0(Path,"n1_m_PoL_Effect_angles_3_18102026.RDS"))
model_4 = readRDS(paste0(Path,"n1_m_PoL_Effect_angles_4_18102026.RDS"))
model_5 = readRDS(paste0(Path,"n1_m_PoL_Effect_angles_5_18102026.RDS"))
model_6 = readRDS(paste0(Path,"n1_m_PoL_Effect_angles_6_18102026.RDS"))

### visualize angle models
samp_2_dis = insight::get_parameters(model_2)%>%mutate(Saccade = 1)
samp_3_dis = insight::get_parameters(model_3)%>%mutate(Saccade = 2)
samp_4_dis = insight::get_parameters(model_4)%>%mutate(Saccade = 3)
samp_5_dis = insight::get_parameters(model_5)%>%mutate(Saccade = 4)
samp_6_dis = insight::get_parameters(model_6)%>%mutate(Saccade = 5)

samp_all <- bind_rows(samp_2_dis, samp_3_dis, samp_4_dis, samp_5_dis,samp_6_dis)

#install.packages("paletteer")
library(paletteer)


library(ggdist)
q1 = ggplot(samp_all, aes(x = b_Intercept,y =as.factor(Saccade), fill =as.factor(Saccade))) +
  stat_halfeye() +
  xlim(-0.2, 0.2) +
  #coord_flip() +
  #facet_wrap(~ Saccade, nrow = 1) +
  theme(legend.position = "none") + xlab("PoL effect (cosine-similarity)") + ylab("Saccade order") + 
  scale_y_discrete(limits = rev,labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) + 
  #  scale_x_discrete(labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) +
  #scale_fill_viridis_d(option = "plasma") + 
  scale_fill_paletteer_d("awtools::a_palette")+ 
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  theme_classic() +
  theme(legend.position = "none", text = element_text(size=25), 
        axis.title.x = element_text(face="bold"), 
        axis.title.y = element_text(face="bold"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour=my_p),
        plot.title = element_text(hjust = 0.5)) 

ggsave(filename = paste0(saving_path,"/n1_cos_m.png"),
       plot = q1,width = 1920, height = 1080, units = "px")


##n-2 
#microwaves
model_2 = readRDS(paste0(Path,"n2_m_PoL_Effect_angles_2_19102026.RDS"))
model_3 = readRDS(paste0(Path,"n2_m_PoL_Effect_angles_3_19102026.RDS"))
model_4 = readRDS(paste0(Path,"n2_m_PoL_Effect_angles_4_19102026.RDS"))
model_5 = readRDS(paste0(Path,"n2_m_PoL_Effect_angles_5_19102026.RDS"))
model_6 = readRDS(paste0(Path,"n2_m_PoL_Effect_angles_6_19102026.RDS"))

### visualize angle models
samp_2_dis = insight::get_parameters(model_2)%>%mutate(Saccade = 1)
samp_3_dis = insight::get_parameters(model_3)%>%mutate(Saccade = 2)
samp_4_dis = insight::get_parameters(model_4)%>%mutate(Saccade = 3)
samp_5_dis = insight::get_parameters(model_5)%>%mutate(Saccade = 4)
samp_6_dis = insight::get_parameters(model_6)%>%mutate(Saccade = 5)

samp_all <- bind_rows(samp_2_dis, samp_3_dis, samp_4_dis, samp_5_dis,samp_6_dis)

#install.packages("paletteer")
library(paletteer)


library(ggdist)
q3 = ggplot(samp_all, aes(x = b_Intercept,y =as.factor(Saccade), fill =as.factor(Saccade))) +
  stat_halfeye() +
  xlim(-0.2, 0.2) +
  #coord_flip() +
  #facet_wrap(~ Saccade, nrow = 1) +
  theme(legend.position = "none") + xlab("PoL effect (cosine-similarity)") + ylab("Saccade order") + 
  scale_y_discrete(limits = rev,labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) + 
  #  scale_x_discrete(labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) +
  #scale_fill_viridis_d(option = "plasma") + 
  scale_fill_paletteer_d("awtools::a_palette")+ 
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  theme_classic() +
  theme(legend.position = "none", text = element_text(size=25), 
        axis.title.x = element_text(face="bold"), 
        axis.title.y = element_text(face="bold"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour=my_p),
        plot.title = element_text(hjust = 0.5)) 

ggsave(filename = paste0(saving_path,"n2_cos_m.png"),
       plot = q3,width = 1920, height = 1080, units = "px")

#n-3

##n-3
#microwaves
model_2 = readRDS(paste0(Path,"n3_m_PoL_Effect_angles_2_19102026.RDS"))
model_3 = readRDS(paste0(Path,"n3_m_PoL_Effect_angles_3_19102026.RDS"))
model_4 = readRDS(paste0(Path,"n3_m_PoL_Effect_angles_4_19102026.RDS"))
model_5 = readRDS(paste0(Path,"n3_m_PoL_Effect_angles_5_19102026.RDS"))
model_6 = readRDS(paste0(Path,"n3_m_PoL_Effect_angles_6_19102026.RDS"))

### visualize angle models
samp_2_dis = insight::get_parameters(model_2)%>%mutate(Saccade = 1)
samp_3_dis = insight::get_parameters(model_3)%>%mutate(Saccade = 2)
samp_4_dis = insight::get_parameters(model_4)%>%mutate(Saccade = 3)
samp_5_dis = insight::get_parameters(model_5)%>%mutate(Saccade = 4)
samp_6_dis = insight::get_parameters(model_6)%>%mutate(Saccade = 5)

samp_all <- bind_rows(samp_2_dis, samp_3_dis, samp_4_dis, samp_5_dis,samp_6_dis)

#install.packages("paletteer")
library(paletteer)


library(ggdist)
q5 = ggplot(samp_all, aes(x = b_Intercept,y =as.factor(Saccade), fill =as.factor(Saccade))) +
  stat_halfeye() +
  xlim(-0.2, 0.2) +
  #coord_flip() +
  #facet_wrap(~ Saccade, nrow = 1) +
  theme(legend.position = "none") + xlab("PoL effect (cosine-similarity)") + ylab("Saccade order") + 
  scale_y_discrete(limits = rev,labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) + 
  #  scale_x_discrete(labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) +
  #scale_fill_viridis_d(option = "plasma") + 
  scale_fill_paletteer_d("awtools::a_palette")+ 
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  theme_classic() +
  theme(legend.position = "none", text = element_text(size=25), 
        axis.title.x = element_text(face="bold"), 
        axis.title.y = element_text(face="bold"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour=my_p),
        plot.title = element_text(hjust = 0.5)) 

ggsave(filename = paste0(saving_path,"n3_cos_m.png"),
       plot = q5,width = 1920, height = 1080, units = "px")




#####clocks ----
#n-1
model_2 = readRDS(paste0(Path,"n1_c_PoL_Effect_angles_2_19102026.RDS"))
model_3 = readRDS(paste0(Path,"n1_c_PoL_Effect_angles_3_19102026.RDS"))
model_4 = readRDS(paste0(Path,"n1_c_PoL_Effect_angles_4_19102026.RDS"))
model_5 = readRDS(paste0(Path,"n1_c_PoL_Effect_angles_5_19102026.RDS"))
model_6 = readRDS(paste0(Path,"n1_c_PoL_Effect_angles_6_19102026.RDS"))

### visualize angle models
samp_2_dis = insight::get_parameters(model_2)%>%mutate(Saccade = 1)
samp_3_dis = insight::get_parameters(model_3)%>%mutate(Saccade = 2)
samp_4_dis = insight::get_parameters(model_4)%>%mutate(Saccade = 3)
samp_5_dis = insight::get_parameters(model_5)%>%mutate(Saccade = 4)
samp_6_dis = insight::get_parameters(model_6)%>%mutate(Saccade = 5)

samp_all <- bind_rows(samp_2_dis, samp_3_dis, samp_4_dis, samp_5_dis,samp_6_dis)

#install.packages("paletteer")
library(paletteer)


library(ggdist)
q2 = ggplot(samp_all, aes(x = b_Intercept,y =as.factor(Saccade), fill =as.factor(Saccade))) +
  stat_halfeye() +
  xlim(-0.2, 0.2) +
  #coord_flip() +
  #facet_wrap(~ Saccade, nrow = 1) +
  theme(legend.position = "none") + xlab("PoL effect (cosine-similarity)") + ylab("Saccade order") + 
  scale_y_discrete(limits = rev,labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) + 
  #  scale_x_discrete(labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) +
  #scale_fill_viridis_d(option = "plasma") + 
  scale_fill_paletteer_d("awtools::a_palette")+ 
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  theme_classic() +
  theme(legend.position = "none", text = element_text(size=25), 
        axis.title.x = element_text(face="bold"), 
        axis.title.y = element_text(face="bold"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour=my_p),
        plot.title = element_text(hjust = 0.5)) 


ggsave(filename = paste0(saving_path,"n1_cos_c.png"),
       plot = q2,width = 1920, height = 1080, units = "px")


##n-2
model_2 = readRDS(paste0(Path,"n2_c_PoL_Effect_angles_2_19102026.RDS"))
model_3 = readRDS(paste0(Path,"n2_c_PoL_Effect_angles_3_19102026.RDS"))
model_4 = readRDS(paste0(Path,"n2_c_PoL_Effect_angles_4_19102026.RDS"))
model_5 = readRDS(paste0(Path,"n2_c_PoL_Effect_angles_5_19102026.RDS"))
model_6 = readRDS(paste0(Path,"n2_c_PoL_Effect_angles_6_19102026.RDS"))

### visualize angle models
samp_2_dis = insight::get_parameters(model_2)%>%mutate(Saccade = 1)
samp_3_dis = insight::get_parameters(model_3)%>%mutate(Saccade = 2)
samp_4_dis = insight::get_parameters(model_4)%>%mutate(Saccade = 3)
samp_5_dis = insight::get_parameters(model_5)%>%mutate(Saccade = 4)
samp_6_dis = insight::get_parameters(model_6)%>%mutate(Saccade = 5)

samp_all <- bind_rows(samp_2_dis, samp_3_dis, samp_4_dis, samp_5_dis,samp_6_dis)

#install.packages("paletteer")
library(paletteer)


library(ggdist)
q4 = ggplot(samp_all, aes(x = b_Intercept,y =as.factor(Saccade), fill =as.factor(Saccade))) +
  stat_halfeye() +
  xlim(-0.2, 0.2) +
  #coord_flip() +
  #facet_wrap(~ Saccade, nrow = 1) +
  theme(legend.position = "none") + xlab("PoL effect (cosine-similarity)") + ylab("Saccade order") + 
  scale_y_discrete(limits = rev,labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) + 
  #  scale_x_discrete(labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) +
  #scale_fill_viridis_d(option = "plasma") + 
  scale_fill_paletteer_d("awtools::a_palette")+ 
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  theme_classic() +
  theme(legend.position = "none", text = element_text(size=25), 
        axis.title.x = element_text(face="bold"), 
        axis.title.y = element_text(face="bold"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour=my_p),
        plot.title = element_text(hjust = 0.5)) 

ggsave(filename = paste0(saving_path,"n2_cos_c.png"),
       plot = q4,width = 1920, height = 1080, units = "px")


#n-3
model_2 = readRDS(paste0(Path,"n3_c_PoL_Effect_angles_2_19102026.RDS"))
model_3 = readRDS(paste0(Path,"n3_c_PoL_Effect_angles_3_19102026.RDS"))
model_4 = readRDS(paste0(Path,"n3_c_PoL_Effect_angles_4_19102026.RDS"))
model_5 = readRDS(paste0(Path,"n3_c_PoL_Effect_angles_5_19102026.RDS"))
model_6 = readRDS(paste0(Path,"n3_c_PoL_Effect_angles_6_19102026.RDS"))

### visualize angle models
samp_2_dis = insight::get_parameters(model_2)%>%mutate(Saccade = 1)
samp_3_dis = insight::get_parameters(model_3)%>%mutate(Saccade = 2)
samp_4_dis = insight::get_parameters(model_4)%>%mutate(Saccade = 3)
samp_5_dis = insight::get_parameters(model_5)%>%mutate(Saccade = 4)
samp_6_dis = insight::get_parameters(model_6)%>%mutate(Saccade = 5)

samp_all <- bind_rows(samp_2_dis, samp_3_dis, samp_4_dis, samp_5_dis,samp_6_dis)

#install.packages("paletteer")
library(paletteer)


library(ggdist)
q6 = ggplot(samp_all, aes(x = b_Intercept,y =as.factor(Saccade), fill =as.factor(Saccade))) +
  stat_halfeye() +
  xlim(-0.2, 0.2) +
  #coord_flip() +
  #facet_wrap(~ Saccade, nrow = 1) +
  theme(legend.position = "none") + xlab("PoL effect (cosine-similarity)") + ylab("Saccade order") + 
  scale_y_discrete(limits = rev,labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) + 
  #  scale_x_discrete(labels = c("1" = "1st", "2" = "2nd", "3" = "3rd", "4" = "4th","5" = "5th")) +
  #scale_fill_viridis_d(option = "plasma") + 
  scale_fill_paletteer_d("awtools::a_palette")+ 
  geom_vline(aes(xintercept = 0), linetype = "dashed")+
  theme_classic() +
  theme(legend.position = "none", text = element_text(size=25), 
        axis.title.x = element_text(face="bold"), 
        axis.title.y = element_text(face="bold"),
        axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour=my_p),
        plot.title = element_text(hjust = 0.5)) 

ggsave(filename = paste0(saving_path,"n3_cos_c.png"),
       plot = q6,width = 1920, height = 1080, units = "px")


## combined - microwaves
library(patchwork)

p1 <- t1 + theme(
  # Increase x-axis TEXT (tick labels) size
  axis.text.x = element_text(colour = "black", size = 40),
  
  # Reduce axis TITLE size (x and y)
  axis.title.x = element_text(face = "bold", size = 35, margin = margin(t = 15)),
  axis.title.y = element_text(face = "bold", size = 35, margin = margin(r = 15)),
  
  # y-axis text (unchanged size, but keeping your color mapping)
  axis.text.y = element_text(colour = my_p, size = 40)
) + xlab("") 
p2 <- t3+ theme(
  # Increase x-axis TEXT (tick labels) size
  axis.text.x = element_text(colour = "black", size = 40),
  
  # Reduce axis TITLE size (x and y)
  axis.title.x = element_text(face = "bold", size = 35, margin = margin(t = 15)),
  axis.title.y = element_text(face = "bold", size = 35, margin = margin(r = 15)),
  
  # y-axis text (unchanged size, but keeping your color mapping)
  axis.text.y = element_text(colour = my_p, size = 40)
)  + theme(axis.title.y = element_blank(),        axis.text.y  = element_blank(),
           axis.ticks.y = element_blank(),axis.line.y = element_blank()) +
  xlab("PoL effect (Euclidean distance, px)") 

p3 <- t5 +  theme(
  # Increase x-axis TEXT (tick labels) size
  axis.text.x = element_text(colour = "black", size = 40),
  
  # Reduce axis TITLE size (x and y)
  axis.title.x = element_text(face = "bold", size = 35, margin = margin(t = 15)),
  axis.title.y = element_text(face = "bold", size = 35, margin = margin(r = 15)),
  
  # y-axis text (unchanged size, but keeping your color mapping)
  axis.text.y = element_text(colour = my_p, size = 40)
) + theme(axis.title.y = element_blank(),        axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),axis.line.y = element_blank()) + 
  xlab("") 




e1 <- q1 +  theme(
  # Increase x-axis TEXT (tick labels) size
  axis.text.x = element_text(colour = "black", size = 40),
  
  # Reduce axis TITLE size (x and y)
  axis.title.x = element_text(face = "bold", size = 35, margin = margin(t = 15)),
  axis.title.y = element_text(face = "bold", size = 35, margin = margin(r = 15)),
  
  # y-axis text (unchanged size, but keeping your color mapping)
  axis.text.y = element_text(colour = my_p, size = 40)
) + 
  xlab("") +  theme(plot.margin = margin(r = 30))
e2 <- q3  +  theme(
  # Increase x-axis TEXT (tick labels) size
  axis.text.x = element_text(colour = "black", size = 40),
  
  # Reduce axis TITLE size (x and y)
  axis.title.x = element_text(face = "bold", size = 35, margin = margin(t = 15)),
  axis.title.y = element_text(face = "bold", size = 35, margin = margin(r = 15)),
  
  # y-axis text (unchanged size, but keeping your color mapping)
  axis.text.y = element_text(colour = my_p, size = 40)
) + theme(axis.title.y = element_blank(),        axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),axis.line.y = element_blank()) + 
  xlab("PoL effect (Cosine similarity)") + theme(plot.margin = margin(r = 30))
e3 <- q5 +  theme(
  # Increase x-axis TEXT (tick labels) size
  axis.text.x = element_text(colour = "black", size = 40),
  
  # Reduce axis TITLE size (x and y)
  axis.title.x = element_text(face = "bold", size = 35, margin = margin(t = 15)),
  axis.title.y = element_text(face = "bold", size = 35, margin = margin(r = 15)),
  
  # y-axis text (unchanged size, but keeping your color mapping)
  axis.text.y = element_text(colour = my_p, size = 40)
) + theme(axis.title.y = element_blank(),        axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),axis.line.y = element_blank()) + 
  xlab("") 

combined <- p1 | p2 | p3

#combined <- combined +plot_annotation(tag_levels = "1", tag_prefix = "N-") &
#  theme(plot.tag = element_text(face = "bold", size = 25))

combined
ggsave(filename = paste0(saving_path,"m.svg"),
       plot = combined,width = 22,#22 
       height = 10, 
       units = "in", 
       dpi = 300)


combined_direction <- e1 | e2 | e3
ggsave(filename = paste0(saving_path,"m_dir.svg"),
       plot = combined_direction,width = 22,#22 
       height = 10, 
       units = "in", 
       dpi = 300)


###combined - clocks ----

cd1 <- t2 + theme(
  # Increase x-axis TEXT (tick labels) size
  axis.text.x = element_text(colour = "black", size = 40),
  
  # Reduce axis TITLE size (x and y)
  axis.title.x = element_text(face = "bold", size = 35, margin = margin(t = 15)),
  axis.title.y = element_text(face = "bold", size = 35, margin = margin(r = 15)),
  
  # y-axis text (unchanged size, but keeping your color mapping)
  axis.text.y = element_text(colour = my_p, size = 40)
) + xlab("") 
cd2 <- t4+ theme(
  # Increase x-axis TEXT (tick labels) size
  axis.text.x = element_text(colour = "black", size = 40),
  
  # Reduce axis TITLE size (x and y)
  axis.title.x = element_text(face = "bold", size = 35, margin = margin(t = 15)),
  axis.title.y = element_text(face = "bold", size = 35, margin = margin(r = 15)),
  
  # y-axis text (unchanged size, but keeping your color mapping)
  axis.text.y = element_text(colour = my_p, size = 40)
)  + theme(axis.title.y = element_blank(),        axis.text.y  = element_blank(),
           axis.ticks.y = element_blank(),axis.line.y = element_blank()) +
  xlab("PoL effect (Euclidean distance, px)") 

cd3 <- t6 +  theme(
  # Increase x-axis TEXT (tick labels) size
  axis.text.x = element_text(colour = "black", size = 40),
  
  # Reduce axis TITLE size (x and y)
  axis.title.x = element_text(face = "bold", size = 35, margin = margin(t = 15)),
  axis.title.y = element_text(face = "bold", size = 35, margin = margin(r = 15)),
  
  # y-axis text (unchanged size, but keeping your color mapping)
  axis.text.y = element_text(colour = my_p, size = 40)
) + theme(axis.title.y = element_blank(),        axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),axis.line.y = element_blank()) + 
  xlab("") 




ce1 <- q2 +  theme(
  # Increase x-axis TEXT (tick labels) size
  axis.text.x = element_text(colour = "black", size = 40),
  
  # Reduce axis TITLE size (x and y)
  axis.title.x = element_text(face = "bold", size = 35, margin = margin(t = 15)),
  axis.title.y = element_text(face = "bold", size = 35, margin = margin(r = 15)),
  
  # y-axis text (unchanged size, but keeping your color mapping)
  axis.text.y = element_text(colour = my_p, size = 40)
) + 
  xlab("") +  theme(plot.margin = margin(r = 30))
ce2 <- q4  +  theme(
  # Increase x-axis TEXT (tick labels) size
  axis.text.x = element_text(colour = "black", size = 40),
  
  # Reduce axis TITLE size (x and y)
  axis.title.x = element_text(face = "bold", size = 35, margin = margin(t = 15)),
  axis.title.y = element_text(face = "bold", size = 35, margin = margin(r = 15)),
  
  # y-axis text (unchanged size, but keeping your color mapping)
  axis.text.y = element_text(colour = my_p, size = 40)
) + theme(axis.title.y = element_blank(),        axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),axis.line.y = element_blank()) + 
  xlab("PoL effect (Cosine similarity)") + theme(plot.margin = margin(r = 30))
ce3 <- q6 +  theme(
  # Increase x-axis TEXT (tick labels) size
  axis.text.x = element_text(colour = "black", size = 40),
  
  # Reduce axis TITLE size (x and y)
  axis.title.x = element_text(face = "bold", size = 35, margin = margin(t = 15)),
  axis.title.y = element_text(face = "bold", size = 35, margin = margin(r = 15)),
  
  # y-axis text (unchanged size, but keeping your color mapping)
  axis.text.y = element_text(colour = my_p, size = 40)
) + theme(axis.title.y = element_blank(),        axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),axis.line.y = element_blank()) + 
  xlab("") 

combined <- cd1 | cd2 | cd3

#combined <- combined +plot_annotation(tag_levels = "1", tag_prefix = "N-") &
#  theme(plot.tag = element_text(face = "bold", size = 25))

combined
ggsave(filename = paste0(saving_path,"c_dis.svg"),
       plot = combined,width = 22,#22 
       height = 10, 
       units = "in", 
       dpi = 300)


combined_direction <- ce1 | ce2 | ce3
ggsave(filename = paste0(saving_path,"c_dir.svg"),
       plot = combined_direction,width = 22,#22 
       height = 10, 
       units = "in", 
       dpi = 300)

