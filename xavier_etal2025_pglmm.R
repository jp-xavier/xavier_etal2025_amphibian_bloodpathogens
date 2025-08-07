### Xavier et al., 2025 - Temperature and ecomorphology linked to blood pathogens incidence in neotropical amphibians
## Phylogenetic generalized linear mixed models (PGLMMs)
### By - Joao Paulo de Oliveira Xavier

# Packages ----

# install.packages(c("ape","phyr","dplyr","car","DHARMa","sjPlot","corrplot","ggcorrplot","ggplot2",
#                   "lme4","gridExtra","tidyr","MuMIn","rr2","ggVennDiagram","ggvenn","patchwork","viridis"))

library(ape)
library(phyr)
library(dplyr)
library(car)
library(DHARMa)
library(sjPlot)
library(corrplot)
library(ggcorrplot)
library(ggplot2)
library(lme4)
library(gridExtra)
library(tidyr)
library(MuMIn)
library(rr2)
library(ggVennDiagram)
library(ggvenn)
library(patchwork)
library(viridis)

#####

# Dataset, phylogenetic tree and covariance matrix----

amphib.blood.par <- read.csv("xavier_etal_2025_dataset.csv", header = T)

# Loading the tree file:
phy <- read.nexus("Arvore_16s_full.tree")

# Calculating phylogenetic covariance matrix:
matrix_cov <- vcv(phy)

# Scalating the matrix (Scalating to an unitary variance unit):
matrix_cov_scale <- matrix_cov/max(matrix_cov) 

# Selecting just the species present in the dataset:
matrix_cov <- matrix_cov[rownames(matrix_cov) %in% amphib.blood.par$sp, 
                                     colnames(matrix_cov) %in% amphib.blood.par$sp]

matrix_cov_scale <- matrix_cov_scale[rownames(matrix_cov_scale) %in% amphib.blood.par$sp, 
                                     colnames(matrix_cov_scale) %in% amphib.blood.par$sp]
#####

# Correlation tests----
cor <- cor(amphib.blood.par[,c("precip","tmean","SVL","mass")], use = "complete.obs")
ggcorrplot(cor, method = "circle", type = "lower", lab = TRUE)
corrplot(cor, method = "color", type = "lower", tl.col = "black", tl.srt = 45,
         addCoef.col = "black", number.cex = 0.7)
kruskal.test(mass ~ habit, data = amphib.blood.par)
jpeg(filename = "corrplot.jpg", width = 150, height = 150, units = "mm", res = 600)
dev.off()
#####

# Modelling and model performance analyses----

# Trypanosomatids:
tryp_mod1 <- pglmm(tryp ~ habit + SVL + precip + tmean + (1 | sp__) + (1 | stream), 
                data = amphib.blood.par, 
                cov_ranef = list(sp = matrix_cov), family = "binomial", bayes = F)
summary(tryp_mod1)
residuals.tryp_mod1 <- simulateResiduals(fittedModel = tryp_mod1, n=1000)
plot(residuals.tryp_mod1) #ok
# AIC:
logLik_model1 <- tryp_mod1$logLik #loglikelihood
k_fixef1 <- length(tryp_mod1$fixef) #fixed parameters
k_ranef1 <- length(tryp_mod1$ranef) #random parameters
k1 <- k_fixef1 + k_ranef1 #number of parameters
AIC_tryp1 <- -2 * logLik_model1 + 2 * k1
AIC_tryp1 #89.7 (best model)
r.squaredGLMM(tryp_mod1)

# Rickettsia:
rick_mod1 <- pglmm(rick ~ habit + SVL + precip + tmean + (1 | sp__) + (1 | stream), 
                   data = amphib.blood.par, 
                   cov_ranef = list(sp = matrix_cov), family = "binomial", bayes = F)
summary(rick_mod1)
residuals.rick_mod1 <- simulateResiduals(fittedModel = rick_mod1, n=1000)
plot(residuals.rick_mod1) #ok
# AIC:
logLik_model2 <- rick_mod1$logLik #loglikelihood
k_fixef2 <- length(rick_mod1$fixef) #fixed parameters
k_ranef2 <- length(rick_mod1$ranef) #random parameters
k2 <- k_fixef2 + k_ranef2 #number of parameters
AIC_rick1 <- -2 * logLik_model2 + 2 * k2
AIC_rick1 #106.2 (best model)

# Hepatozoon:
hepa_mod1 <- pglmm(hepa ~ habit + SVL + precip + tmean + (1 | sp__) + (1 | stream), 
                   data = amphib.blood.par, 
                   cov_ranef = list(sp = matrix_cov), family = "binomial", bayes = F)
summary(hepa_mod1)
residuals.hepa_mod1 <- simulateResiduals(fittedModel = hepa_mod1, n=1000)
plot(residuals.hepa_mod1) #ok
# AIC:
logLik_model3 <- hepa_mod1$logLik #loglikelihood
k_fixef3 <- length(hepa_mod1$fixef) #fixed parameters
k_ranef3 <- length(hepa_mod1$ranef) #random parameters
k3 <- k_fixef3 + k_ranef3 #number of parameters
AIC_hepa1 <- -2 * logLik_model3 + 2 * k3
AIC_hepa1 #69.4 (best model)

# pseudo-R2:
# Tryp:
tryp_R2_marg_cond <- R2(tryp_mod1)
print(tryp_R2_marg_cond)

# Rick:
rick_R2_marg_cond <- R2(rick_mod1)
print(rick_R2_marg_cond)

# Hepa:
hepa_R2_marg_cond <- R2(hepa_mod1)
print(hepa_R2_marg_cond)

#####

# Positives per stream ----

counting_stream <- amphib.blood.par %>%
  group_by(stream) %>%
  summarise(frequencia = n())
print(counting_stream)

tryp_pos <- amphib.blood.par %>%
  group_by(stream) %>%
  summarise(positivos = sum(tryp == 1))
print(tryp_pos)

rick_pos <- amphib.blood.par %>%
  group_by(stream) %>%
  summarise(positivos = sum(rick == 1))
print(rick_pos)

hepa_pos <- amphib.blood.par %>%
  group_by(stream) %>%
  summarise(positivos = sum(hepa == 1))
print(hepa_pos)
#####

# Graphic representation----

# Trypanosomatids:

glmm.tryp1 <- glmer(tryp ~ habit + SVL + precip + tmean + (1|stream), 
                    family = binomial, data = amphib.blood.par)
jpeg(filename = "tryp.jpg", width = 220, height = 140, units = "mm", res = 600)

p.tryp <- plot_model(glmm.tryp1, type = "pred", title = "A)",
                     axis.title = c("Daily average temperature (ºC)", "Trypanosomatidae occurrence (%)"), ci.lvl = 0.95,
                     show.data = T, show.values = T, colors = "orange3", show.intercept = T, dot.size = 5,
                     line.size = 1.5, terms = "tmean")

plot.tryp <- p.tryp +
  theme_sjplot2(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 14),     
    axis.text.y = element_text(size = 14),     
    axis.title.x = element_text(size = 16),    
    axis.title.y = element_text(size = 16)     
  )

plot(plot.tryp)  

dev.off()


# Hepatozoon:

glmm.hepa1 <- glmer(hepa ~ habit + SVL + precip + tmean + (1|stream), 
                    family = binomial, data = amphib.blood.par)
jpeg(filename = "hepa.jpg", width = 220, height = 140, units = "mm", res = 600)

p.hepa <- plot_model(glmm.hepa1, type = "pred", title = "B)",
                     axis.title = c("Daily average temperature (ºC)", "Hepatozoon occurrence (%)"), ci.lvl = 0.95,
                     show.data = T, show.values = T, colors = "purple", show.intercept = T, dot.size = 5,
                     line.size = 1.5, terms = "tmean")

plot.hepa <- p.hepa +
  theme_sjplot2(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 14),     
    axis.text.y = element_text(size = 14),     
    axis.title.x = element_text(size = 16),   
    axis.title.y = element_text(size = 16)     
  )

plot(plot.hepa)  

dev.off()


# Rickettsia:

glmm.rick1 <- glmer(rick ~ habit + SVL + precip + tmean + (1|stream), 
                    family = binomial, data = amphib.blood.par)
habit_colors <- c("Arboreal" = "green3", "Rheophilic" = "blue3", "Terrestrial" = "red3")

jpeg(filename = "rick.jpg", width = 220, height = 140, units = "mm", res = 600)

p.rick <- plot_model(glmm.rick1, type = "pred", title = "",
                     axis.title = c("Ecomorphs", "Rickettsia occurrence (%)"), ci.lvl = 0.95,
                     show.data = F, show.values = T, colors = "darkgrey", show.intercept = T, dot.size = 5,
                     line.size = 1.2, terms = "habit")

plot.rick <- p.rick +
  theme_sjplot2(base_size = 12) +
  theme(
    axis.text.x = element_text(color = habit_colors[levels(amphib.blood.par$habit)], 
                               size = 14), 
    axis.title.y = element_text(size = 16) 
  )

plot(plot.rick)  

dev.off()


jpeg(filename = "tryp_hepa.jpg", width = 250, height = 350, units = "mm", res = 600) #graphs combined
grid.arrange(plot.tryp,plot.hepa, ncol = 1, nrow = 2)
dev.off()

#####

# Venn Diagrams - Co-infections----

preparar_lista_venn <- function(df) {   # Filter the lines where there is (1) presence of parasite and extract the sample ID
  list(
    Tryp = df$ZUFABC[df$tryp == 1],
    Rick = df$ZUFABC[df$rick == 1],
    Hepa = df$ZUFABC[df$hepa == 1]
  )
}

# Filter data and prepare lists for each habit from the dataframe 'amphib.blood.par'
dados_arboreal <- subset(amphib.blood.par, habit == "Arboreal")
lista_arboreal <- preparar_lista_venn(dados_arboreal)

dados_rheophilic <- subset(amphib.blood.par, habit == "Rheophilic")
lista_rheophilic <- preparar_lista_venn(dados_rheophilic)

dados_terrestrial <- subset(amphib.blood.par, habit == "Terrestrial")
lista_terrestrial <- preparar_lista_venn(dados_terrestrial)

# Venn Diagram - Arboreal
p_arboreal <- ggVennDiagram(
  lista_arboreal,
  category.names = c("Tryp", "Rick", "Hepa"),
  label_alpha = 0,
  set_size = 4.5,
  label_color = "black"
) +
  scale_fill_gradient(low = "lightgreen", high = "darkgreen") +
  labs(
    title = "Arboreal",
    fill = "N"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  )

print(p_arboreal)

# Venn Diagram - Rheophilic
p_rheophilic <- ggVennDiagram(
  lista_rheophilic,
  category.names = c("Tryp", "Rick", "Hepa"),
  label_alpha = 0,
  set_size = 4.5,
  label_color = "black"
) +
  scale_fill_gradient(low = "lightblue", high = "blue3") + 
  labs(
    title = "Rheophilic",
    fill = "N"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  )

print(p_rheophilic)

# Venn Diagram - Terrestrial
p_terrestrial <- ggVennDiagram(
  lista_terrestrial,
  category.names = c("Tryp", "Rick", "Hepa"),
  label_alpha = 0,
  set_size = 4.5,
  label_color = "black"
) +
  scale_fill_gradient(low = "lightcoral", high = "darkred") + # <-- MUDANÇA AQUI
  labs(
    title = "Terrestrial",
    fill = "N"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  )

print(p_terrestrial)

# Combine the diagrams ensuring the same size for all

first_line <- p_arboreal + p_rheophilic
second_line <- plot_spacer() + p_terrestrial + plot_spacer() +
  plot_layout(widths = c(0.55, 1, 0.45)) 

combined_image <- first_line / second_line

# 5. Exportar a figura final com o layout corrigido
ggsave("venn_diagrams_habits.png",
       plot = combined_image,
       device = "png",
       width = 10,
       height = 8,
       units = "in",
       dpi = 600)
 #####
