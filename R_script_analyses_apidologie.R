# Versão 19/03/2025
# Bibliotecas -------------------------------------------------------------

library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(car)
library(lme4)
library(lmerTest)
library(emmeans)
library(DHARMa)
library(dunn.test)
library(glmmTMB)
library(DHARMa)
library(ggridges)
library(patchwork)
library(effectsize)
library(fitdistrplus) 
library(sjPlot)
library(performance)
library(glossary)
library(introdataviz)
library(ggeffects)
library(see)
library(remotes)

# Correção da tabela ------------------------------------------------------
dados_aquisicao <- read_csv("G:/Meu Drive/Mestrado biodiversidade animal/Análises estatísticas/Analises_mestrado/analise_agrotoxicos4/Oficial Dados Experimentais - Mestrado - agrotoxicos_aquisicao.csv")


dados_aquisicao$Group <- factor(dados_aquisicao$Group, levels = (c('controle','delta','imida','imida_delta')))
dados_aquisicao$Group <- recode_factor(dados_aquisicao$Group,
                                       "controle" = "Control",
                                       "delta" = "DEL",
                                       "imida" = "IMI",
                                       "imida_delta" = "COM")


dados_aquisicao$Trial <- factor(dados_aquisicao$Trial, ordered = T)
dados_aquisicao$bloco <- as.factor(dados_aquisicao$bloco)
dados_aquisicao$individuo <- as.factor(dados_aquisicao$individuo)
str (dados_aquisicao)


# Adicição da coluna aprendeu
dados_aquisicao <- dados_aquisicao %>%
  group_by(individuo) %>%
  mutate(aprendeu = if_else(sum(resposta) > 0 | memoria_1h == 1, 1, 0)) %>%
  ungroup()

# GLMM a partir dos dados agrupados, exclui medidas repetidas e trials ------------------------------------
agrup_aprend <- dados_aquisicao %>% 
  group_by (Group, bloco, individuo) %>% 
  summarise (aprendeu = mean (aprendeu)) 


# Seleção de modelos ------------------------------------------------------
modelo_glmm <- glmmTMB(aprendeu ~ Group + 
                         (1| bloco),
                       family = binomial,
                       data = agrup_aprend)

modelo_glmm.1 <- glmmTMB(aprendeu ~ Group,
                       family = binomial,
                       data = agrup_aprend)

# anova pacote Car - Seleção por AIC
anova(modelo_glmm, modelo_glmm.1, test = 'Chisq') # Será utilizado o modelo modelo_glmm.1 por ser mais simples

#Visualização do modelo
plot(simulacao_residuos <- simulateResiduals(fittedModel = modelo_glmm.1))

summary(modelo_glmm.1)

# Extração de tamanhos de efeito e intervalos de confiança ----------------
predic_prop <-ggpredict(modelo_glmm.1, type ="fixed" , terms = c("Group"), ci_level = 0.95)
predic_prop$predicted
predic_prop$conf.low
predic_prop$conf.high 
predic_prop$std.error

# Precisei colocar essa função pra extrair corretamente os SE, porque estava o SE do coef original
predic_prop <- predic_prop %>%
  mutate(
    std.error_trans = std.error * (predicted * (1 - predicted)))%>%
  na.omit ()

tab_model(modelo_glmm.1, transform = "exp", digits = 4)




# Plot - Curva de aquisição -----------------------------------------------
predic_prop <- predic_prop %>% group_by (x) %>% mutate (porcentagem = 1)
# média da resposta agrupada por Trials e Groups
resposta_agregada <- dados_aquisicao %>%
  group_by( Trial,Group) %>%
  summarise( media = mean(resposta, na.rm=T),
             desvio_padrao = sd(resposta, na.rm=T),
             n = n(),
             erro_padrao = desvio_padrao / sqrt (n))
print(resposta_agregada)

  

# Transformação da coluna em numérica
resposta_agregada$Trial_numeric <- as.numeric (resposta_agregada$Trial)

## GGplot - Curva aquisição

plot_curva_aquisicao <- ggplot(resposta_agregada, aes(x = Trial_numeric, y = media,shape = Group, fill = Group,))+
  geom_line(aes(linetype = Group), size = 0.5, alpha = 1, position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = media - erro_padrao, ymax = media+erro_padrao),
                size = 0.5,
                width = 0.3,
                alpha = 1,
                position = position_dodge(width = 0.3))+
  geom_point(size = 4, alpha = 1, colour = "black", position = position_dodge(width = 0.3) ) +
  scale_fill_manual(values = c("#fff", "#bbb", "#bbb", "#444")) +
  scale_shape_manual(values = c(21:25)) +
  scale_linetype_manual(values = c(1,3,5,10))+
  guides(color = guide_legend(override.aes = list(size = 3, shape = NA)),
    linetype = guide_legend(override.aes = list(size = 2.5))) +
  labs(title = "",
       x = "Trial",
       y = "Conditioned PER",
       color = NULL)+
  # guides(
  #   color = guide_legend(position = 'top'),  # Mover a legenda de cor
  #   shape = guide_legend(position = 'top')) +
  # theme_minimal(base_size = 13) +
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format())+
  theme_classic()+
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 12, color = 'black'),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = 'black'),
    axis.title.x = element_text(face = "bold", color = "black", size = 12),
    axis.title.y = element_text(face = "bold", color = "black", size = 12),
    axis.text.x = element_text(face = 'bold', angle = 0, hjust = 0.5, size = 12, color = 'black'),
    axis.text.y = element_text(face = 'bold', size = 12, color = 'black'),
    legend.position = c(0.15,0.87),
    legend.title = element_blank(),
    legend.text = element_text(size = 12, color = 'black'),
    panel.border = element_blank(),  # Remove todas as bordas
    panel.grid.major = element_line(color = "white", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    axis.line.y.left = element_line(color = 'black'),  
    axis.line.x = element_line(color = 'black'),  
    axis.line.y = element_line(color = 'black'),
    axis.line.x.top = element_line(color = "white"))




# Plot de barras --------------------------------------------------------

# Cria um agrupamento para cada individuo com sua média de respostas nos 5 Trials
agrup_aqui <- dados_aquisicao%>%
  group_by(Group,individuo, bloco)%>%
  summarise( resposta = mean(resposta),
             aprendeu = mean(aprendeu))

# A partir da média da resposta pode-se criar uma relação de desvio da média geral do grupo
agrup_aqui_aprendeu <- agrup_aqui%>%
  group_by (Group)%>%
  summarise( media_aprendeu = mean(aprendeu),
             n = n(),
             desvio_padrao = sd(aprendeu),
             erro_padrao = desvio_padrao/ sqrt(n),
             total_aprendeu = sum( aprendeu ==1),
             porcentagem = n()/n())

# Grafico de barras proporção treinados x condicionados -------------------
library(ggpattern)

plot_proporcao <-  ggplot(predic_prop, aes(x = x)) +
  geom_bar( stat = 'identity', aes( y = porcentagem), 
            fill = 'white',
            color = 'black', 
            alpha = 1)+
  geom_text(data = agrup_aqui_aprendeu, aes(x = Group, y = porcentagem, label = n),
            vjust = -0.3,
            color = "black", 
            size = 4.5,
            inherit.aes = F) +
  geom_bar (stat = 'identity', aes(y = predicted, fill = x, x = x),
            color = 'black')+
   # geom_bar_pattern(stat = 'identity',
   #                  aes(y = predicted, fill = x, x = x, pattern = x),
   #                  color = 'black',
   #   pattern_fill = "black",
   #   pattern_density = 0.05,
   #   pattern_spacing = 0.05,
   #   alpha = 1
   # ) +
   # scale_pattern_manual(values = c('circle','stripe','stripe','crosshatch'))+
  geom_text(aes(y = 0, 
                label = total_aprendeu, x = Group),
            vjust = -0.4, 
            color = "black", 
            size = 4.5,
            data = agrup_aqui_aprendeu,
            inherit.aes = F) +
  geom_errorbar(data = predic_prop, aes(ymin = conf.low, ymax = conf.high, x = x),
                linewidth = 0.7, 
                width = 0.2, 
                alpha = 1,
                inherit.aes = F)+
  geom_text(data = predic_prop, aes(y = conf.low,
                                    x = x,
                label = c('a','b','b','b')),
            vjust = 1.5, 
            color = "black", 
            size = 5,
            inherit.aes = F) +
  scale_fill_manual(values = c("#fff", "#bbb", "#bbb", "#777")) +
  labs(title = "",
       x = '',
       y = "Proportion of bees (%)",
       color = NULL )+
  theme_minimal(base_size = 12) +
  theme_classic()+
  theme(
    # text = element_text(family = "Comic Sans MS"),
    # panel.background = element_rect(fill = "white"),
    plot.title = element_text(hjust = 0.5, face = "plain", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title.x = element_text(face = "bold", color = "black"),
    axis.title.y = element_text(face = "bold", color = "black", size = 12),
    axis.text.x = element_text(face = 'bold', angle = 0, hjust = 0.5, size = 12, color = 'black'),
    axis.text.y = element_text(face = 'bold' , size = 12, color = 'black'),
    legend.position = 'NONE',
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),  # Remove todas as bordas
    panel.grid.major = element_line(color = "white", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    axis.line.y.left = element_line(color = 'black'),  # Remove a borda do eixo y esquerda
    axis.line.x = element_line(color = 'black'),  # Remove a borda do eixo x inferior
    axis.line.y = element_line(color = 'red'),  # Adiciona borda ao eixo y direita
    axis.line.x.top = element_line(color = "white"))+
  scale_y_continuous( labels = scales::percent_format())+
  scale_x_discrete(labels = c(('Control'), ('DEL'),('IMI'),('COM ')))


# Acquisition Rate --------------------------------------------------------
#Agrupamento de dados por individuos
dados_velocidade <- dados_aquisicao%>%
  group_by(individuo, Group,bloco)%>%
  summarise(score = sum(resposta),
            .groups = 'drop')%>%
  filter(score !=0)%>%
  mutate (qual_aprendeu = 6-score)
str(dados_velocidade)

tapply(dados_velocidade$qual_aprendeu, dados_velocidade$Group, mean)
tapply(dados_aquisicao$aprendeu, dados_aquisicao$Group, mean)


descdist(dados_velocidade$qual_aprendeu, discrete = T, boot = 1000)

# Seleção de Modelos
modelo_glmm_ratio <- glmmTMB(qual_aprendeu ~ Group + 
                         (1| bloco),
                       family = Gamma (link = 'log'),
                       data = dados_velocidade)

modelo_glmm_ratio.1 <- glmmTMB(qual_aprendeu ~ Group,
                             family = Gamma (link = 'log'),
                             data = dados_velocidade)
# anova (Car) seleção por AIC e chisq
anova(modelo_glmm_ratio, modelo_glmm_ratio.1)

dados_velocidade$Group <- factor(dados_velocidade$Group, levels = c('Control','DEL','IMI','COM'))
summary(modelo_glmm_ratio.1)
exp(1.40691) -(1.96 * 0.08466) 
(1.40691) +(1.96 * 0.08466) 

simulacao_residuos <- simulateResiduals(fittedModel = modelo_glmm_ratio.1)
plot(simulacao_residuos)
check_model (modelo_glmm_ratio.1)

# Tamanhos de efeitos
  standardize_parameters(modelo_glmm_ratio.1)
  tab_model(modelo_glmm_ratio.1,  digits = 4, show.ci = 0.95)
# Coeficientes e intervalos de confiança ----------------------------------

#Intervalos de confiança
predic_inter <-ggpredict(modelo_glmm_ratio.1, type ="fixed" , terms = c("Group"), ci_level = 0.95)
predic_inter$predicted
predic_inter$conf.low
predic_inter$conf.high
predic_inter$std.error


library(emmeans)
# Extrair previsões para os níveis de 'Group'
emm <- emmeans(modelo_glmm_ratio.1, ~ Group)

# Visualizar as previsões com ICs
summary(emm)

# Comparar grupos (opcional)
pairs(emm)

library(effects)
# Extrair previsões para o modelo
eff <- effect("Group", modelo_glmm_ratio.1)
# Visualizar as previsões com ICs
summary(eff)
plot(eff)

library(marginaleffects)
# Extrair previsões para os níveis de 'Group'
pred <- predictions(modelo_glmm_ratio.1, newdata = datagrid(Group = unique))
# Visualizar as previsões com ICs
summary(pred)



#Plotagem em barras
ggplot(velocidade_agrup, aes(x = Group, y = media, fill = Group)) +
  geom_bar (stat = 'identity', , color = 'black')+
  geom_errorbar(aes(ymin = media - erro_padrao, ymax = media +erro_padrao),
                size = 0.8, width = 0.2, alpha = 0.8)+
  geom_text(data = agrup_aqui_aprendeu, aes(y = 0.2, 
                                            label = total_aprendeu),
            vjust = 0.2, 
            color = "black", 
            size = 5) +
  scale_fill_manual(values = c("#fff", "#bbb", "#bbb", "#444")) +
  labs(title = "Means",
       x = '',
       y = "Trials",
       color = NULL )+
  theme_minimal(base_size = 12) +
  theme_classic()+
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(face = "plain", color = "black"),
    axis.title.y = element_text(face = "plain", color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, color = 'black', face = 'bold'),
    axis.text.y = element_text(size = 12, color = 'black', face = "plain"),
    legend.position = 'NONE',
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.border = element_blank(),  # Remove todas as bordas
    panel.grid.major = element_line(color = "white", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    axis.line.y.left = element_line(color = 'black'),  # Remove a borda do eixo y esquerda
    axis.line.x = element_line(color = 'black'),  # Remove a borda do eixo x inferior
    axis.line.y = element_line(color = 'red'),  # Adiciona borda ao eixo y direita
    axis.line.x.top = element_line(color = "white"))+
  scale_x_discrete(labels = c(('Control'), ('Deltamethrin'),('Imidacloprid'),('Combination')))+
  coord_flip ()


###### Grafico de Densidade
plot_densidade <- ggplot(dados_velocidade, aes(x = qual_aprendeu, y = Group, fill = Group)) +
  geom_density_ridges(alpha = 0.9,
                      rel_min_height = 0,
                      scale = 1,
                      color = "black",
                      linewidth = 0.7) +
  scale_fill_manual(values = c("#D3D3D3", "#767F8B", "#B3B7B8", "#5C6068")) +
  labs(title = "Distribuition",
       x = 'Trials',
       y = NULL,
       color = NULL )+
  scale_x_continuous( breaks = 2:5)+
  coord_cartesian(ylim = c(1.6, 4.5),
                  xlim = c(1,6))+
  geom_vline(xintercept = c(2, 3, 4, 5), color = "black", linetype = "dashed") +
  theme_minimal(base_size = 12) +
  theme_classic()+
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(face = "plain", color = "black"),
    axis.title.y = element_text(face = "bold", color = "red"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, color = 'black'),
    axis.text.y = element_blank(),
    legend.position = 'NONE',
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.border = element_blank(),  # Remove todas as bordas
    panel.grid.major = element_line(color = "white", linetype = "dashed"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    axis.line.y.left = element_line(color = 'black'),  # Remove a borda do eixo y esquerda
    axis.line.x = element_line(color = 'black'),  
    axis.line.y = element_line(color = 'red'),  
    axis.line.x.top = element_line(color = "white"))
#Intervalos de confiança

predic_inter <-ggpredict(modelo_glmm_ratio.1, type ="fixed" , terms = c("Group"), ci_level = 0.95)
predic_inter$predicted
predic_inter$conf.low
predic_inter$conf.high
exp(predic_inter$std.error)
summary(modelo_glmm_ratio.1)
exp(0.0938)
plot(predic_inter)

velocidade_agrup <- dados_velocidade %>%
  group_by(Group) %>%
  summarise( media = mean(qual_aprendeu, na.rm=T),
             desvio_padrao = sd(qual_aprendeu, na.rm=T),
             n = n(),
             erro_padrao = desvio_padrao / sqrt (n))

# Graficos combinados

ggplot(data = dados_velocidade, aes(x = Group, y =qual_aprendeu , fill = Group)) +
  introdataviz::geom_flat_violin(trim=F, alpha = 1,
                                 position = position_nudge(x = 0),
                                 scale = 'width',
                                 width = 1.5,
                                 adjust = 1)+
  geom_point (data = predic_inter, aes(y = predicted, x = x, fill = x),
              stat = 'identity',
              color = 'black',
              inherit.aes = FALSE,
              position = position_nudge(x = -0.14),
              size = 3)+
  geom_errorbar(data = predic_inter, aes(ymin = conf.low, ymax = conf.high, x = x),
                stat = 'identity',
                size = 0.8,
                width = 0.2,
                alpha = 0.8,
                inherit.aes = FALSE,
                position = position_nudge(x = -0.14))+
  geom_text(data = predic_inter, aes(y = conf.high,
                                     x = x,
            label = c('b','ab','a','ab')),
            position = position_nudge(x = -0.14),
            vjust = -0.5, 
            color = "black",
            inherit.aes = FALSE,
            size = 5) +
  scale_fill_manual(values = c("#FFF", "#bbb", "#bbb", "#777")) +
  labs(title = "Olfactory Learning Rate",
       x = '',
       y = 'Trial',
       color = NULL )+
  scale_y_continuous( breaks = 2:5)+
  theme_minimal(base_size = 12) +
  theme_classic()+
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(face = "plain", color = "black"),
    axis.title.y = element_text(face = "plain", color = "black", size = 14),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, color = 'black'),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = 12, color = 'black'),
    legend.position = 'NONE',
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.border = element_blank(),  # Remove todas as bordas
    panel.grid.major = element_blank(),  # Remove todas as linhas principais (temporário)
    panel.grid.minor = element_blank(),  # Remove todas as linhas menores
    panel.grid.major.y = element_line(color = "gray", linetype = "dashed"),  # Mostra apenas as linhas horizontais principais
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    axis.line.y.left = element_line(color = 'black'),  # Linha do eixo Y
    axis.line.x = element_line(color = 'black'),  # Linha do eixo X
    axis.line.x.top = element_blank(),  # Remove a linha do eixo X no topo
    axis.line.y = element_line(color = 'black')  # Mantém a linha do eixo Y
  )+
  scale_x_discrete(labels = c(('Control'), ('DEL'),('IMI'),('COM')))


# Memory Retention 1h -----------------------------------------------------

retencao_1h <- dados_aquisicao%>%
  group_by(Group,individuo,bloco)%>%
  summarise(media_aprendeu = mean(aprendeu),
            memoria_1h = mean(memoria_1h))

# Dados filtrados apenas para individuos condicionados
retencao_1h_fil <- retencao_1h%>%
  filter (media_aprendeu != 0)

# Seleção de modelos
glmm_retencao_1h <- glmmTMB(memoria_1h ~ Group +
                              (1| bloco),
                            family = binomial,
                            data = retencao_1h_fil)

glmm_retencao_1h.1 <- glmmTMB(memoria_1h ~ Group,
                            family = binomial,
                            data = retencao_1h_fil)
#Selecao anova (car) AIC e chisq
anova(glmm_retencao_1h,glmm_retencao_1h.1)

#Modelo mais simples é mais adequado
summary(glmm_retencao_1h.1)

#Extração de tamanhos de efeito
tab_model(glmm_retencao_1h.1, transform = "exp", digits = 4)


# Validação do modelo de retenção 1h
simulacao_residuos <- simulateResiduals(fittedModel = glmm_retencao_1h)
plot(simulacao_residuos)

# Extração de intervalos de confiança e predição
predic_1h <-ggpredict(glmm_retencao_1h.1, type ="fixed" , terms = c("Group"), ci_level = 0.95)
predic_1h$predicted
predic_1h$conf.low
predic_1h$conf.high

predic_1h <- predic_1h %>% group_by (x) %>% mutate (porcentagem = 1)

# Plot barras retenção 1h
retencao_1h_agrup <- retencao_1h_fil %>%
  group_by(Group)%>%
  summarise(   media_lembrou_1h = mean(memoria_1h),
               n_aprendeu = n(),
               n_lembrou_1h = sum(memoria_1h),
               sd_lembrou_1h = sd(memoria_1h),
               se_lembrou_1h = sd_lembrou_1h/sqrt(n_aprendeu),
               porcentagem_1h = round(n_lembrou_1h/n_aprendeu,2),
               porce = 1)


#Plot 
plot_barras_1h <- ggplot(predic_1h, aes(x = x)) +
  geom_bar( stat = 'identity', aes( y = porcentagem, fill = 'white', x=x), 
            color = 'black', 
            alpha = 1)+
  geom_bar( stat = 'identity', aes( y = predicted,fill = x), 
            color = 'black', 
            alpha = 1)+
  geom_text(data = retencao_1h_agrup, aes(y = porce,
                label = n_aprendeu, x= Group),
            vjust = -0.3,
            # hjust = 0,
            color = 'black', 
            size = 5,
            inherit.aes = F) +
  geom_text(data = retencao_1h_agrup,aes(y = 0, 
                label = n_lembrou_1h, x = Group),
            vjust = -0.4, 
            hjust = 0,
            color = "black", 
            size = 5,
            inherit.aes = F) +
  geom_errorbar( aes(ymin = conf.low, ymax = conf.high, x = x),
                stat = 'identity',
                size = 0.8,
                width = 0.2,
                alpha = 0.8,
                inherit.aes = FALSE,
                position = position_nudge(x = -0.07))+
  scale_fill_manual(values = c("#fff", "#bbb", "#bbb", "#777",'#fff')) +
  labs(title = "1 Hour",
       x = '',
       y = "Proportion of bees (%)",
       color = NULL )+
  theme_minimal(base_size = 12) +
  theme_classic()+
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(face = "plain", color = "black"),
    axis.title.y = element_text(face = "plain", color = "black", size = 14),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, color = 'black'),
    axis.text.y = element_text(size = 12, color = 'black'),
    legend.position = 'NONE',
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.border = element_blank(),  # Remove todas as bordas
    panel.grid.major = element_line(color = "white", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    axis.line.y.left = element_line(color = 'black'),  # Remove a borda do eixo y esquerda
    axis.line.x = element_line(color = 'black'),  # Remove a borda do eixo x inferior
    axis.line.y = element_line(color = 'red'),  # Adiciona borda ao eixo y direita
    axis.line.x.top = element_line(color = "white"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c(('Control'), ('DEL'),('IMI'),('COM')))

# Importação dados de extinção --------------------------------------------


dados_extincao <- read_csv("Oficial Dados Experimentais - Mestrado - agrotoxicos_extincao.csv")

dados_extincao$Group <- factor(dados_extincao$Group, levels = (c("controle","delta","imida","imida_delta")))
dados_extincao$Trial <- factor(dados_extincao$Trial, ordered = T)
dados_extincao$bloco <- as.factor(dados_extincao$bloco)
dados_extincao$individuo <- as.factor(dados_extincao$individuo)
str(dados_aquisicao)
  
# GLMM memoria 24h --------------------------------------------------------

# Agrupamento da tabela de extinção e filtragem
retencao_24h <- dados_extincao%>%
  group_by(Group,individuo,bloco,aprendeu)%>%
  filter (aprendeu == 1)%>%
  summarise(lembrou_24h = mean(memoria_24h))
str(retencao_24h)

## Modelo a filtrado a partir dos que aprenderam
glmm_retencao_24h <- glmmTMB(lembrou_24h ~ Group +
                               (1| bloco),
                             family = binomial,
                             data = retencao_24h)

glmm_retencao_24h.1 <- glmmTMB(lembrou_24h ~ Group,
                             family = binomial,
                             data = retencao_24h)
#seleção de modelos anova (Car)
anova(glmm_retencao_24h, glmm_retencao_24h.1)
Anova(glmm_retencao_24h)
#extração de dados
summary (glmm_retencao_24h)

plot(simulacao_residuos <- simulateResiduals(fittedModel = glmm_retencao_24h))

# Extração de tamanhos de efeito
tab_model(glmm_retencao_24h, transform = "exp", digits = 4)
tapply(retencao_24h$lembrou_24h, retencao_24h$Group, mean)

predic_24h <-ggpredict(glmm_retencao_24h, type ="fixed" , terms = c("Group"), ci_level = 0.95)
predic_24h$predicted
predic_24h$conf.low
predic_24h$conf.high

predic_24h <- predic_24h %>% group_by (x) %>% mutate (porcentagem = 1)

library(emmeans)
emmeans(glmm_retencao_24h, pairwise ~ Group, type = "response")


# Plot memoria 24h --------------------------------------------------------


# Agrupamento dos dados para criação dos graficos

retencao_24h_agrup <- dados_extincao %>%
  filter(aprendeu == 1)%>%
  group_by(Group)%>%
  summarise(   media_memoria_24h = mean(memoria_24h, na.rm = TRUE),
               n_aprendeu = n()/10,
               n_memoria_24h = sum(memoria_24h, na.rm = TRUE)/10,
               sd_memoria_24h = sd(memoria_24h, na.rm = TRUE),
               se_memoria_24h = sd_memoria_24h/sqrt(n_aprendeu),
               porcentagem_24h = round(n_memoria_24h/n_aprendeu,2),
               porce = 1)

#Plot 
plot_barras_24h <- ggplot(predic_24h, aes(x = x)) +
  geom_bar( stat = 'identity', aes( y = porce,  x = Group), 
            fill = 'white',
            color = 'black', 
            alpha = 1,
            data = retencao_24h_agrup,
            inherit.aes = F)+
  geom_bar( stat = 'identity', aes( y = predicted,fill = x), 
            color = 'black', 
            alpha = 1)+
  geom_text(data = retencao_24h_agrup, aes(y = porce,
                label = n_aprendeu, x = Group),
            vjust = -0.3,
            # hjust = -0.3,
            color = 'black', 
            size = 5,
            inherit.aes = F) +
  geom_text(aes(y = 0, 
                label = n_memoria_24h, x= Group),
            vjust = -0.4, 
            hjust = -0.3,
            color = "black", 
            size = 5,
            data = retencao_24h_agrup,
            inherit.aes = F) +
  geom_errorbar( aes(ymin = conf.low, ymax = conf.high, x = x),
                 stat = 'identity',
                 size = 0.8,
                 width = 0.2,
                 alpha = 0.8,
                 inherit.aes = FALSE,
                 position = position_nudge(x = -0.07))+
  geom_text(data = predic_24h, aes(y = conf.high,
                                     x = x,
                                     label = c('a','a','b','b')),
            position = position_nudge(x = -0.07),
            vjust = -0.5, 
            color = "black",
            inherit.aes = FALSE,
            size = 5) +
  scale_fill_manual(values = c("#fff", "#bbb", "#bbb", "#777")) +
  labs(title = "24 Hours",
       x = '',
       y = "Proportion of bees (%)",
       color = NULL )+
  theme_minimal(base_size = 12) +
  theme_classic()+
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    axis.title.x = element_text(face = "plain", color = "black"),
    axis.title.y = element_text(face = "plain", color = "black", size = 14),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 12, color = 'black'),
    axis.text.y = element_text(size = 12, color = 'black'),
    legend.position = 'NONE',
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.border = element_blank(),  # Remove todas as bordas
    panel.grid.major = element_line(color = "white", linetype = "dashed"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    axis.line.y.left = element_line(color = 'black'),  # Remove a borda do eixo y esquerda
    axis.line.x = element_line(color = 'black'),  # Remove a borda do eixo x inferior
    axis.line.y = element_line(color = 'red'),  # Adiciona borda ao eixo y direita
    axis.line.x.top = element_line(color = "white"))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels = c(('Control'), ('DEL'),('IMI '),('COM ')))

# Graficos combinados -----------------------------------------------------

  #Plot de aquisição e proporção
 plot_proporcao +plot_curva_aquisicao +  plot_annotation(tag_levels = 'a') 

  #Barras e densidade de velocidade
plot_barras_velocidade + plot_densidade +
  plot_annotation(title = "Trials Needed For Acquisition")+
  plot_layout(ncol = 2,nrow = 1, heights = c(1, 1),widths = c(1,1.2), guides = "collect") 

  #Barras de memoria 1h e 24h

plot_barras_1h + plot_barras_24h +
  plot_annotation (title = "", tag_levels = 'a')
# Análise de mortalidade --------------------------------------------------

mortalidade <- dados_extincao%>%
  group_by (Group, individuo) %>%
  summarise (resposta = mean (resposta))%>%
  mutate ( Morta = ifelse (is.na (resposta), 1, 0))
View(mortalidade)

mortalidade_agrup <- mortalidade%>%
  group_by (Group)%>%
  summarise (media_morta = mean (Morta))

anova_mortalidade <- aov(Morta ~ Group, data = mortalidade)
summary(anova_mortalidade)
