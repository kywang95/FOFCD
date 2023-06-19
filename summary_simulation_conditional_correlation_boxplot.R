library(dplyr)
library(ggplot2)
library(showtext)
showtext_auto(enable = TRUE)
font_add("newrom", regular = "C:/Windows/Fonts/times.ttf")

############################################################################
## corelevant simulation
path <- "./"
iteration <- 100
method = c('FSCA.1', 'FCCA1', 'FunNCC1.osb'
           , 'FCCA12', 'FunNCC12.osb'
)
method_name = c('FSCA.DC1', 'FCCA1', 
                'FunNCC1'
                , 'FCCA12', 'FunNCC1|2'
)

pathname <- paste0(path, "/")
data_out <- read.csv("result_partial_correlation_simu_corelevant_corr.csv")
data <- data_out %>% 
  dplyr::select(method) 
result <- data.frame(data = unlist(data), 
                     method = rep(method_name, each = iteration))
result$method <- factor(result$method, levels = method_name, ordered = TRUE)

result_1 <- result[which(result$method %in% method_name[1:3]),]
P1 <- result_1 %>% 
  ggplot(aes(x = method, y = data, fill=method)) +
  stat_boxplot(geom="errorbar",alpha=0.7, width=0.3, position=position_dodge(0.8)) +
  geom_boxplot(alpha=1, outlier.color = "grey50",                 
               outlier.fill = "grey50",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               position=position_dodge(0.8), 
               width=0.6) +
  geom_hline(aes(yintercept = 0), linetype = 5, col = "grey50") + 
  scale_y_continuous(name = NULL,
                     limits = c(-0.25, 1), breaks = seq(-0.25, 1, by = 0.25)
  ) +
  scale_x_discrete(name = NULL, 
                   labels = c(expression(R*'('*'Y'*','*'X'[1]*')'),
                              expression(rho*'('*'Y'*','*'X'[1]*')'),
                              expression(gamma*'('*'Y'*','*'X'[1]*')'))) +
  scale_fill_manual(
    values=c(
      "FCCA1" = "#990066", 
      "FunNCC1" = "#FFCC00", 
      "FSCA.DC1" = "#CC0033"
    )) + 
  theme_bw() +
  theme(panel.background=element_rect(fill='transparent'),  
        panel.grid=element_blank(),                          
        axis.title.y = element_text(size = 18, family="newrom"), 
        axis.title.x = element_text(size = 18, family="newrom"), 
        axis.text.x = element_text(size = 15, family="newrom"), 
        axis.text.y = element_text(size = 15, family="newrom"),
        legend.position ="none"
        )

ggsave(P1, filename = paste0("boxplot_partial_correlation_simu_corelevant_a_230526.pdf"),
       width = 10, height = 12,
       units = "cm",
       dpi = 1200, device = cairo_pdf)
ggsave(P1, filename = paste0("boxplot_partial_correlation_simu_corelevant_a_230526.eps"),
       width = 10, height = 12,
       units = "cm",
       dpi = 1200, device = cairo_ps)

result_2 <- result[which(result$method %in% method_name[4:5]),]
P2 <- result_2 %>% 
  ggplot(aes(x = method, y = data, fill=method)) +
  stat_boxplot(geom="errorbar",alpha=0.7, width=0.3, position=position_dodge(0.8)) +
  geom_boxplot(alpha=1, outlier.color = "grey50",                 
               outlier.fill = "grey50",
               outlier.shape = 19,
               outlier.size = 1,
               outlier.stroke = 0.5,
               outlier.alpha = 45,
               position=position_dodge(0.8), 
               width=0.6) +
  geom_hline(aes(yintercept = 0), linetype = 5, col = "grey50") + 
  scale_y_continuous(name = NULL,
                     limits = c(-0.25, 1), breaks = seq(-0.25, 1, by = 0.25)
  ) +
  scale_x_discrete(name = NULL, 
                   labels = c(expression(gamma*'('*'Y'*','*'X'[1]*'|'*'X'[2]*')'),
                              expression(rho*'('*'Y'*','*'X'[1]*'|'*'X'[2]*')')
                              )) +
  scale_fill_manual(
    values=c("FunNCC1|2" = "#FFCC00",
      'FCCA12' = "#990066"
    )) + 
  theme_bw() +
  theme(panel.background=element_rect(fill='transparent'),  
        panel.grid=element_blank(),                         
        axis.title.y = element_text(size = 18, family="newrom"), 
        axis.title.x = element_text(size = 18, family="newrom"), 
        axis.text.x = element_text(size = 15, family="newrom"), 
        axis.text.y = element_text(size = 15, family="newrom"),
        legend.position ="none"
  )

ggsave(P2, filename = paste0("boxplot_partial_correlation_simu_corelevant_b_230526.pdf"),
       width = 7.5, height = 12,
       units = "cm",
       dpi = 1200, device = cairo_pdf)
ggsave(P2, filename = paste0("boxplot_partial_correlation_simu_corelevant_b_230526.eps"),
       width = 7.5, height = 12,
       units = "cm",
       dpi = 1200, device = cairo_ps)
