
#########################################################作图#######################################################
head(pi_matrix)
#不同条件下基因 Pausing Index（PI）的 log 转换后的箱线图 + 统计检验（Wilcoxon）。
pi_long <- pi_matrix %>%
  pivot_longer(
    cols = -gene_id,
    names_to = "sample",
    values_to = "PI"
  ) %>%
  filter(PI > 0) %>%                 # 防止 log 出问题
  mutate(logPI = log10(PI))
head(pi_long)

pi_long <- pi_long %>%
  mutate(group = str_remove(sample, "_sorted_rmdup$"))


pi_long$group<-factor(pi_long$group,levels = c("B16-Ser2","P-Ser2","35-Ser2","36-Ser2",
                                               "B16-Ser5","P-Ser5","35-Ser5","36-Ser5"))

p <- ggplot(pi_long, aes(x = group, y = logPI, fill = group)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA
  ) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  # 设置颜色 - 为每个样本单独指定颜色
  scale_fill_manual(values = c(
    "B16-Ser2" = "#F7F7F7",    # 浅灰
    "P-Ser2" = "#E0E0E0",      # 稍深灰
    "35-Ser2" = "#FDD9C4",     # 浅橙色
    "36-Ser2" = "#F4A582",     # 橙色
    "B16-Ser5" = "#D1E5F0",    # 浅蓝色
    "P-Ser5" = "#92C5DE",      # 蓝色
    "35-Ser5" = "#F4C3C2",     # 浅粉色
    "36-Ser5" = "#E78AC3"      # 粉色
  )) +
  labs(
    y = "log2(pausing index)",
    x = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 12)
  )

p


library(tidyverse)
library(ggpubr)

p +
  stat_compare_means(
    comparisons = list(
      c("B16-Ser2", "35-Ser2"),
      c("P-Ser2", "35-Ser2"),
      c("B16-Ser2", "36-Ser2"),
      c("P-Ser2", "36-Ser2"),
      c("B16-Ser5", "35-Ser5"),
      c("P-Ser5", "35-Ser5"),
      c("B16-Ser5", "36-Ser5"),
      c("P-Ser5", "36-Ser5")
    ),
    method = "wilcox.test",
    label = "p.format"
  )



library(extrafont)
font_import()  # 导入系统字体（第一次需要运行）
loadfonts()    # 加载字体
fonts()        # 查看可用字体


ggsave("箱式图-ser2+ser5.pdf",width = 8,height = 6, family = "Arial")



#Pausing index (PI) was calculated as the ratio of RNAPII ChIP–seq signal density at promoter regions (−200 bp to +200 bp relative to the TSS) to that across the gene body (from +400 bp downstream of the TSS to the TES). Cumulative distribution plots of log-transformed PI values were generated separately for RNAPII Ser2- and Ser5-phosphorylated forms.Pausing index (PI) was calculated as the ratio of RNAPII ChIP–seq signal density at promoter regions (−200 bp to +200 bp relative to the TSS) to that across the gene body (from +400 bp downstream of the TSS to the TES). Cumulative distribution plots of log-transformed PI values were generated separately for RNAPII Ser2- and Ser5-phosphorylated forms.Pausing index (PI) was calculated as the ratio of RNAPII ChIP–seq signal density at promoter regions (−200 bp to +200 bp relative to the TSS) to that across the gene body (from +400 bp downstream of the TSS to the TES). Cumulative distribution plots of log-transformed PI values were generated separately for RNAPII Ser2- and Ser5-phosphorylated forms.Pausing index (PI) was calculated as the ratio of RNAPII ChIP–seq signal density at promoter regions (−200 bp to +200 bp relative to the TSS) to that across the gene body (from +400 bp downstream of the TSS to the TES). Cumulative distribution plots of log-transformed PI values were generated separately for RNAPII Ser2- and Ser5-phosphorylated forms.Pausing index (PI) was calculated as the ratio of RNAPII ChIP–seq signal density at promoter regions (−200 bp to +200 bp relative to the TSS) to that across the gene body (from +400 bp downstream of the TSS to the TES). Cumulative distribution plots of log-transformed PI values were generated separately for RNAPII Ser2- and Ser5-phosphorylated forms.Pausing index (PI) was calculated as the ratio of RNAPII ChIP–seq signal density at promoter regions (−200 bp to +200 bp relative to the TSS) to that across the gene body (from +400 bp downstream of the TSS to the TES). Cumulative distribution plots of log-transformed PI values were generated separately for RNAPII Ser2- and Ser5-phosphorylated forms.Pausing index (PI) was calculated as the ratio of RNAPII ChIP–seq signal density at promoter regions (−200 bp to +200 bp relative to the TSS) to that across the gene body (from +400 bp downstream of the TSS to the TES). Cumulative distribution plots of log-transformed PI values were generated separately for RNAPII Ser2- and Ser5-phosphorylated forms.Pausing index (PI) was calculated as the ratio of RNAPII ChIP–seq signal density at promoter regions (−200 bp to +200 bp relative to the TSS) to that across the gene body (from +400 bp downstream of the TSS to the TES). Cumulative distribution plots of log-transformed PI values were generated separately for RNAPII Ser2- and Ser5-phosphorylated forms.

library(ggplot2)
library(dplyr)



colnames(pi_matrix) <-str_remove(colnames(pi_matrix), "_sorted_rmdup$")

pi_colors <- c(
  "B16-Ser2" = "#F7F7F7",    # 浅灰
  "P-Ser2"   = "#E0E0E0",    # 稍深灰
  "35-Ser2"  = "#FDD9C4",    # 浅橙色
  "36-Ser2"  = "#F4A582",    # 橙色
  "B16-Ser5" = "#D1E5F0",    # 浅蓝色
  "P-Ser5"   = "#92C5DE",    # 蓝色
  "35-Ser5"  = "#F4C3C2",    # 浅粉色
  "36-Ser5"  = "#E78AC3"     # 粉色
)

plot_pi_ecdf <- function(pi_matrix, pattern, title = NULL) {
  
  # 选取 Ser2 或 Ser5
  cols <- grep(pattern, colnames(pi_matrix), value = TRUE)
  
  df <- pi_matrix %>%
    dplyr::select(all_of(cols)) %>%
    pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = "pi"
    ) %>%
    filter(is.finite(pi), pi > 0) %>%
    mutate(log_pi = log10(pi))
  
  ggplot(df, aes(x = log_pi, color = sample)) +
    stat_ecdf(size = 1.1) +
    scale_color_manual(values = pi_colors, drop = FALSE) +
    labs(
      x = "log10(pausing index)",
      y = "Fraction of genes",
      title = title
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.title = element_blank(),
      legend.position = "right"
    )
}
p_ser2 <- plot_pi_ecdf(
  pi_matrix,
  pattern = "Ser2",
  title = "RNAPII Ser2 pausing index"
)

print(p_ser2)


ggsave("RNAPII 暂停指数的累计图-ser2.pdf",width = 5,height = 5, family = "Arial")

p_ser5 <- plot_pi_ecdf(
  pi_matrix,
  pattern = "Ser5",
  title = "RNAPII Ser5 pausing index"
)

print(p_ser5)


ggsave("RNAPII 暂停指数的累计图-ser5.pdf",width = 5,height = 5, family = "Arial")

p_ser2+p_ser5

ggsave("RNAPII 暂停指数的累计图-ser2+ser5.pdf",width = 9,height = 5, family = "Arial")
