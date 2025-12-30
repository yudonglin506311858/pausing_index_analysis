

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
  mutate(logPI = log2(PI))
head(pi_long)

pi_long <- pi_long %>%
  mutate(group = str_remove(sample, "_sorted_rmdup$"))

head(pi_long)
head(pi_long)
scale_fill_manual(values = c(
  "35" = "#F7F7F7",
  "36" = "#FDD9C4",
  "DPY30–mAID" = "#F4A582"
)
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


