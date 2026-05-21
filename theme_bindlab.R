# theme_bindlab.R — Nature-style ggplot2 publication theme

library(ggplot2)

# ---- 字体设置（Windows中文环境兼容） ----
# 使用sans（系统默认无衬线字体），避免Windows字体注册问题
.bindlab_font <- "sans"

theme_bindlab <- function(base_size = 12, base_family = .bindlab_font) {
  theme_classic(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold",
                                hjust = 0, margin = ggplot2::margin(b = 10)),
      plot.subtitle = element_text(size = base_size, color = "grey30",
                                   hjust = 0, margin = ggplot2::margin(b = 8)),
      axis.title = element_text(size = base_size, face = "bold", color = "black"),
      axis.title.x = element_text(margin = ggplot2::margin(t = 8)),
      axis.title.y = element_text(angle = 90, vjust = 0.5, margin = ggplot2::margin(r = 6)),
      axis.text = element_text(size = base_size - 1, color = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = unit(0.15, "cm"),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text = element_text(size = base_size - 2),
      legend.key.size = unit(0.4, "cm"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.position = "right",
      panel.background = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", color = "black"),
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
}

theme_bindlab_box <- function(base_size = 12, base_family = .bindlab_font) {
  theme_bindlab(base_size, base_family) %+replace%
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.line = element_blank()
    )
}

theme_bindlab_minimal <- function(base_size = 12, base_family = .bindlab_font) {
  theme_bindlab(base_size, base_family) %+replace%
    theme(
      axis.title = element_text(size = base_size - 1, face = "plain", color = "grey30")
    )
}

theme_volcano <- function(base_size = 12) {
  theme_bindlab(base_size) %+replace%
    theme(
      legend.position = c(0.85, 0.85),
      legend.background = element_rect(fill = "white", colour = NA)
    )
}

theme_set(theme_bindlab())

# ---- Treatment标签映射函数 ----
# 在图表中自动将内部"Induced"/"Control"显示为Adi/Ost/Cho等
scale_color_treatment <- function(...) {
  labels <- COLORS$treatment_labels  # c(Induced="Adi", Control="Adi-Control")
  scale_color_manual(values = COLORS$treatment, labels = labels, ...)
}

scale_fill_treatment <- function(...) {
  labels <- COLORS$treatment_labels
  scale_fill_manual(values = COLORS$treatment, labels = labels, ...)
}

cat("[THEME] Nature-style theme loaded (font:", .bindlab_font, ")\n")
