
# reading KE-DHS data & trying vaccine timeliness on Measles
pacman::p_load(ggpubr, dplyr, haven, ggplot2, janitor, lubridate, stringr, scales,
               patchwork)
mvs <- naniar::miss_var_summary

theme_util <- function(base_size = 12,
                       base_family = "Times New Roman") {

  theme_bw(base_size = base_size, base_family = base_family) +

    theme(

      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),

      panel.border = element_rect(
        colour = "black",
        linewidth = 0.6
      ),

      axis.line = element_blank(),

      axis.title = element_text(
        face = "bold",
        size = base_size + 1
      ),

      axis.text = element_text(
        colour = "black",
        size = base_size
      ),

      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.5
      ),

      plot.title = element_text(
        face = "bold",
        hjust = 0,
        size = base_size + 2
      ),

      plot.subtitle = element_text(
        hjust = 0,
        size = base_size
      ),

      plot.caption = element_text(
        hjust = 1,
        size = base_size - 2
      ),

      strip.background = element_blank(),

      strip.text = element_text(
        face = "bold",
        size = base_size
      ),

      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(face = "bold"),

      plot.margin = margin(8, 8, 8, 8)
    )
}

# data --------------------------------------------------------------------

vax_data <- readRDS("data/processed/vaxdata-components.rds")

use <- c('vita1', 'penta1', 'penta2', 'penta3', 'bcg', 'mcv1', 'mcv2')
for (i in 1:length(use)) {
  message(use[i])

  x <- vax_data[[use[i]]] |> filter(!is.na(vaxx_date)) |>
    mutate(
      vaxx_date = ifelse(vaxx_date < birth_date, (birth_date), (vaxx_date)),
      vaxx_date = vaxx_date + ymd('1970-01-01'),
      dvar = vaxx_date - due_date,
    )
  sq <- seq(from = min(x$dvar), to = max(x$dvar), length = 10)
  sq_use <- log(abs(sq), base = 10) * sign(sq)

  ggplot(
    x |> filter(dvar != 0) |> mutate(dvar = as.integer(dvar))
  ) +
    geom_histogram(
      aes(x = dvar),
      colour = "white",
      fill = "black",
      bins = 100
    ) +
    labs(
      x = "Days relative to the recommended vaccination age",
      y = "Number of children"
    ) +
    theme_util() |

    ggplot(
      x |> filter(dvar != 0) |> mutate(dvar = as.integer(dvar))
    ) +
    geom_histogram(
      aes(dvar),
      bins = 100,
      fill = "black",
      colour = "white"
    ) +
    scale_x_continuous(
      trans = pseudo_log_trans(base = 10),
      breaks = c(-365,-180,-90,-30,-14,-7,-1,
                 1,7,14,30,90,180,365,730)
    ) +
    labs(
      x = "Days relative to the recommended vaccination age",
      y = "Number of children",
      caption = 'Log(base=10) transformation used\nNegative values signify early vaccination, while positive values signify vaccine delay'
    ) +
    theme_util()
  ggsave(filename = paste0('output/img/', use[i], '/vaccine-delay-histogram.png'),
         height = 6, width = 14, dpi = 1000)
}


# -------------------------------------------------------------------------

x <- vax_data$bcg
dim(x)
