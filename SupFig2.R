###############################################################################
# Install Packages
###############################################################################

options(repos = "http://cran.us.r-project.org")
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pacman::p_load(vroom, tidyr, tidyfast, dplyr, ggplot2, ggrepel)

###############################################################################
# Download and format the data
###############################################################################

tbl_raw <- vroom("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3069nnn/GSM3069461/suppl/GSM3069461_SER6_DGE.txt.gz")

tbl_tidy <-
    tbl_raw %>%
    mutate(across(where(is.double), ~ .x / sum(.x) * 10000)) %>%
    mutate(across(where(is.double), ~ log1p(.x))) %>%
    dt_pivot_longer(
        cols = -GENE,
        names_to = "cell",
        values_to = "expression"
    ) %>%
    as_tibble()

tbl_mean_expression <-
    tbl_tidy %>%
    group_by(GENE) %>%
    summarize(mean = log(mean(expression, na.rm = TRUE))) %>%
    mutate(size = if_else(str_detect(GENE, "Exoc[1-9]$"), 2, 1)) %>%
    mutate(label = if_else(str_detect(GENE, "Exoc[1-9]$"), GENE, ""))

###############################################################################
# Violin plot
###############################################################################

ggplot(tbl_mean_expression,
    aes(x = 1, y = mean, color = size, size = size, label = label)) +
    geom_violin() +
    geom_point() +
    geom_text_repel(hjust = 1) +
    labs(x = "SER6", y = "log of average gene expressions", legend = NULL) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
        text = element_text(size = 20),
        legend.position = "none")

ggsave("SER6_Exoc.pdf", width = 6, height = 14)
