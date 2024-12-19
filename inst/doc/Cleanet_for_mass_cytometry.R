## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(Cleanet)
library(readr)
library(dplyr)
library(ggplot2)

## ----read_data----------------------------------------------------------------
path <- system.file("extdata", "df_mdipa.csv", package="Cleanet")
df_mdipa <- read_csv(path, col_types=cols())
print(df_mdipa)

## ----cleanet_basic------------------------------------------------------------
cols <- c("CD45", "CD123", "CD19", "CD11c", "CD16",
          "CD56", "CD294", "CD14", "CD3", "CD20",
          "CD66b", "CD38", "HLA-DR", "CD45RA",
          "DNA1", "DNA2")
cleanet_res <- cleanet(df_mdipa, cols, cofactor=5)

## ----cleanet_basic_status-----------------------------------------------------
print(table(cleanet_res$status))

## ----cleanet_basic_sensitivity------------------------------------------------
print(cleanet_res$sensitivity)

## ----bivariate_basic----------------------------------------------------------
ggplot(df_mdipa, aes(x=asinh(DNA1/5), y=asinh(CD45/5), color=cleanet_res$status)) +
  geom_point(size=0.2) +
  scale_color_discrete(name="Status") +
  theme_bw()

## ----debris_default-----------------------------------------------------------
is_debris <- filter_debris_cytof(df_mdipa, cols)

## ----debris_custom------------------------------------------------------------
is_debris <- filter_debris_cytof(df_mdipa, cols, threshold = 0.35)

## ----cleanet_filtered---------------------------------------------------------
cleanet_res <- cleanet(df_mdipa, cols, cofactor=5, is_debris=is_debris)
print(cleanet_res$sensitivity)

ggplot(df_mdipa, aes(x=asinh(DNA1/5), y=asinh(CD45/5), color=cleanet_res$status)) +
  geom_point(size=0.2) +
  scale_color_discrete(name="Status") +
  theme_bw()

## ----CD294--------------------------------------------------------------------
ggplot(df_mdipa, aes(x=asinh(DNA1/5), y=asinh(CD45/5), color=asinh(CD294/5))) +
  geom_point(size=0.2) +
  scale_color_gradient(low="black", high="red") +
  theme_bw()

## ----label--------------------------------------------------------------------
print(table(df_mdipa$label))

## ----classify_doublets--------------------------------------------------------
singlet_clas <- df_mdipa$label[which(cleanet_res$status!="Doublet")]
doublet_clas <- classify_doublets(cleanet_res, singlet_clas)
sort(table(doublet_clas))

## ----compare------------------------------------------------------------------
df_exp_obs <- compare_doublets_exp_obs(doublet_clas, singlet_clas, cleanet_res)
arrange(df_exp_obs, -Expected)

## ----compare_plot-------------------------------------------------------------
ggplot(df_exp_obs, aes(x=Expected, y=Observed)) +
  geom_point() +
  geom_abline(slope=1, yintercept=0, linetype="dotted") +
  theme_bw()

