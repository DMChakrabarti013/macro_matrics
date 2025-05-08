## continuous treatment effects

# continuous data table

datatable <- read_excel("/Users/dwaipayanchakrabarti/Downloads/cont_datatable.xlsx")

pdf <- pdata.frame(datatable, index = c("State", "YearQuarter"))

install.packages("fixest")
library(fixest)

# additive model
mod_fixest <- feols(
  Y_w ~ M_w + X1_w + X2_w + X3_w 
  | State + YearQuarter, 
  data = pdf
)

summary(mod_fixest, vcov = "HC1") 
summary(mod_fixest, cluster = ~State)

# lag model
mod_fixest2 <- feols(
  Y_w ~ M_w + l(M_w, 1:4) + X1_w + X2_w + X3_w 
  | State + YearQuarter,
  panel.id = c("State", "YearQuarter"),
  data = pdf
)
summary(mod_fixest2, cluster = ~State)
summary(mod_fixest2, vcov = "HC1")
