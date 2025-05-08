#### G-scores test

# condition 2

# import file

df <- read_csv("/Users/dwaipayanchakrabarti/Downloads/final_data-3.csv")
df2 <- read_csv("/Users/dwaipayanchakrabarti/Downloads/baseline_data.csv")
gdp_g_lm <- lm(df$gdp_growth ~ df$G_score)
summary(gdp_g_lm)

gdp_g_lm_all <- lm(df$gdp_growth ~ df$G_score + df2$low_diff...8 
                   + df2$assetquality_diff + df2$profitability_diff)
summary(gdp_g_lm_all)

plot(df$G_score, df$gdp_growth)

# condition 3 done in excel

# condition 4 - logit regression
library(aod)

df$G_score <- factor(df$G_score)

propensity_model <- glm(df$G_score ~ df2$low_diff...8 + df2$assetquality_diff + df2$profitability_diff, 
                 family = "binomial")

df$propensity_score <- predict(propensity_model, type = "response")

plot_data <- data.frame(
  propensity = df$propensity_score,
  group = ifelse(df$G_score > median(df$G_score), "Treated", "Donor")
)

library(ggplot2)
overlap_plot <- ggplot(plot_data, aes(x = propensity, fill = group)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("Donor" = "blue", "Treated" = "red")) +
  labs(x = "Propensity Score", y = "Density", fill = "Group") +
  theme_minimal() +
  theme(legend.position = "bottom")
