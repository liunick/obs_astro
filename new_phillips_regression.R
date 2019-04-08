library(lmtest)
library(sandwich)

b_data <- read.csv("b_band.csv")
v_data <- read.csv("v_band.csv")
i_data <- read.csv("i_band.csv")

new_b_data <- data.frame(peak_mag = b_data$peak_mag, delta_mag = (b_data$X15_mag - b_data$peak_mag))
plot(new_b_data)
b_reg <- lm(new_b_data$peak_mag~new_b_data$delta_mag)
summary(b_reg)

new_v_data <- data.frame(peak_mag = v_data$peak_mag, delta_mag = (v_data$X15_mag - v_data$peak_mag))
plot(new_v_data)
v_reg <- lm(new_v_data$peak_mag~new_v_data$delta_mag)
summary(v_reg)

new_i_data <- data.frame(peak_mag = i_data$peak_mag, delta_mag = (i_data$X15_mag - i_data$peak_mag))
plot(new_i_data)
i_reg <- lm(new_i_data$peak_mag~new_i_data$delta_mag)
summary(i_reg)

