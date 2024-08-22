file_path = "/Users/joshnghe/Desktop/Code/AMS315/project 2/data_analysis_2/P2_908487.csv"
dataTable = read.csv(file_path)

rm(data)

library(MASS)

cor_matrix <- cor(dataTable)
high_cor <- cor_matrix[upper.tri(cor_matrix, diag = TRUE) & abs(cor_matrix) > 0.05]
print(high_cor)

iv_cor <- cor(dataTable)[,"Y"]
print(iv_cor)

library(ggplot2)
y <- dataTable$Y
environmental_variables <- dataTable$E1
plot(y,environmental_variables, col='darkgreen', pch=19)

# Adding scatterplot of gfg_x2 vs gfg_y2
points(y, dataTable$E2, col='red', pch=19)

# Adding scatterplot of gfg_x3 vs gfg_y3
points(dataTable$Y, dataTable$E3, col='blue', pch=19)

# Adding scatterplot of gfg_x4 vs gfg_y4
points(dataTable$Y, dataTable$E4, col='orange', pch=19)

legend('topleft',legend=c('E1', 'E2','E3','E4'),
       pch=c(19, 19), col=c('darkgreen', 'red','blue','orange'))

plot(Y ~ (E1+E2+E3+E4), data=dataTable)

M_E <- lm(Y ~ E1+E2+E3+E4, data=dataTable)
plot(resid(M_E) ~ fitted(M_E), main='Environmental Residual Plot')
plot(M_E)
summary(M_E)
me_rsquared <- summary(M_E)$adj.r.squared


M_raw <- lm(Y ~ (E1+E2+E3+E4+G1+G2+G3+G4+G5+G6+G7+G8+G9+G10+G11+G12+G13+G14+G15+G16+G17+G18+G19+G20)^2, data=dataTable)
plot(resid(M_raw) ~ fitted(M_raw), main='Residual Plot')
bc = boxcox(M_raw, lambda = seq(-1, 5, by = 0.5))
lambda = bc$x[which.max(bc$y)]
summary(M_raw)

mraw_rsquare <- summary(M_raw)$adj.r.square

options(max.print = 9999)

Y_trans = (((dataTable$Y^2.03)-1)/2.03)

M_trans <- lm( Y_trans ~ (E1+E2+E3+E4+G1+G2+G3+G4+G5+G6+G7+G8+G9+G10+G11+G12+G13+G14+G15+G16+G17+G18+G19+G20)^2, data=dataTable)
summary(M_trans)

plot(resid(M_trans) ~ fitted(M_raw), main="Transform Residual Plot")
transformation_rsquare <- summary(M_trans)$adj.r.square


library(leaps)
M <- regsubsets( model.matrix(M_trans)[,-1], Y_trans,
                 nbest = 1 , nvmax=5, 
                 method = 'forward', intercept = TRUE )
temp <- summary(M)

library(knitr)
Var <- colnames(model.matrix(M_trans))
M_select <- apply(temp$which, 1, function(x) paste0(Var[x], collapse='+'))
kable(data.frame(cbind( model = M_select, adjR2 = temp$adjr2, BIC = temp$bic)),
      caption='Model Summary')

##

M_main <- lm(Y_trans ~ (E1+E2+E3+E4+G1+G2+G3+G4+G5+G6+G7+G8+G9+G10+G11+G12+G13+G14+G15+G16+G17+G18+G19+G20), data=dataTable)
temp <- summary(M_main)
summary(M_main)
kable(temp$coefficients[ abs(temp$coefficients[,4]) <= 0.001, ], caption='Sig Coefficients')

M_2stage <- lm(Y_trans ~ (E1+E2+E3)^2, data=dataTable)
temp <- summary(M_2stage)
summary(M_2stage)
kable(temp$coefficients[ abs(temp$coefficients[,3]) >= 1, ])


M_final <- lm(Y_trans ~ E2 + E1:E3 + E2:E3, data = dataTable) #note when including + G4:G20 + G7:G14, pr>t value are small, (<0.01), which may indicate statistical significance, however the effect size of including the variables is negligible.
plot(resid(M_final) ~ fitted(M_final), main="Final Model Residual Plot")
summary(M_final)
plot(M_final)


anova_table <- anova(M_final)
print(anova_table)

plot(Y_trans ~ E2 + E1:E3 + E2:E3, data = dataTable)
#we perform a transformation on each E value ? 