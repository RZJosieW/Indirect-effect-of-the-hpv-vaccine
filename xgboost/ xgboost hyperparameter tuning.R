library(readxl)
library(xgboost)
library(Metrics)
library(dplyr)

set.seed(123)

hpvdata <- read_excel("/Users/V/Desktop/hpvfinal.xlsx")
colnames(hpvdata)[colnames(hpvdata) == "Pre HPV-16, -18"] <- "pre"
colnames(hpvdata)[colnames(hpvdata) == "Post HPV-16, -18"] <- "post"
colnames(hpvdata)[colnames(hpvdata) == "vaccine coverage"] <- "coverage"
colnames(hpvdata)[colnames(hpvdata) == "age"] <- "age1"
colnames(hpvdata)[colnames(hpvdata) == "agemid"] <- "age"

hpv_filtered <- hpvdata %>% select(indirect, coverage, age, year)

X <- as.matrix(hpv_filtered[, c("coverage", "age", "year")])
y <- hpv_filtered$indirect
n <- nrow(X)

learning_rates <- c(0.01, 0.05, 0.1)
max_depths <- c(2, 3, 4)
subsamples <- c(0.7, 0.9)
reg_alphas <- c(0, 0.1)

results <- data.frame()
best_mae <- Inf
best_params <- list()


for (lr in learning_rates) {
  for (depth in max_depths) {
    for (sub in subsamples) {
      for (alpha in reg_alphas) {
        
        preds <- numeric(n)
        
        for (i in 1:n) {
          X_train <- X[-i, ]
          y_train <- y[-i]
          X_test <- X[i, , drop = FALSE]
          
          dtrain <- xgb.DMatrix(data = X_train, label = y_train)
          dtest <- xgb.DMatrix(data = X_test)
          
          model <- xgb.train(
            params = list(
              objective = "reg:squarederror",
              learning_rate = lr,
              max_depth = depth,
              subsample = sub,
              reg_alpha = alpha,
              verbosity = 0
            ),
            data = dtrain,
            nrounds = 100,
            verbose = 0
          )
          
          preds[i] <- predict(model, dtest)
        }
        
        mae_score <- mae(y, preds)
        r2_score <- 1 - sum((y - preds)^2) / sum((y - mean(y))^2)
        
        cat(sprintf("lr = %.2f | depth = %d | subsample = %.1f | alpha = %.2f --> MAE = %.4f | R² = %.4f\n",
                    lr, depth, sub, alpha, mae_score, r2_score))
        
        results <- rbind(results, data.frame(
          learning_rate = lr,
          max_depth = depth,
          subsample = sub,
          reg_alpha = alpha,
          MAE = mae_score,
          R2 = r2_score
        ))
        
        if (mae_score < best_mae) {
          best_mae <- mae_score
          best_params <- list(
            learning_rate = lr,
            max_depth = depth,
            subsample = sub,
            reg_alpha = alpha
          )
        }
      }
    }
  }
}

print(best_params)




dtrain_all <- xgb.DMatrix(data = X, label = y)

final_model <- xgb.train(
  params = list(
    objective = "reg:squarederror",
    learning_rate = 0.05,
    max_depth = 2,
    subsample = 0.7,
    reg_alpha = 0.10,
    verbosity = 0
  ),
  data = dtrain_all,
  nrounds = 100
)


library(SHAPforxgboost)

shap_values <- shap.values(xgb_model = final_model, X_train = X)
shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = X)

shap.plot.summary(shap_long)
shap.plot.dependence(data_long = shap_long, x = "age")
shap.plot.dependence(data_long = shap_long, x = "coverage")
shap.plot.dependence(data_long = shap_long, x = "year")
