# 加载必要的R包
library(readr)
library(dplyr)
library(ggpubr)
library(stats)
library(nortest)  # 确保加载了 nortest 包以使用 ad.test
library(ggsci)

# 读取数据集
data <- read_csv("C:/Users/拯救者233/Desktop/TCGA_InfoWithGrade.csv")

# 区分连续变量和分类变量
continuous_var <- "Age_at_diagnosis"
categorical_vars <- setdiff(names(data), c(continuous_var, "Grade"))

# 保存结果的列表
results_list <- list()

# 设置显著水平
alpha <- 0.05

# 针对连续变量：Age_at_diagnosis 进行正态性检验
# 先按Grade分组，进行Anderson-Darling检验
group_0 <- data$Age_at_diagnosis[data$Grade == 0]
group_1 <- data$Age_at_diagnosis[data$Grade == 1]

# 进行 Anderson-Darling 检验
ad_test_result_0 <- ad.test(group_0)
ad_test_result_1 <- ad.test(group_1)

# 保存 Anderson-Darling 检验结果
results_list[["Age_at_diagnosis_anderson_darling_test"]] <- data.frame(
  Variable = "Age_at_diagnosis",
  Method = "Anderson-Darling test",
  p_value = c(ad_test_result_0$p.value, ad_test_result_1$p.value),
  statistic = c(ad_test_result_0$statistic, ad_test_result_1$statistic),
  Significance = ifelse(c(ad_test_result_0$p.value, ad_test_result_1$p.value) < alpha, "Significant", "Not Significant"),
  Additional_Info = NA
)

# 进行相关性检验
if (all(c(ad_test_result_0$p.value, ad_test_result_1$p.value) >= alpha)) {
  # 先进行方差齐性检验
  var_test_result <- var.test(Age_at_diagnosis ~ Grade, data = data)
  
  # 根据方差齐性检验结果决定t检验方式
  if (var_test_result$p.value > alpha) {
    # 方差齐性成立，执行标准 t 检验
    age_test_result <- t.test(Age_at_diagnosis ~ Grade, data = data, var.equal = TRUE)
    method_used <- "t-test (equal variance)"
  } else {
    # 方差齐性不成立，执行 Welch t 检验
    age_test_result <- t.test(Age_at_diagnosis ~ Grade, data = data, var.equal = FALSE)
    method_used <- "Welch t-test (unequal variance)"
  }
  
  results_list[["Age_at_diagnosis_t_test"]] <- data.frame(
    Variable = "Age_at_diagnosis",
    Method = method_used,
    p_value = age_test_result$p.value,
    statistic = age_test_result$statistic,
    Significance = ifelse(age_test_result$p.value < alpha, "Significant", "Not Significant"),
    Additional_Info = NA
  )
} else {
  # 否则进行Mann-Whitney U检验
  age_test_result <- wilcox.test(Age_at_diagnosis ~ Grade, data = data)
  results_list[["Age_at_diagnosis_mann_whitney"]] <- data.frame(
    Variable = "Age_at_diagnosis",
    Method = "Mann-Whitney U test",
    p_value = age_test_result$p.value,
    statistic = age_test_result$statistic,
    Significance = ifelse(age_test_result$p.value < alpha, "Significant", "Not Significant"),
    Additional_Info = NA
  )
}

# 针对分类变量进行卡方检验或Fisher精确检验
for (var in categorical_vars) {
  # 创建变量与标签的列联表
  tbl <- table(data[[var]], data$Grade)
  
  # 检查期望频数是否低于5
  expected_counts <- chisq.test(tbl, simulate.p.value = TRUE)$expected
  if (any(expected_counts < 5)) {
    # 期望频数低于5，进行Fisher精确检验
    fisher_result <- fisher.test(tbl)
    results_list[[paste("Fisher_test", var, sep = "_")]] <- data.frame(
      Variable = var,
      Method = "Fisher exact test",
      p_value = fisher_result$p.value,
      statistic = NA,
      Significance = ifelse(fisher_result$p.value < alpha, "Significant", "Not Significant"),
      Additional_Info = NA
    )
  } else {
    # 期望频数均不低于5，进行不带校正的卡方检验
    chi_square_result <- chisq.test(tbl, correct = FALSE)
    results_list[[paste("Chi_square_test", var, sep = "_")]] <- data.frame(
      Variable = var,
      Method = "Chi-square test",
      p_value = chi_square_result$p.value,
      statistic = chi_square_result$statistic,
      Significance = ifelse(chi_square_result$p.value < alpha, "Significant", "Not Significant"),
      Additional_Info = NA
    )
  }
}

# 将所有结果合并为一个数据框
final_results <- do.call(rbind, results_list)

# 格式化p值，保留三位小数，<0.001用"<0.001"表示
final_results$p_value <- sapply(final_results$p_value, function(p) {
  if (p < 0.001) {
    return("<0.001")
  } else {
    return(formatC(p, format = "f", digits = 3))
  }
})

# 在R控制台中显示所有结果
print(final_results)

# 将结果保存到CSV文件
write.csv(final_results, "C:/Users/拯救者233/Desktop/TCGA_Analysis_Results.csv", row.names = FALSE)
