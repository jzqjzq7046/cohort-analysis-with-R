setwd("E:/YFS/workspace/Microbiome_CrossSectional")

# 数据路径
data_dir <- "data"
result_dir <- "result"
figure_dir <- "figure"
table_dir <- "table"

# 随机种子
set.seed(2026)

# 协变量
covariates <- c(
  "IKA18",
  "SP",
  "BMI18"
)

# 结局变量
outcome <- "HYPERLIPIDEMIA_18"

# FDR阈值
alpha <- 0.05

# Create directories
dir.create(result_dir, showWarnings = FALSE)
dir.create(figure_dir, showWarnings = FALSE)
dir.create(table_dir, showWarnings = FALSE)

