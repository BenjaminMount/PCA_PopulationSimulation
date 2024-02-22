# #1、构造Balding-Nichols模型
# 
# n <- 100000 #对每一个亚种群抽取100000次等位基因频率
# n_sub <- 3 #亚种群群数量
# fk <- c(0.05,0.15,0.25)  
# 
# # 初始化一个空矩阵用于存储所有的等位基因频率
# # 其中第一列对应原祖先种群，每二、三、四列对应三个亚种群
# allele_freqs <- matrix(NA, nrow = n, ncol = n_sub + 1)
# set.seed(1)
# # 循环生成等位基因频率
# for (i in 1:n) 
# {
#   ps <- runif(1, 0.1, 0.9)  
#   allele_freqs[i, 1] <- ps #对locus s处，祖先等位基因频率ps从一个均匀分布U(0.1,0.9)中抽取
#   
#   for (j in 1:n_sub) 
#   {
#     allele_freqs[i, j + 1] <- rbeta(1, ps * (1 - fk[j]) / fk[j], (1 - ps) * (1 - fk[j]) / fk[j])
#     #Balding-Nicholes model 第2、3、4列分别是亚种群2、3、4在locus s处的祖先等位基因频率
#   }
# }
# head(allele_freqs)
# #至此，等位基因频率生成完毕，其中每一个s处的等位基因为一个常量
# 
# #write.csv(allele_freqs, "F:/hust/Doct/Thesis/5.0：模拟/Code/allele_ancestral_freqs.csv", row.names=FALSE)
# #2、接着生成三个不同种群结构的祖先向量
# #install.packages("DirichletReg")
# library(DirichletReg)
# 
# #Population Structure I
# set.seed(2)
# N <- 140 # 假设有140个个体，在简化的家谱中每一个家谱有14个创始人，生成140个个体可以生成包含10个家谱的信息。
# alpha_params <- c(1,1,1)
# # 为每个个体生成一个祖先向量
# ancestry_vectors_I <- rdirichlet(N, alpha_params)
# # 查看前几个祖先向量
# head(ancestry_vectors_I)
# 
# 
# #生成家谱创始人基因型数值
# #种群结构Ⅰ的家谱创始人
# allele_freqs_k <- allele_freqs[,-1] #三个亚种群的等位基因频率矩阵
# genotype_matrix_I <- matrix(NA, nrow = N, ncol = 100000) #基因型数值矩阵
# set.seed(4)
# for (s in 1:100000) 
# {
#   allele_freq_vector <- allele_freqs_k[s,]  #s处的三个亚种群的等位基因频率向量
#   for (i in 1:N) 
#   {
#     ancestry_vector <- ancestry_vectors_I[i,] #第i个家谱创始人的祖先向量
#     binom_param <- sum(allele_freq_vector * ancestry_vector) #二项分布参数
#     genotype_matrix_I[i, s] <- rbinom(1, 2, binom_param) #从二项分布中抽取基因型数值并且存入基因型数值矩阵中
#   }
# }

# write.csv(genotype_matrix_I, "F:/hust/Doct/Thesis/5.0：模拟/Code/Data_Pc_relate_founders_genotype/genotype_matrix_I.csv", row.names=FALSE)
# write.csv(ancestry_vectors_I, "F:/hust/Doct/Thesis/5.0：模拟/Code/Data_Pc_relate_founders_genotype/genotype_matrix_I.csv", row.names=FALSE)

#以上部分在复现PC-Relate时均已完成，在PC-Relate部分中我们按照三个不同的种群结构设置
#生成了三个不同种群结构设置下的祖先等位基因频率、祖先向量以及基于这两者的创始人的等位基因频率矩阵。
#在此部分我们首先研究种群结构Ⅰ下的表现，导入数据即可
#不同的是为了达到等权重假设，我们的谱系不同
library(data.table)

genotype_matrix_I <- fread("F:/hust/Doct/Thesis/5.0：模拟/Code/Data_Pc_relate_founders_genotype/genotype_matrix_I.csv", header = TRUE, sep = ",")
ancestry_vectors_I <- fread("F:/hust/Doct/Thesis/5.0：模拟/Code/Data_Pc_relate_founders_genotype/ancestry_vectors_I.csv", header = TRUE, sep = ",")
genotype_matrix_I <- as.data.frame(genotype_matrix_I)
ancestry_vectors_I <- as.data.frame(ancestry_vectors_I)

#生成我们预期的家族谱系

# 定义一个函数用于根据孟德尔遗传定律从父母基因型数据生成后代基因型数据，不考虑突变
simulate_child_genotype <- function(mother, father) {
  child <- integer(length(mother)) #分配内存
  same_homozygous <- mother == father & (mother %in% c(0, 2)) #找到都是相同纯合子的位置
  child[same_homozygous] <- mother[same_homozygous] #case1 父母都是纯合子
  
  both_heterozygous <- mother == 1 & father == 1 # 找到都是杂合子的位置
  child[both_heterozygous] <- sample(c(0, 1, 1, 2), sum(both_heterozygous), replace = TRUE) #case2 父母都是杂合子
  
  one_AA_one_aa <- (mother == 2 & father == 0) | (mother == 0 & father == 2) #找到一个是AA一个是aa的位置
  child[one_AA_one_aa] <- 1 #case3 父母一方是AA 一方是aa
  
  one_Aa_one_AA <- (mother == 1 & father == 2) | (mother == 2 & father == 1)#找到父母位置
  child[one_Aa_one_AA] <- sample(c(1, 2), sum(one_Aa_one_AA), replace = TRUE) #case4 父母一方是AA一方是Aa
  
  one_Aa_one_aa <- (mother == 1 & father == 0) | (mother == 0 & father == 1)
  child[one_Aa_one_aa] <- sample(c(0, 1), sum(one_Aa_one_aa), replace = TRUE) #case5 父母一方是aa一方是Aa
  
  return(child)
}


#种群结构Ⅰ中操作
set.seed(7)
remaining_genotypes_I = as.matrix(genotype_matrix_I) #复制种群矩阵用于模拟
remaining_ancestry_I = as.matrix(ancestry_vectors_I) # 复制祖先向量用于模拟

# 模拟生成一个Figure_Simplified的家谱生成过程
# 初始家谱基因型矩阵以及祖先向量矩阵用于存储
pedigree <- matrix(NA, nrow = 28, ncol = ncol(remaining_genotypes_I)) 
pedigree_ancestry <- matrix(NA, nrow = 28, ncol = ncol(remaining_ancestry_I)) #每次循环的祖先向量导入值

# 模拟第一代(抽取14个作为创始人 indx:1-14 n=14)
selected_founders_idx <- sample(1:nrow(remaining_genotypes_I), 14)
pedigree[1:14, ] <- remaining_genotypes_I[selected_founders_idx, ] #两个创始人的基因型数据信息
pedigree_ancestry[1:14, ] <- remaining_ancestry_I[selected_founders_idx, ]
remaining_genotypes_I <- remaining_genotypes_I[-selected_founders_idx, , drop = FALSE] #不放回抽样，更新种群，即为剩下的种群
remaining_ancestry_I <- remaining_ancestry_I[-selected_founders_idx, , drop = FALSE] #更新祖先向量
# 模拟第二代 (交配得到8个后代  indx:15-22 n=8)
pedigree[15, ] <- simulate_child_genotype(pedigree[1, ], pedigree[2, ])
pedigree_ancestry[15, ] <- (pedigree_ancestry[1, ] + pedigree_ancestry[2, ]) / 2
pedigree[16, ] <- simulate_child_genotype(pedigree[3, ], pedigree[4, ])
pedigree_ancestry[16, ] <- (pedigree_ancestry[3, ] + pedigree_ancestry[4, ]) / 2
pedigree[17, ] <- simulate_child_genotype(pedigree[5, ], pedigree[6, ])
pedigree_ancestry[17, ] <- (pedigree_ancestry[5, ] + pedigree_ancestry[6, ]) / 2
pedigree[18, ] <- simulate_child_genotype(pedigree[7, ], pedigree[8, ])
pedigree_ancestry[18, ] <- (pedigree_ancestry[7, ] + pedigree_ancestry[8, ]) / 2
pedigree[19, ] <- simulate_child_genotype(pedigree[7, ], pedigree[8, ])
pedigree_ancestry[19, ] <- (pedigree_ancestry[7, ] + pedigree_ancestry[8, ]) / 2
pedigree[20, ] <- simulate_child_genotype(pedigree[9, ], pedigree[10, ])
pedigree_ancestry[20, ] <- (pedigree_ancestry[9, ] + pedigree_ancestry[10, ]) / 2
pedigree[21, ] <- simulate_child_genotype(pedigree[11, ], pedigree[12, ])
pedigree_ancestry[21, ] <- (pedigree_ancestry[11, ] + pedigree_ancestry[12, ]) / 2
pedigree[22, ] <- simulate_child_genotype(pedigree[13, ], pedigree[14, ])
pedigree_ancestry[22, ] <- (pedigree_ancestry[13, ] + pedigree_ancestry[14, ]) / 2
# 模拟第三代 （交配得到4个后代 indx:23-26 n=4）
pedigree[23, ] <- simulate_child_genotype(pedigree[15, ], pedigree[16, ])
pedigree_ancestry[23, ] <- (pedigree_ancestry[15, ] + pedigree_ancestry[16, ]) / 2
pedigree[24, ] <- simulate_child_genotype(pedigree[17, ], pedigree[18, ])
pedigree_ancestry[24, ] <- (pedigree_ancestry[17, ] + pedigree_ancestry[18, ]) / 2
pedigree[25, ] <- simulate_child_genotype(pedigree[19, ], pedigree[20, ])
pedigree_ancestry[25, ] <- (pedigree_ancestry[19, ] + pedigree_ancestry[20, ]) / 2
pedigree[26, ] <- simulate_child_genotype(pedigree[21, ], pedigree[22, ])
pedigree_ancestry[26, ] <- (pedigree_ancestry[21, ] + pedigree_ancestry[22, ]) / 2
# 模拟第四代  （交配得到2个后代 indx:27-28 n=2）
pedigree[27, ] <- simulate_child_genotype(pedigree[23, ], pedigree[24, ])
pedigree_ancestry[27, ] <- (pedigree_ancestry[23, ] + pedigree_ancestry[24, ]) / 2
pedigree[28, ] <- simulate_child_genotype(pedigree[25, ], pedigree[26, ])
pedigree_ancestry[28, ] <- (pedigree_ancestry[25, ] + pedigree_ancestry[26, ]) / 2

#创建一个Half-siblings的情形
selected_founders_idx_half_siblings <- sample(1:nrow(remaining_genotypes_I), 1)
individual_7_new_spouse_genotype <- remaining_genotypes_I[selected_founders_idx_half_siblings, ] # 或者是实际的基因型数据
remaining_genotypes_I <- remaining_genotypes_I[-selected_founders_idx_half_siblings, , drop = FALSE] #不放回抽样，更新种群，即为剩下的种群
# 向 pedigree 矩阵添加第29行用于存储Half_siblings情形的个体7的另外一个配偶
pedigree <- rbind(pedigree, individual_7_new_spouse_genotype)
individual_7_new_spouse_ancestry <- remaining_ancestry_I[selected_founders_idx_half_siblings, ]#对应的祖先向量
pedigree_ancestry <- rbind(pedigree_ancestry, individual_7_new_spouse_ancestry)
remaining_ancestry_I <- remaining_ancestry_I[-selected_founders_idx_half_siblings, , drop = FALSE]

#模拟出一个个体作为个体7与个体29的后代，这意味着个体7进行了二婚，与个体29产生了后代30，个体30与个体18、19为Half_sibilings的关系
half_siblings_of_individual_18_19_genotype <- simulate_child_genotype(pedigree[7, ], pedigree[29, ])
pedigree <- rbind(pedigree, half_siblings_of_individual_18_19_genotype)
half_siblings_of_individual_18_19_ancestry <- (pedigree_ancestry[1, ] + pedigree_ancestry[2, ]) / 2
pedigree_ancestry <- rbind(pedigree_ancestry, half_siblings_of_individual_18_19_ancestry)


# write.csv(pedigree, "F:/hust/Doct/Thesis/5.0：模拟/Code/code_simplified/pedigree_simplified.csv", row.names=FALSE)
# write.csv(pedigree_ancestry, "F:/hust/Doct/Thesis/5.0：模拟/Code/code_simplified/pedigree_simplified_ancestry_vector_I.csv", row.names=FALSE)
#导出后不必每次再随机生成，确保数据一致性

# 函数计算N*3矩阵，每行对应个体的N_Aa_Aa, N_AA_aa, N_Aa_i
calculate_king_robust_matrix <- function(genotype_matrix) {
  # 计算每个个体的N_Aa (即基因型值为1的位点数)
  N_Aa <- rowSums(genotype_matrix == 1)
  
  # 初始化N*N矩阵来存储所有个体之间的KING-ROBUST估计量
  kinship_matrix <- matrix(NA, nrow = nrow(genotype_matrix), ncol = nrow(genotype_matrix))
  
  # 对于矩阵中的每对个体，计算KING-ROBUST估计量
  for (i in 1:(nrow(genotype_matrix) - 1)) {
    for (j in (i + 1):nrow(genotype_matrix)) {
      N_Aa_Aa <- sum(genotype_matrix[i, ] == 1 & genotype_matrix[j, ] == 1)
      N_AA_aa <- sum((genotype_matrix[i, ] == 2 & genotype_matrix[j, ] == 0) | (genotype_matrix[i, ] == 0 & genotype_matrix[j, ] == 2))
      kinship_matrix[i, j] <- (N_Aa_Aa - 2 * N_AA_aa) / (N_Aa[i] + N_Aa[j])
      kinship_matrix[j, i] <- kinship_matrix[i, j] # kinship is symmetric
    }
  }
  
  # 对角线元素，个体与自己的亲缘系数设为0.5
  diag(kinship_matrix) <- 0.5
  
  return(kinship_matrix)
}

# 计算KING-ROBUST估计量矩阵
king_robust_matrix <- calculate_king_robust_matrix(pedigree)
print(king_robust_matrix)

# 函数识别与给定个体没有亲缘关系的个体集合
find_non_related_individuals <- function(kinship_matrix, target_individual_index, threshold = c(-0.025, 0.025)) {
  # 识别与目标个体没有亲缘关系的个体
  # 这里我们检查亲缘系数矩阵中目标个体对应行的所有元素，以判断每个个体是否与目标个体相关
  non_related_indices <- which(kinship_matrix[target_individual_index, ] >= threshold[1] & 
                                 kinship_matrix[target_individual_index, ] <= threshold[2])
  
  # 移除目标个体自身的索引
  non_related_indices <- non_related_indices[non_related_indices != target_individual_index]
  
  return(non_related_indices)
}

#明确我们接下来的步骤，我们分别利用我们的思路来估计个体18、19之间的亲缘系数，个体18、19是来自个体7、8的Full_siblings，其
#亲缘系数在理论上应当是0.25，然后我们要研究个体18与个体7之间的Parent/offsprings的亲缘系数，其亲缘系数在理论上应当是0.25，
#接着我们要研究个体18与个体30之间的亲缘系数，他们之间在理论上是Half_siblings的关系，其亲缘系数理论值为0.125.
#首先我们要利用PC-AiR的思想来估计参考等位基因频率的大小。

# 找出与个体18和19没有亲缘关系的个体集合
non_related_to_18 <- find_non_related_individuals(king_robust_matrix, 18)
non_related_to_19 <- find_non_related_individuals(king_robust_matrix, 19)

pc_air_founders <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,29)

#对于互不相关的个体，我们首先封装一个函数使得能够基于一个互不相关的集合来估计个体的等位基因频率

calculate_hatD_M <- function(genotype_matrix) {
  delta_M_i <- rowSums(genotype_matrix * (2 - genotype_matrix)) / ncol(genotype_matrix)
  hatD_M <- diag(delta_M_i)
  return(hatD_M)
}
# 示例：使用基因型矩阵 Y 调用这个函数
# hatD_M <- calculate_hatD_M(Y)
# print(hatD_M)

calculate_hatQ <- function(genotype_matrix, hatD_M, K) {
  M <- ncol(genotype_matrix) 
  # 特征值分解
  eigen_decomp <- eigen((1/M) * t(genotype_matrix) %*% genotype_matrix - hatD_M)
  # 选择与最大的K个正特征值相对应的特征向量
  Q_hat <- eigen_decomp$vectors[, 1:K]
  return(Q_hat)
}

#如何确定K？
# 示例：使用基因型矩阵Y和对角矩阵hatD_M调用这个函数
# Q_hat <- calculate_hatQ(Y, hatD_M, K)
# print(Q_hat)
calculate_phi_hat <-function(genotype_matrix, K){
  hatD_M <-calculate_hatD_M(genotype_matrix)
  Q_hat <- calculate_hatQ(genotype_matrix, hatD_M, K)
  Pi_hat <- (1/2) * genotype_matrix %*% t(Q_hat) %*% solve(Q_hat %*% t(Q_hat)) %*% Q_hat
  return(Pi_hat)
}

#对于个体18我们需要得到与其相关，但又互不相关的个体集合，算法设计暂且不讨论，我们在模拟部分已知个体7和8最符合条件。
founders_18 <- seq(1, 14)
founders_19 <- c(7, 8) #个体19也是同理的
#接着我们计算个体18和19在各个位点上的等位基因频率
Pi_hat_18 <- calculate_phi_hat(t(pedigree[founders_18, ]) ,2)


