data = t(data_denoised_1[, 4:103])  # 矩阵转置

st <- function(data){
  dimension = dim(data)  # 矩阵维数
  Si = 1/dimension[1] * colSums(data)  # 原矩阵每列均值
  data_rand = matrix(0, nrow = dimension[1], ncol = dimension[2])  # 随机排序矩阵初始化
  loop1 = 50  # 多次实验便于画出平滑的ROC
  loop2 = 1000  # 循环次数
  Smax = rep(0, loop2)  # 每个矩阵的列均值最大初始化
  p_value = rep(0, dimension[2])  # 每个特征的判决值初始化
  a = 0.05  # 阈值
  # 突变范围
  ap_s1 = 100
  ap_e1 = 149
  ap_s2 = 500
  ap_e2 = 529
  ap_s3 = 900
  ap_e3 = 919
  change_range = (ap_e1 - ap_s1 + 1) + (ap_e2 - ap_s2 + 1) + (ap_e2 - ap_s2 + 1)  # 突变个数
  x_y = array(0, dim=c(2, dimension[2], loop1))  # 多组ROc坐标点初始化
  x_y_result = matrix(0, nrow=2, ncol=dimension[2])  # 平滑ROC坐标点初始化
  step = 1  # 步长
  for(k in c(1:loop1)){
    #
    # 检测突变数据
    #
    cat("start ", k, "\n")
    for(i in c(1:loop2)){  # 构建1000个按行随机排列的矩阵
      for(j in c(1:dimension[1])){  # 行随机排列矩阵产生
        rand = sample(1:dimension[2], dimension[2])  # 矩阵行随机排序顺序产生
        data_rand[j,] = data[j, rand]  # 按产生的随机顺序赋值矩阵每一行
        }
      Si_rand = 1/dimension[1] * colSums(data_rand)  # 行随机排列矩阵每列均值
      Smax[i] = max(Si_rand)  # 找出列均值最大
    }
    for(i in c(1:dimension[2])){  # 计算每个特征的判决值
      I_sum = 0
      for(j in c(1:loop2)){
        if(Si[i]<=Smax[j]){
          I_sum = I_sum + 1
        }
      }
      p_value[i] = I_sum / 1000
    }
    # result = p_value <= a  # 检测到的突变数据 T or F
    #
    # ROC
    #
    p_value_argsort = order(p_value)  # p-value从小到大排序的索引
    # result_argsort = result[p_value_argsort]  # 检测结果排序
    x_start = 0
    y_start = 0
    for(i in c(1:dimension[2])){
      if((p_value_argsort[i]>=ap_s1 && p_value_argsort[i]<=ap_e1) || (p_value_argsort[i]>=ap_s2 && p_value_argsort[i]<=ap_e2) || (p_value_argsort[i]>=ap_s3 && p_value_argsort[i]<=ap_e3)){
        y_start = y_start + 1
        }
      else{
        x_start = x_start + 1
        }
      x_y[1, i, k] = x_start / (dimension[2] - change_range)
      x_y[2, i, k] = y_start / change_range
      }
    }
  #
  # 绘制平均ROC
  #
  cat("ROC\n")
  for(i in c(1:dimension[2])){
    dim3_x_sum = 0
    dim3_y_sum = 0
    for(j in c(1:loop1)){
      dim3_x_sum = dim3_x_sum + x_y[1, i, j]
      dim3_y_sum = dim3_y_sum + x_y[2, i, j]
    }
    x_y_result[1, i] = dim3_x_sum / loop1
    x_y_result[2, i] = dim3_y_sum / loop1
  }
  plot(x_y_result[1,], x_y_result[2,], type = 'l')
}