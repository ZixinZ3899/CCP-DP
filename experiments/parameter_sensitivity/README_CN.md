# CCP-DP 参数敏感性实验

本目录的正式复现说明请参见英文文档：

- [README.md](README.md)

该实验采用单因素变化设计：基准配置使用正式 CCP-DP 默认参数，每个
非基准配置只改变一个参数族，其余算法设置保持不变。

常用命令：

    mkdir -p build

    g++ \
      -std=c++17 \
      -O3 \
      -fopenmp \
      experiments/parameter_sensitivity/main_sensitivity.cpp \
      -o build/ccpdp_sensitivity

    JOBS=2 \
    OUTPUT_DIR=/tmp/ccpdp_parameter_sensitivity_smoke \
    experiments/parameter_sensitivity/run_sensitivity.sh smoke

完整实验、输入输出、参数范围和论文图表重建方式请以 `README.md` 为准。
