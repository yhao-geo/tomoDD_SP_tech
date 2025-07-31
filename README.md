# tomoDD-SP：双波层析成像软件包

tomoDD-SP 是一个基于双差层析成像方法的地震定位与速度结构反演软件包，支持 P 波和 S 波联合反演。本软件包专为教学和研究目的设计，提供完整的数据处理流程和示例。

## 功能特点

- 支持 P 波和 S 波走时联合反演
- 提供完整的数据预处理工具链
- 包含多个实际应用案例
- 支持并行计算
- 提供详细的示例数据和教程

## 目录结构

```
tomoDD_SP_tech/
├── src/                    # 源代码
│   ├── src_ph2dt/         # 震相数据预处理模块
│   └── src_tomoDD_SP_parallel/  # 主程序源代码（并行版本）
├── cases/                  # 示例案例
│   ├── 01_case_data_prepare/    # 数据准备示例
│   ├── 02_case_field_demo_chuandian/  # 川滇地区实例
│   ├── 03_case_synthetic_checkerboard/ # 棋盘格测试
│   └── 04_case_L-Curve_selection/      # 参数选择示例
├── tools/                  # 辅助工具
│   ├── CoordTransformTools/     # 坐标转换工具
│   └── ModelingTools/           # 模型处理工具
└── env/                    # 环境配置文件
```

## 安装说明

### 环境要求
- Linux/Unix 操作系统
- Fortran 编译器（gfortran 推荐）
- Python 3.7+

### 安装步骤

1. 克隆仓库
```bash
git clone https://github.com/yhao-geo/tomoDD_SP_tech.git
cd tomoDD_SP_tech
```

2. 编译源代码
```bash
cd src/src_ph2dt
make clean; make
cd ../src_tomoDD_SP_parallel
make clean; make
```

3. 配置 Python 环境
```bash
conda env create -f env/environment.yml
conda activate tomoddsp-tech
```

## 使用流程

### 1. 数据预处理
- 震相数据选择与格式化
- 使用 ph2dt 构建事件对
- 提取 S-P 走时数据

### 2. 参数配置
- 编辑 `ph2dt.inp` 设置预处理参数
- 编辑 `tomoDD_SP.inp` 设置反演参数
- 准备初始速度模型

### 3. 运行程序
```bash
cd cases/your_case_directory
../src/src_tomoDD_SP_parallel/tomoDD_SP
```

### 4. 结果分析
- 使用 `Script` 目录下的绘图脚本可视化结果
- 检查残差文件评估反演质量
- 分析速度模型输出

## 示例数据

本软件包提供了完整的示例数据集：

1. **数据准备示例**
   - 展示完整的数据预处理流程
   - 包含参数设置说明

2. **川滇地区实例**
   - 实际应用案例
   - 包含完整的输入输出数据
   - 提供结果可视化脚本

3. **棋盘格测试**
   - 用于评估反演分辨率
   - 包含合成数据生成脚本

4. **L曲线参数选择**
   - 展示关键参数选择方法
   - 包含参数测试脚本

## 参考文献

如果您在研究中使用了本软件，请引用以下文献：

Guo H, Zhang H, Froment B. Structural control on earthquake behaviors revealed by high-resolution Vp/Vs imaging along the Gofar transform fault, East Pacific Rise[J]. Earth and Planetary Science Letters, 2018, 499: 243-255.

Zhang H, Thurber C, Bedrosian P. Joint inversion for vp, vs, and vp/vs at SAFOD, Parkfield, California[J]. Geochemistry, Geophysics, Geosystems, 2009, 10(11).

Zhang H, Thurber C H. Double-difference tomography: The method and its application to the Hayward fault, California[J]. Bulletin of the Seismological Society of America, 2003, 93(5): 1875-1889.


## 联系方式

youngh_geo@mail.ustc.edu.cn
