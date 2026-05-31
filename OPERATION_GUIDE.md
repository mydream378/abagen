# PD vs HC 胆碱能受体分析 - 操作指南

## 📋 分析概述

本分析用于比较**帕金森病(PD)患者**和**健康对照组(HC)**在T1加权MRI影像中胆碱能受体相关基因的皮层表达差异。

## 🎯 分析内容

### 分析的胆碱能相关基因（10个）：

| 基因 | 全名 | 功能 |
|------|------|------|
| **CHAT** | 胆碱乙酰转移酶 | 合成乙酰胆碱 |
| **SLC5A7** | 胆碱转运体 | 胆碱再摄取 |
| **CHRNA4** | nAChR α4 | 烟碱受体亚基 |
| **CHRNA5** | nAChR α5 | 烟碱受体亚基 |
| **CHRNA7** | nAChR α7 | 烟碱受体亚基 |
| **CHRNB2** | nAChR β2 | 烟碱受体亚基 |
| **CHRM1** | mAChR M1 | 毒蕈碱受体 |
| **CHRM2** | mAChR M2 | 毒蕈碱受体 |
| **CHRM3** | mAChR M3 | 毒蕈碱受体 |
| **CHRM4** | mAChR M4 | 毒蕈碱受体 |

## 🚀 快速开始

### 方法一：使用模拟数据演示

```bash
# 1. 创建示例数据结构
python pd_hc_cholinergic_simple.py --create-template ./example_data

# 2. 运行分析
python pd_hc_cholinergic_simple.py --data-dir ./example_data
```

### 方法二：使用真实数据

```bash
# 1. 准备数据目录
python pd_hc_cholinergic_simple.py --create-template ./my_study

# 2. 放入你的数据后运行
python pd_hc_cholinergic_simple.py --data-dir ./my_study
```

## 📁 数据准备

### 1. 受试者信息表（必需）

创建 `subjects.csv` 文件：

```csv
subject_id,group,age,sex
PD001,PD,65,M
PD002,PD,68,F
PD003,PD,72,M
PD004,PD,61,F
PD005,PD,70,M
PD006,PD,66,F
PD007,PD,69,M
PD008,PD,64,F
PD009,PD,71,M
PD010,PD,67,F
HC001,HC,63,M
HC002,HC,66,F
HC003,HC,70,M
HC004,HC,59,F
HC005,HC,72,M
HC006,HC,68,F
HC007,HC,65,M
HC008,HC,69,F
HC009,HC,67,M
HC010,HC,71,F
HC011,HC,64,M
HC012,HC,68,F
HC013,HC,70,M
HC014,HC,62,F
HC015,HC,66,M
HC016,HC,69,F
HC017,HC,73,M
HC018,HC,67,F
HC019,HC,65,M
HC020,HC,71,F
```

**注意**：
- `subject_id`: 受试者唯一标识
- `group`: 必须是 "PD" 或 "HC"
- `age`: 年龄（用于协变量分析）
- `sex`: 性别（M/F，用于协变量分析）

### 2. 基因表达数据（可选）

如果已经有基因表达数据，创建 `cholinergic_expression.csv`：

```csv
subject_id,lh_bankssts,lh_caudal_acfb,lh_caudal_mfg,...
PD001,1.234,0.892,1.456,...
PD002,1.456,0.923,1.678,...
...
```

**格式要求**：
- 第一列为受试者ID（与subjects.csv一致）
- 其他列为脑区名称（72个Desikan-Killiany脑区）
- 数值为基因表达值

### 3. T1图像和皮层厚度（可选）

```
data_dir/
├── subjects.csv              # 受试者信息（必需）
├── cholinergic_expression.csv # 基因表达数据（可选，已有则跳过）
├── t1_images/                # T1图像（可选）
│   ├── PD001_T1.nii.gz
│   └── ...
└── cortical_thickness/       # 皮层厚度数据（可选）
    ├── PD001_thickness.csv
    └── ...
```

## 📊 分析流程

```
1. 数据加载
   └─> subjects.csv
   └─> cholinergic_expression.csv (如无则生成模拟数据)

2. 统计分析
   └─> 独立样本t检验 (PD vs HC)
   └─> FDR多重比较校正
   └─> Cohen's d 效应量计算

3. 结果输出
   ├─> group_comparison_results.csv  (详细统计结果)
   ├─> gene_expression_patterns.csv  (基因表达模式)
   └─> summary_report.txt            (可读报告)
```

## 📈 结果解读

### CSV文件字段说明

**group_comparison_results.csv**：

| 字段 | 说明 |
|------|------|
| `region` | 脑区名称 |
| `n_pd`, `n_hc` | 各组样本量 |
| `mean_pd`, `mean_hc` | 各组平均值 |
| `std_pd`, `std_hc` | 各组标准差 |
| `mean_diff` | 组间差异 |
| `t_statistic` | t检验统计量 |
| `p_value` | 原始p值 |
| `cohens_d` | 效应量 |
| `effect_size_interpretation` | 效应量解释 |
| `p_fdr_corrected` | FDR校正后p值 |
| `significant` | 是否显著（TRUE/FALSE） |

### 效应量解释 (Cohen's d)

- **|d| < 0.2**: 微小效应（实际意义很小）
- **0.2 ≤ |d| < 0.5**: 小效应
- **0.5 ≤ |d| < 0.8**: 中等效应（有一定临床意义）
- **|d| ≥ 0.8**: 大效应（有重要临床意义）

### FDR校正

FDR (False Discovery Rate) 校正用于控制多重比较导致的假阳性：
- **p_fdr < 0.05**: 认为组间存在显著差异
- **p_fdr ≥ 0.05**: 组间差异不显著

## 🎨 扩展分析建议

### 1. 添加协变量（年龄、性别）

```python
# 在 perform_group_comparison() 函数中添加
from scipy import stats
import statsmodels.api as sm

def compare_with_covariates(pd_values, hc_values, pd_ages, hc_ages):
    """控制年龄后的组间比较"""
    # 合并数据
    all_values = np.concatenate([pd_values, hc_values])
    all_ages = np.concatenate([pd_ages, hc_ages])
    groups = np.concatenate([np.ones(len(pd_values)), np.zeros(len(hc_values))])
    
    # 回归分析
    X = sm.add_constant(np.column_stack([groups, all_ages]))
    model = sm.OLS(all_values, X).fit()
    
    return model
```

### 2. 非参数检验

如果数据不符合正态分布，使用Mann-Whitney U检验：

```python
from scipy.stats import mannwhitneyu

stat, p_value = mannwhitneyu(pd_values, hc_values, alternative='two-sided')
```

### 3. 置换检验

对于小样本数据：

```python
from sklearn.utils import shuffle

def permutation_test(pd_values, hc_values, n_permutations=10000):
    observed_diff = np.mean(pd_values) - np.mean(hc_values)
    combined = np.concatenate([pd_values, hc_values])
    n = len(pd_values)
    
    count = 0
    for _ in range(n_permutations):
        shuffled = shuffle(combined)
        perm_diff = np.mean(shuffled[:n]) - np.mean(shuffled[n:])
        
        if abs(perm_diff) >= abs(observed_diff):
            count += 1
    
    return count / n_permutations
```

### 4. 脑区可视化

使用 nilearn 进行结果可视化：

```python
import nilearn.plotting as plotting
import numpy as np

# 创建一个示例统计图
plotting.plot_stat_map(
    stat_img,
    title='PD vs HC Cholinergic Expression',
    cut_coords=7,
    display_mode='mosaic'
)
```

## ❓ 常见问题

### Q1: 需要多少样本？
**A**: 每组至少3个样本才能进行t检验。建议每组至少10-20个样本以获得可靠的统计结果。

### Q2: 数据缺失怎么办？
**A**: 脚本会自动跳过缺失值。也可以使用线性插值填充：
```python
expression_df.interpolate(method='linear', inplace=True)
```

### Q3: 如何使用其他脑区模板？
**A**: 可以修改代码中的 `DESIKAN_KILLIANY_REGIONS` 列表，使用其他模板（如Schaefer、AAL等）。

### Q4: 如何处理纵向数据？
**A**: 需要使用重复测量ANOVA或混合效应模型：
```python
import statsmodels.formula.api as smf

model = smf.mixedlm("expression ~ group * time", 
                     data=df, 
                     groups=df["subject_id"]).fit()
```

### Q5: 可以分析其他基因吗？
**A**: 可以！修改代码中的 `CHOLINERGIC_GENES` 列表：
```python
MY_GENES = [
    'TH',        # 酪氨酸羟化酶
    'DAT',       # 多巴胺转运体
    'SNCA',      # α-突触核蛋白
    # 添加你想分析的基因
]
```

## 📚 相关资源

1. **Allen Human Brain Atlas**: https://human.brain-map.org/
2. **abagen工具箱**: https://github.com/rmarkello/abagen
3. **FreeSurfer**: https://surfer.nmr.mgh.harvard.edu/
4. **FSL**: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki

## 🔬 参考文献

1. Hawrylycz et al. (2012). An anatomically comprehensive atlas of the adult human brain transcriptome. Nature.
2. Arnatkevičiūtė et al. (2019). The importance of data processing at the whole-brain level. NeuroImage.
3. Markello et al. (2021). Misalignment of AHBA data with MNI space. bioRxiv.

## 📞 技术支持

如有问题，请检查：
1. 数据文件格式是否正确
2. 依赖包是否正确安装
3. 路径是否正确

运行测试：
```bash
python pd_hc_cholinergic_simple.py --help
```
