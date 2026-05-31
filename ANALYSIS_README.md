# PD vs HC T1影像胆碱能受体皮层差异分析

## 概述

本分析脚本用于比较帕金森病(PD)患者和健康对照组(HC)在T1加权MRI影像中胆碱能受体相关基因的皮层表达差异。

## 分析流程

```
1. 数据准备
   ├── T1加权MRI图像
   ├── 受试者信息表
   └── 皮层厚度数据（可选）
   
2. 影像预处理
   ├── 皮层分割（使用FreeSurfer）
   ├── 脑区划分（Desikan-Killiany模板）
   └── 形态学指标提取

3. 基因表达分析
   ├── 从Allen Human Brain Atlas获取胆碱能受体基因表达
   └── 使用abagen包进行标准化处理

4. 统计分析
   ├── 独立样本t检验
   ├── FDR多重比较校正
   └── 效应量计算（Cohen's d）

5. 结果报告
   ├── 组间差异结果表
   ├── 相关分析结果（可选）
   └── 可视化报告
```

## 胆碱能相关基因列表

本分析涵盖以下胆碱能受体基因：

- **CHAT**: 胆碱乙酰转移酶
- **SLC5A7**: 胆碱转运体
- **CHRNA4, CHRNA5, CHRNA7**: 烟碱型受体α亚基
- **CHRNB2**: 烟碱型受体β亚基
- **CHRM1-CHRM4**: 毒蕈碱型受体M1-M4

## 安装依赖

```bash
pip install abagen numpy pandas scipy statsmodels nibabel
```

## 数据准备

### 1. 创建受试者信息表

创建CSV文件（如 `subjects.csv`），包含以下列：

```csv
subject_id,group,age,sex
PD001,PD,65,M
PD002,PD,68,F
HC001,HC,63,M
HC002,HC,66,F
...
```

### 2. 准备T1图像

将T1加权MRI图像按以下结构组织：

```
data_dir/
├── subjects.csv
├── t1_images/
│   ├── PD001_T1.nii.gz
│   ├── PD002_T1.nii.gz
│   └── ...
└── cortical_thickness/
    ├── PD001_thickness.csv
    ├── PD002_thickness.csv
    └── ...
```

### 3. 皮层厚度数据格式

皮层厚度CSV文件应包含脑区名称和厚度值：

```csv
region,thickness
lh_bankssts,2.45
lh_caudal_acfb,2.38
lh_caudal_mfg,2.51
...
```

## 使用方法

### 基本用法

```bash
python pd_hc_cholinergic_analysis.py \
    --data-dir /path/to/your/data \
    --subjects-file subjects.csv \
    --output ./results
```

### 完整参数

```bash
python pd_hc_cholinergic_analysis.py \
    --data-dir /path/to/your/data \
    --subjects-file subjects.csv \
    --output ./results \
    --atlas desikan_killiany
```

### 参数说明

- `--data-dir`: 数据目录路径（必需）
- `--subjects-file`: 受试者信息文件名，默认 `subjects.csv`
- `--output`: 输出目录，默认 `./results`
- `--atlas`: 脑区划分模板，可选 `desikan_killiany` 或 `schaefer`

## 输出文件

分析完成后会在输出目录生成以下文件：

1. **group_comparison.csv**: 组间比较的详细统计结果
   - 包含：脑区名称、样本量、均值、标准差、t统计量、p值、效应量等

2. **correlations.csv**: 基因表达与皮层厚度的相关分析（如果提供皮层数据）

3. **analysis_report.txt**: 可读的分析报告摘要

## 结果解释

### 统计结果字段

- `region`: 脑区名称
- `n_pd`, `n_hc`: 各组样本量
- `mean_pd`, `mean_hc`: 各组平均值
- `mean_diff`: 组间差异
- `t_statistic`: t检验统计量
- `p_value`: 原始p值
- `p_fdr_corrected`: FDR校正后的p值
- `cohen_d`: 效应量
- `significant`: 是否显著（FDR < 0.05）

### 效应量解释 (Cohen's d)

- **|d| < 0.2**: 微小效应
- **|d| ≈ 0.5**: 中等效应
- **|d| ≈ 0.8**: 大效应

## 注意事项

1. **样本量要求**: 每组至少3个样本才能进行有效的t检验
2. **数据质量**: 确保T1图像质量良好，皮层分割准确
3. **协变量**: 如有需要，可以在代码中添加年龄、性别等协变量的回归
4. **多重比较**: 使用FDR校正控制假阳性率

## 扩展分析

### 添加更多协变量

可以在 `perform_group_analysis()` 函数中添加协变量回归：

```python
from statsmodels.api import OLS

# 控制年龄和性别
X = subjects_df[['age', 'sex_numeric']]
X = sm.add_constant(X)

for region in regions:
    y = expression_df[region]
    model = OLS(y, X).fit()
    # 提取组间效应
```

### 置换检验

对于小样本数据，可以使用置换检验：

```python
from sklearn.utils import shuffle

def permutation_test(pd_values, hc_values, n_permutations=10000):
    observed_diff = np.mean(pd_values) - np.mean(hc_values)
    combined = np.concatenate([pd_values, hc_values])
    n_pd = len(pd_values)
    
    count = 0
    for _ in range(n_permutations):
        shuffled = shuffle(combined)
        perm_pd = shuffled[:n_pd]
        perm_hc = shuffled[n_pd:]
        perm_diff = np.mean(perm_pd) - np.mean(perm_hc)
        
        if abs(perm_diff) >= abs(observed_diff):
            count += 1
    
    return count / n_permutations
```

### 脑区可视化

可以使用nilearn或brainrender进行结果可视化：

```python
import nilearn.plotting as plotting

# 显示显著差异脑区
display = plotting.plot_stat_map(
    stat_map_img,
    title='PD vs HC Cholinergic Expression Differences',
    cut_coords=7,
    display_mode='mosaic'
)
```

## 常见问题

**Q: 如何处理缺失数据？**  
A: 脚本会自动跳过缺失值，也可以使用插值方法填充。

**Q: 可以使用其他脑区模板吗？**  
A: 可以修改 `atlas` 参数或使用自定义模板。

**Q: 如何进行纵向分析？**  
A: 需要添加时间点变量，并使用重复测量 ANOVA 或混合效应模型。

## 参考文献

1. Arnatkevičiūtė et al. (2019). The importance of data processing at the
   whole-brain level for microarray expression data.
   NeuroImage, 184, 291-308.

2. Hawrylycz et al. (2012). An anatomically comprehensive atlas of the
   adult human brain transcriptome. Nature, 489(7416), 391-399.

3. Markello et al. (2021). Misalignment of AHBA data with MNI space
   impacts interpretations of brain organization. bioRxiv.
