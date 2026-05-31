#!/usr/bin/env python
"""
PD vs HC T1影像胆碱能受体皮层差异分析

分析目标：比较帕金森病(PD)患者和健康对照组(HC)在T1影像中
胆碱能受体相关基因的皮层表达差异

使用方法：
    python pd_hc_cholinergic_analysis.py --data-dir /path/to/data --output results/
"""

import os
import argparse
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

try:
    import abagen
except ImportError:
    print("错误：请先安装abagen包：pip install abagen")
    raise

CHOLINERGIC_GENES = [
    'CHAT',      # 胆碱乙酰转移酶
    'SLC5A7',    # 胆碱转运体
    'CHRNA4',    # 烟碱型乙酰胆碱受体α4亚基
    'CHRNA5',    # 烟碱型乙酰胆碱受体α5亚基
    'CHRNA7',    # 烟碱型乙酰胆碱受体α7亚基
    'CHRNB2',    # 烟碱型乙酰胆碱受体β2亚基
    'CHRM1',     # 毒蕈碱型乙酰胆碱受体M1
    'CHRM2',     # 毒蕈碱型乙酰胆碱受体M2
    'CHRM3',     # 毒蕈碱型乙酰胆碱受体M3
    'CHRM4',     # 毒蕈碱型乙酰胆碱受体M4
]


def load_mri_data(data_dir, subjects_file):
    """
    加载MRI数据
    
    参数:
        data_dir: 数据目录路径
        subjects_file: 受试者信息文件（CSV格式，包含subject_id, group列）
    
    返回:
        subjects_df: 受试者信息DataFrame
    """
    subjects_file = os.path.join(data_dir, subjects_file)
    subjects_df = pd.read_csv(subjects_file)
    
    print(f"\n加载了 {len(subjects_df)} 个受试者")
    print(f"  - PD组: {sum(subjects_df['group'] == 'PD')} 人")
    print(f"  - HC组: {sum(subjects_df['group'] == 'HC')} 人")
    
    return subjects_df


def extract_cortical_thickness(data_dir, subjects_df, atlas='desikan_killiany'):
    """
    从T1图像提取皮层厚度
    
    这是一个示例函数，实际使用时需要调用FreeSurfer等工具
    """
    print("\n提取皮层厚度数据...")
    
    cortical_data = {}
    for _, row in subjects_df.iterrows():
        subject_id = row['subject_id']
        group = row['group']
        
        thickness_file = os.path.join(data_dir, 'cortical_thickness', 
                                     f'{subject_id}_thickness.csv')
        
        if os.path.exists(thickness_file):
            thickness_data = pd.read_csv(thickness_file)
            cortical_data[subject_id] = {
                'group': group,
                'thickness': thickness_data
            }
        else:
            print(f"  警告：未找到 {subject_id} 的皮层厚度数据")
    
    return cortical_data


def get_cholinergic_expression(atlas='desikan_killiany'):
    """
    使用abagen获取胆碱能受体基因在各脑区的表达数据
    
    参数:
        atlas: 脑区划分模板名称
    
    返回:
        expression_df: 各脑区的基因表达数据
    """
    print("\n从Allen Human Brain Atlas获取胆碱能受体基因表达数据...")
    
    atlas_file = abagen.fetch_desikan_killiany()
    
    print(f"  使用模板: {atlas}")
    print(f"  分析基因: {', '.join(CHOLINERGIC_GENES)}")
    
    expression_data = abagen.get_expression_data(
        atlas_file,
        probes=CHOLINERGIC_GENES,
        return_counts=True,
        exact_probe=False,
        sample_norm=True,
        probe_norm='spline',
        donors='all',
        tolerance=2.0,
        agg_metric='mean'
    )
    
    print(f"  成功获取 {len(expression_data)} 个脑区的表达数据")
    
    return expression_data


def perform_group_analysis(expression_df, subjects_df, cortical_data=None):
    """
    执行PD vs HC组间分析
    
    参数:
        expression_df: 胆碱能基因表达数据
        subjects_df: 受试者信息
        cortical_data: 皮层厚度数据（可选）
    
    返回:
        results_df: 统计分析结果
    """
    print("\n执行组间统计分析...")
    print("  PD组 vs HC组")
    
    results = []
    
    regions = expression_df.columns
    
    for region in regions:
        if region == 'label' or region == 'structure_id':
            continue
            
        pd_values = []
        hc_values = []
        
        for subject_id in expression_df.index:
            if subject_id in subjects_df['subject_id'].values:
                group = subjects_df[subjects_df['subject_id'] == subject_id]['group'].values[0]
                value = expression_df.loc[subject_id, region]
                
                if not np.isnan(value):
                    if group == 'PD':
                        pd_values.append(value)
                    elif group == 'HC':
                        hc_values.append(value)
        
        if len(pd_values) >= 3 and len(hc_values) >= 3:
            t_stat, p_value = stats.ttest_ind(pd_values, hc_values)
            mean_diff = np.mean(pd_values) - np.mean(hc_values)
            
            effect_size = mean_diff / np.sqrt((np.std(pd_values)**2 + np.std(hc_values)**2) / 2)
            
            results.append({
                'region': region,
                'n_pd': len(pd_values),
                'n_hc': len(hc_values),
                'mean_pd': np.mean(pd_values),
                'mean_hc': np.mean(hc_values),
                'std_pd': np.std(pd_values),
                'std_hc': np.std(hc_values),
                'mean_diff': mean_diff,
                't_statistic': t_stat,
                'p_value': p_value,
                'cohen_d': effect_size
            })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        reject, p_corrected, _, _ = multipletests(
            results_df['p_value'].values, 
            method='fdr_bh'
        )
        results_df['p_fdr_corrected'] = p_corrected
        results_df['significant'] = reject
        
        print(f"\n分析完成！")
        print(f"  总共分析了 {len(results_df)} 个脑区")
        print(f"  显著差异脑区 (FDR < 0.05): {sum(reject)} 个")
    
    return results_df


def correlate_with_cortical_measures(expression_df, cortical_data, subjects_df):
    """
    将基因表达与皮层厚度进行相关分析
    
    这可以帮助理解基因表达与形态学指标的关系
    """
    print("\n执行基因表达与皮层厚度相关分析...")
    
    correlations = []
    
    for gene in CHOLINERGIC_GENES:
        if gene in expression_df.columns:
            for region in expression_df.columns:
                if region in cortical_data:
                    
                    expression_values = []
                    thickness_values = []
                    
                    for subject_id in expression_df.index:
                        if subject_id in cortical_data:
                            expr = expression_df.loc[subject_id, gene]
                            thick = cortical_data[subject_id]['thickness'].get(region)
                            
                            if not np.isnan(expr) and not np.isnan(thick):
                                expression_values.append(expr)
                                thickness_values.append(thick)
                    
                    if len(expression_values) >= 5:
                        r, p = stats.pearsonr(expression_values, thickness_values)
                        
                        correlations.append({
                            'gene': gene,
                            'region': region,
                            'r': r,
                            'p_value': p,
                            'n': len(expression_values)
                        })
    
    corr_df = pd.DataFrame(correlations)
    
    if len(corr_df) > 0:
        reject, p_corrected, _, _ = multipletests(
            corr_df['p_value'].values, 
            method='fdr_bh'
        )
        corr_df['p_fdr_corrected'] = p_corrected
        corr_df['significant'] = reject
        
        print(f"  完成 {len(corr_df)} 个相关分析")
        print(f"  显著相关: {sum(reject)} 个")
    
    return corr_df


def generate_report(results_df, output_dir):
    """
    生成分析报告
    """
    print("\n" + "="*60)
    print("分析报告")
    print("="*60)
    
    report_file = os.path.join(output_dir, 'analysis_report.txt')
    
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("PD vs HC 胆碱能受体皮层差异分析报告\n")
        f.write("="*60 + "\n\n")
        
        f.write("分析方法:\n")
        f.write("-" * 40 + "\n")
        f.write("1. 独立样本t检验比较PD和HC组间差异\n")
        f.write("2. FDR校正进行多重比较\n")
        f.write("3. Cohen's d计算效应大小\n\n")
        
        f.write(f"分析的胆碱能相关基因:\n")
        f.write(", ".join(CHOLINERGIC_GENES) + "\n\n")
        
        if len(results_df) > 0:
            sig_results = results_df[results_df['significant'] == True]
            
            f.write(f"\n显著差异脑区 (FDR < 0.05): {len(sig_results)} 个\n")
            f.write("-" * 40 + "\n")
            
            if len(sig_results) > 0:
                f.write("\n脑区\t\t\t\tPD均值\tHC均值\t差异\tp值\tCohen's d\n")
                f.write("-" * 80 + "\n")
                
                for _, row in sig_results.iterrows():
                    f.write(f"{row['region'][:20]:<20}\t{row['mean_pd']:.4f}\t{row['mean_hc']:.4f}\t"
                           f"{row['mean_diff']:.4f}\t{row['p_fdr_corrected']:.4e}\t{row['cohen_d']:.4f}\n")
            else:
                f.write("未发现显著差异脑区 (FDR < 0.05)\n")
            
            f.write("\n\n效应大小解释 (Cohen's d):\n")
            f.write("-" * 40 + "\n")
            f.write("  |d| < 0.2: 微小效应\n")
            f.write("  |d| ≈ 0.5: 中等效应\n")
            f.write("  |d| ≈ 0.8: 大效应\n\n")
            
            f.write("top 10效应量脑区:\n")
            f.write("-" * 40 + "\n")
            top_effect = results_df.nlargest(10, 'cohen_d', keep='first')
            for _, row in top_effect.iterrows():
                direction = "↑" if row['mean_diff'] > 0 else "↓"
                f.write(f"  {row['region']}: d={row['cohen_d']:.3f} {direction} "
                       f"(PD vs HC, p={row['p_value']:.4e})\n")
        else:
            f.write("未完成有效分析\n")
        
        f.write("\n" + "="*60 + "\n")
        f.write("分析完成\n")
    
    print(f"\n报告已保存至: {report_file}")
    
    return report_file


def main():
    parser = argparse.ArgumentParser(
        description='PD vs HC T1影像胆碱能受体皮层差异分析'
    )
    parser.add_argument('--data-dir', type=str, required=True,
                       help='数据目录路径')
    parser.add_argument('--subjects-file', type=str, default='subjects.csv',
                       help='受试者信息文件')
    parser.add_argument('--output', type=str, default='./results',
                       help='输出目录')
    parser.add_argument('--atlas', type=str, default='desikan_killiany',
                       choices=['desikan_killiany', 'schaefer'],
                       help='脑区划分模板')
    
    args = parser.parse_args()
    
    os.makedirs(args.output, exist_ok=True)
    
    print("\n" + "="*60)
    print("PD vs HC 胆碱能受体皮层差异分析")
    print("="*60)
    
    subjects_df = load_mri_data(args.data_dir, args.subjects_file)
    
    expression_df = get_cholinergic_expression(atlas=args.atlas)
    
    cortical_data = extract_cortical_thickness(args.data_dir, subjects_df)
    
    results_df = perform_group_analysis(expression_df, subjects_df, cortical_data)
    
    if cortical_data:
        corr_df = correlate_with_cortical_measures(expression_df, cortical_data, subjects_df)
        corr_file = os.path.join(args.output, 'correlations.csv')
        corr_df.to_csv(corr_file, index=False)
        print(f"\n相关分析结果已保存至: {corr_file}")
    
    results_file = os.path.join(args.output, 'group_comparison.csv')
    results_df.to_csv(results_file, index=False, encoding='utf-8')
    print(f"组间比较结果已保存至: {results_file}")
    
    report_file = generate_report(results_df, args.output)
    
    print("\n" + "="*60)
    print("分析完成！")
    print("="*60)


if __name__ == '__main__':
    main()
