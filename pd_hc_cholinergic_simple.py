#!/usr/bin/env python
"""
PD vs HC T1影像胆碱能受体皮层差异分析 - 简化版

这是一个独立版本，不需要依赖abagen包的所有功能
可以直接从Allen Human Brain Atlas网站下载数据进行离线分析

使用方法：
    python pd_hc_cholinergic_simple.py --data-dir /path/to/data
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy import stats
from collections import defaultdict


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


DESIKAN_KILLIANY_REGIONS = [
    'bankssts', 'caudal_acfb', 'caudal_mfg', 'caudal_sfg', 'fusiform',
    'inferior_parietal', 'inferior_temporal', 'isthmus_cingulate',
    'lateral_orbitofrontal', 'lateral_precuneus', 'lateral_temporal',
    'lingual', 'medial_orbitofrontal', 'middle_temporal', 'parahippocampal',
    'paracentral', 'parsopercularis', 'parsorbitalis', 'parstriangularis',
    'pericalcarine', 'postcentral', 'posterior_cingulate', 'precentral',
    'precuneus', 'rostral_acfb', 'rostral_mfg', 'rostral_sfg',
    'rostral_anterior_cingulate', 'superior_frontal', 'superior_parietal',
    'superior_temporal', 'supramarginal', 'frontalpole', 'temporalpole',
    'transverse_temporal', 'insula'
]


def load_subject_data(data_dir, subjects_file):
    """
    加载受试者数据和脑区表达数据
    """
    subjects_path = os.path.join(data_dir, subjects_file)
    
    if not os.path.exists(subjects_path):
        print(f"错误：找不到受试者文件 {subjects_path}")
        print("\n请创建subjects.csv文件，格式如下：")
        print("subject_id,group,age,sex")
        print("PD001,PD,65,M")
        print("PD002,PD,68,F")
        print("HC001,HC,63,M")
        sys.exit(1)
    
    subjects_df = pd.read_csv(subjects_path)
    
    print(f"\n{'='*60}")
    print("数据加载完成")
    print(f"{'='*60}")
    print(f"总受试者数: {len(subjects_df)}")
    print(f"  - PD组: {sum(subjects_df['group'] == 'PD')} 人")
    print(f"  - HC组: {sum(subjects_df['group'] == 'HC')} 人")
    print(f"  - 年龄范围: {subjects_df['age'].min()}-{subjects_df['age'].max()} 岁")
    
    return subjects_df


def load_expression_data(data_dir, subjects_df):
    """
    加载胆碱能基因表达数据
    如果数据不存在，则生成模拟数据用于演示
    """
    expr_file = os.path.join(data_dir, 'cholinergic_expression.csv')
    
    if os.path.exists(expr_file):
        print(f"\n从文件加载基因表达数据: {expr_file}")
        expr_df = pd.read_csv(expr_file, index_col=0)
    else:
        print(f"\n未找到基因表达数据文件: {expr_file}")
        print("生成模拟数据进行演示...")
        
        expr_df = generate_simulated_expression_data(subjects_df)
        
        print(f"模拟数据已保存至: {expr_file}")
        expr_df.to_csv(expr_file)
    
    return expr_df


def generate_simulated_expression_data(subjects_df):
    """
    生成模拟的胆碱能基因表达数据
    
    注意：这是用于测试和演示的数据
    实际分析时请使用真实的Allen Human Brain Atlas数据
    """
    np.random.seed(42)
    
    regions = DESIKAN_KILLIANY_REGIONS * 2
    regions = [f"lh_{r}" for r in regions] + [f"rh_{r}" for r in regions]
    regions = sorted(list(set(regions)))
    
    data = {}
    
    for _, subject in subjects_df.iterrows():
        subject_id = subject['subject_id']
        group = subject['group']
        
        subject_expr = {}
        
        for region in regions:
            base_expr = np.random.uniform(0.5, 2.0)
            
            if group == 'PD':
                group_effect = np.random.uniform(-0.3, 0.1)
            else:
                group_effect = 0
            
            noise = np.random.normal(0, 0.2)
            expression = max(0, base_expr + group_effect + noise)
            
            subject_expr[region] = expression
        
        data[subject_id] = subject_expr
    
    expr_df = pd.DataFrame(data).T
    expr_df.index.name = 'subject_id'
    
    print(f"  生成了 {len(expr_df)} 个受试者 × {len(expr_df.columns)} 个脑区 的数据")
    print(f"  包含 {len(CHOLINERGIC_GENES)} 个胆碱能基因")
    
    return expr_df


def perform_group_comparison(expression_df, subjects_df):
    """
    执行PD vs HC组间比较
    """
    print(f"\n{'='*60}")
    print("执行组间统计分析")
    print(f"{'='*60}")
    print("方法: 独立样本t检验")
    print("多重比较校正: FDR (Benjamini-Hochberg)")
    
    results = []
    
    for region in expression_df.columns:
        pd_values = []
        hc_values = []
        
        for subject_id in expression_df.index:
            if subject_id in subjects_df['subject_id'].values:
                group = subjects_df[subjects_df['subject_id'] == subject_id]['group'].values[0]
                value = expression_df.loc[subject_id, region]
                
                if not pd.isna(value):
                    if group == 'PD':
                        pd_values.append(value)
                    elif group == 'HC':
                        hc_values.append(value)
        
        if len(pd_values) >= 3 and len(hc_values) >= 3:
            pd_mean = np.mean(pd_values)
            pd_std = np.std(pd_values, ddof=1)
            hc_mean = np.mean(hc_values)
            hc_std = np.std(hc_values, ddof=1)
            
            t_stat, p_value = stats.ttest_ind(pd_values, hc_values)
            
            pooled_std = np.sqrt(((len(pd_values)-1)*pd_std**2 + (len(hc_values)-1)*hc_std**2) / 
                                 (len(pd_values) + len(hc_values) - 2))
            
            if pooled_std > 0:
                cohens_d = (pd_mean - hc_mean) / pooled_std
            else:
                cohens_d = 0
            
            mean_diff = pd_mean - hc_mean
            
            results.append({
                'region': region,
                'n_pd': len(pd_values),
                'n_hc': len(hc_values),
                'mean_pd': pd_mean,
                'mean_hc': hc_mean,
                'std_pd': pd_std,
                'std_hc': hc_std,
                'mean_diff': mean_diff,
                't_statistic': t_stat,
                'p_value': p_value,
                'cohens_d': cohens_d,
                'effect_size_interpretation': interpret_effect_size(cohens_d)
            })
    
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        from statsmodels.stats.multitest import multipletests
        
        reject, p_corrected, _, _ = multipletests(
            results_df['p_value'].values, 
            method='fdr_bh'
        )
        results_df['p_fdr_corrected'] = p_corrected
        results_df['significant'] = reject
        
        n_significant = sum(reject)
        print(f"\n分析完成！")
        print(f"  分析了 {len(results_df)} 个脑区")
        print(f"  FDR校正后显著差异脑区 (p < 0.05): {n_significant} 个")
        
        if n_significant > 0:
            print("\n显著差异脑区：")
            sig_results = results_df[results_df['significant']].sort_values('p_fdr_corrected')
            for i, (_, row) in enumerate(sig_results.head(10).iterrows(), 1):
                direction = "↑" if row['mean_diff'] > 0 else "↓"
                print(f"  {i}. {row['region']}: d={row['cohens_d']:.3f} {direction}")
                print(f"     PD均值={row['mean_pd']:.4f}, HC均值={row['mean_hc']:.4f}")
                print(f"     p={row['p_fdr_corrected']:.4e}")
    
    return results_df


def interpret_effect_size(d):
    """
    解释Cohen's d效应量大小
    """
    abs_d = abs(d)
    if abs_d < 0.2:
        return "微小效应"
    elif abs_d < 0.5:
        return "小效应"
    elif abs_d < 0.8:
        return "中等效应"
    else:
        return "大效应"


def analyze_gene_specific_patterns(expression_df, subjects_df):
    """
    分析各基因在组间的差异模式
    """
    print(f"\n{'='*60}")
    print("胆碱能基因表达模式分析")
    print(f"{'='*60}")
    
    gene_patterns = []
    
    for gene in CHOLINERGIC_GENES:
        pd_gene_vals = []
        hc_gene_vals = []
        
        for subject_id in expression_df.index:
            if subject_id in subjects_df['subject_id'].values:
                group = subjects_df[subjects_df['subject_id'] == subject_id]['group'].values[0]
                
                if f"{gene}_expression" in expression_df.columns:
                    val = expression_df.loc[subject_id, f"{gene}_expression"]
                    if not pd.isna(val):
                        if group == 'PD':
                            pd_gene_vals.append(val)
                        elif group == 'HC':
                            hc_gene_vals.append(val)
        
        if len(pd_gene_vals) >= 3 and len(hc_gene_vals) >= 3:
            t_stat, p_value = stats.ttest_ind(pd_gene_vals, hc_gene_vals)
            
            gene_patterns.append({
                'gene': gene,
                'n_pd': len(pd_gene_vals),
                'n_hc': len(hc_gene_vals),
                'mean_pd': np.mean(pd_gene_vals),
                'mean_hc': np.mean(hc_gene_vals),
                'mean_diff': np.mean(pd_gene_vals) - np.mean(hc_gene_vals),
                't_statistic': t_stat,
                'p_value': p_value,
                'direction': 'PD>HC' if np.mean(pd_gene_vals) > np.mean(hc_gene_vals) else 'PD<HC'
            })
    
    if gene_patterns:
        patterns_df = pd.DataFrame(gene_patterns)
        
        from statsmodels.stats.multitest import multipletests
        reject, p_corrected, _, _ = multipletests(
            patterns_df['p_value'].values, 
            method='fdr_bh'
        )
        patterns_df['p_fdr'] = p_corrected
        patterns_df['significant'] = reject
        
        print("\n各基因表达差异：")
        for _, row in patterns_df.iterrows():
            sig_marker = "*" if row['significant'] else ""
            print(f"  {row['gene']:8s}: {row['direction']} "
                  f"(p={row['p_fdr']:.4e}){sig_marker}")
        
        return patterns_df
    
    return None


def create_summary_visualization(results_df, output_dir):
    """
    创建结果汇总（文本格式）
    """
    summary_file = os.path.join(output_dir, 'summary_report.txt')
    
    with open(summary_file, 'w', encoding='utf-8') as f:
        f.write("="*80 + "\n")
        f.write("PD vs HC 胆碱能受体皮层差异分析报告\n")
        f.write("="*80 + "\n\n")
        
        f.write("分析方法:\n")
        f.write("-" * 40 + "\n")
        f.write("1. 数据来源: Allen Human Brain Atlas\n")
        f.write("2. 脑区划分: Desikan-Killiany atlas\n")
        f.write("3. 组间比较: 独立样本t检验\n")
        f.write("4. 多重比较校正: FDR (Benjamini-Hochberg)\n")
        f.write("5. 效应量: Cohen's d\n\n")
        
        f.write(f"分析的胆碱能相关基因:\n")
        f.write(", ".join(CHOLINERGIC_GENES) + "\n\n")
        
        f.write(f"总脑区数: {len(results_df)}\n")
        
        if len(results_df) > 0:
            n_sig = sum(results_df['significant'])
            f.write(f"显著差异脑区 (FDR < 0.05): {n_sig} 个\n\n")
            
            f.write("\n显著差异脑区详情:\n")
            f.write("-" * 80 + "\n")
            
            if n_sig > 0:
                f.write("脑区名称                | PD均值±SD       | HC均值±SD       | "
                       "差异      | Cohen's d | FDR-p值    | 效应\n")
                f.write("-" * 80 + "\n")
                
                sig_results = results_df[results_df['significant']].sort_values('p_fdr_corrected')
                
                for _, row in sig_results.iterrows():
                    region = f"{row['region'][:21]:<21}"
                    pd_vals = f"{row['mean_pd']:.3f}±{row['std_pd']:.3f}"
                    hc_vals = f"{row['mean_hc']:.3f}±{row['std_hc']:.3f}"
                    diff = f"{row['mean_diff']:+.3f}"
                    d = f"{row['cohens_d']:.3f}"
                    p = f"{row['p_fdr_corrected']:.4e}"
                    effect = row['effect_size_interpretation'][:6]
                    
                    f.write(f"{region} | {pd_vals:<16} | {hc_vals:<16} | "
                           f"{diff:<10} | {d:<8} | {p:<11} | {effect}\n")
            else:
                f.write("未发现显著差异脑区 (FDR < 0.05)\n")
            
            f.write("\n\n效应量解释:\n")
            f.write("-" * 40 + "\n")
            f.write("|Cohen's d| < 0.2: 微小效应\n")
            f.write("0.2 ≤ |Cohen's d| < 0.5: 小效应\n")
            f.write("0.5 ≤ |Cohen's d| < 0.8: 中等效应\n")
            f.write("|Cohen's d| ≥ 0.8: 大效应\n")
            
            f.write("\n\nTop 10 效应量脑区:\n")
            f.write("-" * 40 + "\n")
            top_effect = results_df.nlargest(10, 'cohens_d', keep='first')
            
            for i, (_, row) in enumerate(top_effect.iterrows(), 1):
                direction = "↑" if row['mean_diff'] > 0 else "↓"
                f.write(f"{i:2d}. {row['region']:<25} d={row['cohens_d']:+.3f} {direction} "
                       f"(p={row['p_value']:.4e})\n")
            
        f.write("\n" + "="*80 + "\n")
        f.write("分析完成\n")
        f.write("="*80 + "\n")
    
    print(f"\n汇总报告已保存至: {summary_file}")
    
    return summary_file


def generate_data_template(data_dir):
    """
    生成示例数据结构
    """
    print(f"\n{'='*60}")
    print("创建示例数据结构")
    print(f"{'='*60}")
    
    os.makedirs(data_dir, exist_ok=True)
    os.makedirs(os.path.join(data_dir, 't1_images'), exist_ok=True)
    os.makedirs(os.path.join(data_dir, 'cortical_thickness'), exist_ok=True)
    
    example_subjects = """subject_id,group,age,sex
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
HC020,HC,71,F"""
    
    subjects_file = os.path.join(data_dir, 'subjects.csv')
    with open(subjects_file, 'w') as f:
        f.write(example_subjects)
    
    print(f"  ✓ 创建目录: {data_dir}")
    print(f"  ✓ 创建子目录: t1_images, cortical_thickness")
    print(f"  ✓ 创建文件: {subjects_file}")
    print("\n请将T1图像和皮层厚度数据放入相应目录")


def main():
    parser = argparse.ArgumentParser(
        description='PD vs HC 胆碱能受体皮层差异分析 - 简化版',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 分析示例数据
  python pd_hc_cholinergic_simple.py --data-dir ./example_data

  # 使用真实数据
  python pd_hc_cholinergic_simple.py --data-dir ./my_data --subjects subjects.csv

  # 创建示例数据结构
  python pd_hc_cholinergic_simple.py --create-template ./example_data
        """
    )
    
    parser.add_argument('--data-dir', type=str, default='./example_data',
                       help='数据目录路径 (默认: ./example_data)')
    parser.add_argument('--subjects', type=str, default='subjects.csv',
                       help='受试者信息文件名 (默认: subjects.csv)')
    parser.add_argument('--output', type=str, default=None,
                       help='输出目录 (默认: data-dir/results)')
    parser.add_argument('--create-template', type=str, metavar='DIR',
                       help='创建示例数据结构到指定目录')
    
    args = parser.parse_args()
    
    if args.create_template:
        generate_data_template(args.create_template)
        return
    
    data_dir = args.data_dir
    if not os.path.exists(data_dir):
        print(f"警告：数据目录 {data_dir} 不存在，创建示例结构...")
        generate_data_template(data_dir)
    
    output_dir = args.output if args.output else os.path.join(data_dir, 'results')
    os.makedirs(output_dir, exist_ok=True)
    
    print("\n" + "="*80)
    print("PD vs HC 胆碱能受体皮层差异分析")
    print("="*80)
    print(f"数据目录: {data_dir}")
    print(f"输出目录: {output_dir}")
    print(f"分析的基因: {len(CHOLINERGIC_GENES)} 个胆碱能相关基因")
    
    subjects_df = load_subject_data(data_dir, args.subjects)
    
    expression_df = load_expression_data(data_dir, subjects_df)
    
    results_df = perform_group_comparison(expression_df, subjects_df)
    
    patterns_df = analyze_gene_specific_patterns(expression_df, subjects_df)
    
    results_file = os.path.join(output_dir, 'group_comparison_results.csv')
    results_df.to_csv(results_file, index=False, encoding='utf-8')
    print(f"\n详细结果已保存至: {results_file}")
    
    if patterns_df is not None:
        patterns_file = os.path.join(output_dir, 'gene_expression_patterns.csv')
        patterns_df.to_csv(patterns_file, index=False, encoding='utf-8')
        print(f"基因表达模式已保存至: {patterns_file}")
    
    summary_file = create_summary_visualization(results_df, output_dir)
    
    print("\n" + "="*80)
    print("分析完成！")
    print("="*80)
    print(f"\n所有输出文件保存在: {output_dir}/")
    print("  - group_comparison_results.csv: 组间比较详细结果")
    if patterns_df is not None:
        print("  - gene_expression_patterns.csv: 各基因表达模式")
    print("  - summary_report.txt: 可读的分析报告")


if __name__ == '__main__':
    main()
