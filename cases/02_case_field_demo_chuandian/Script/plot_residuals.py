#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
走时残差分析绘图程序 (Travel Time Residual Analysis)
==========================================

本程序用于分析和可视化反演前后的走时残差，包括：
1. P波走时残差（绝对和相对走时）
2. S波走时残差（绝对和相对走时）
3. S-P走时差残差（绝对和相对走时差）

功能特点：
---------
1. 区分绝对和相对走时残差
2. 反演前后残差对比
3. 统计分析（均值、标准差等）
4. 核密度估计曲线
5. 美观的学术风格绘图
6. 输出高质量矢量图和位图

输入文件：
--------
- initial.res：反演前P波残差
- initial_SP.res：反演前S波和S-P残差
- res：反演后P波残差
- res_SP：反演后S波和S-P残差

输出文件：
--------
- residuals_P.pdf/png：P波残差统计图
- residuals_S.pdf/png：S波残差统计图
- residuals_SP.pdf/png：S-P残差统计图

作者信息：
--------
作者：杨浩
单位：中国科学技术大学 地球和空间科学学院
邮箱：youngh_geo@mail.ustc.edu.cn

License:
--------
Copyright (c) 2025 Yang Hao. All rights reserved.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import os

# ========== 参数设置 ==========
outdir = './figure/residuals/'
os.makedirs(outdir, exist_ok=True)

# 颜色设置
colors = {
    'initial_abs': '#FFA07A',  # 浅珊瑚色
    'initial_rel': '#FFB6C1',  # 浅粉色
    'final_abs': '#CD5C5C',    # 印度红
    'final_rel': '#DC143C'     # 深红色
}

# 文件路径
initial_file = "../Output_Files/initial.res"
initial_sp_file = "../Output_Files/initial_SP.res"
final_file = "../Output_Files/res"
final_sp_file = "../Output_Files/res_SP"

def read_residual_file(filename):
    """读取残差文件
    
    Args:
        filename: 残差文件路径
        
    Returns:
        DataFrame: 包含残差数据的DataFrame
    """
    try:
        # 先读取第一行获取列名
        with open(filename, 'r') as f:
            header = f.readline().strip()
            # 处理标题行可能的多空格情况
            header = ' '.join(header.split())
            header = header.split(' ')
        
        # 读取数据行
        data_lines = []
        skipped_lines = 0
        with open(filename, 'r') as f:
            next(f)  # 跳过标题行
            for line in f:
                line = line.strip()
                if not line:  # 跳过空行
                    continue
                    
                # 处理多个空格的情况
                parts = ' '.join(line.split()).split(' ')
                
                if len(parts) >= 9:  # 确保行有足够的列
                    try:
                        # 处理残差值
                        res_val = parts[6]
                        if '*' in res_val:  # 处理 ************ 的情况
                            skipped_lines += 1
                            continue
                            
                        # 创建数据行
                        row = {
                            'STA': parts[0],
                            'DT': float(parts[1]),
                            'C1': int(parts[2]),
                            'C2': int(parts[3]),
                            'IDX': int(parts[4]),
                            'QUAL': float(parts[5]),
                            'RES': float(res_val),
                            'WT': float(parts[7]),
                            'OFFS': float(parts[8])
                        }
                        data_lines.append(row)
                    except (ValueError, IndexError) as e:
                        skipped_lines += 1
                        continue
                else:
                    skipped_lines += 1
        
        if not data_lines:
            print(f"警告: {filename} 中没有有效数据")
            return None
            
        # 创建DataFrame
        df = pd.DataFrame(data_lines)
        print(f"成功读取 {len(df)} 行有效数据从 {filename} (跳过 {skipped_lines} 行无效数据)")
        return df
        
    except Exception as e:
        print(f"读取文件 {filename} 时出错: {str(e)}")
        return None

def process_residuals(df, wave_type=None):
    """处理残差数据，区分绝对和相对走时
    
    Args:
        df: 残差数据DataFrame
        wave_type: 波类型，'P'、'S'或'SP'
        
    Returns:
        tuple: (绝对走时残差, 相对走时残差)
    """
    if df is None or len(df) == 0:
        return np.array([]), np.array([])
        
    try:
        # 根据波类型筛选数据
        if wave_type == 'P':
            df = df[df['IDX'] == 3]  # IDX=3 表示P波
        elif wave_type == 'S':
            df = df[df['IDX'] == 4]  # IDX=4 表示S波
        # SP文件不需要筛选，因为本身就是S-P数据
        
        # 过滤异常值
        df = df[np.abs(df['RES']) < 10000]  # 移除过大的残差值
        
        # 区分绝对和相对走时
        abs_mask = df['C1'] == df['C2']
        abs_res = df[abs_mask]['RES'].values.astype(float)
        rel_res = df[~abs_mask]['RES'].values.astype(float)
        
        if len(abs_res) > 0 or len(rel_res) > 0:
            print(f"{wave_type}波: 绝对走时残差 {len(abs_res)} 个, 相对走时残差 {len(rel_res)} 个")
        
        return abs_res, rel_res
    except Exception as e:
        print(f"处理残差数据时出错: {str(e)}")
        return np.array([]), np.array([])

def plot_residual_histogram(ax, data_initial, data_final, title, xlabel='Residual (ms)'):
    """绘制残差直方图
    
    Args:
        ax: matplotlib轴对象
        data_initial: 反演前残差数据 (abs_res, rel_res)
        data_final: 反演后残差数据 (abs_res, rel_res)
        title: 图标题
        xlabel: x轴标签
    """
    initial_abs, initial_rel = data_initial
    final_abs, final_rel = data_final
    
    # 确保所有数据都是数值型
    all_data = np.concatenate([
        initial_abs.astype(float),
        initial_rel.astype(float),
        final_abs.astype(float),
        final_rel.astype(float)
    ])
    
    # 检查是否有有效数据
    if len(all_data) == 0:
        print(f"警告: {title} 没有有效数据")
        return
    
    # 计算合适的直方图区间
    data_range = np.percentile(all_data[~np.isnan(all_data)], [1, 99])
    bins = np.linspace(data_range[0], data_range[1], 50)
    
    # 绘制直方图和核密度估计曲线
    for data, label, color, alpha in [
        (initial_abs, 'Initial (Absolute)', colors['initial_abs'], 0.6),
        (initial_rel, 'Initial (Relative)', colors['initial_rel'], 0.6),
        (final_abs, 'Final (Absolute)', colors['final_abs'], 0.7),
        (final_rel, 'Final (Relative)', colors['final_rel'], 0.7)
    ]:
        if len(data) > 0:
            # 过滤无效值
            valid_data = data[~np.isnan(data)]
            if len(valid_data) > 0:
                # 计算统计量
                mean = np.mean(valid_data)
                std = np.std(valid_data)
                
                # 直方图
                ax.hist(valid_data, bins=bins, density=True, alpha=alpha, color=color,
                       label=f'{label}\nN={len(valid_data)}\nμ={mean:.2f}, σ={std:.2f}',
                       edgecolor='black', linewidth=0.5)
                
                # 核密度估计
                if len(valid_data) > 1:
                    kde = stats.gaussian_kde(valid_data)
                    x_range = np.linspace(data_range[0], data_range[1], 200)
                    ax.plot(x_range, kde(x_range), '-', color=color, linewidth=1.5)
                    
                    # 添加均值线
                    ax.axvline(mean, color=color, linestyle='--', linewidth=1,
                             alpha=0.8)
    
    # 设置图形样式
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel('Density', fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.legend(fontsize=8, loc='upper right', bbox_to_anchor=(1.02, 1.02))
    
    # 添加零线
    ax.axvline(0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)

def create_residual_plots():
    """创建所有残差图"""
    # 读取数据
    initial_data = read_residual_file(initial_file)
    initial_sp_data = read_residual_file(initial_sp_file)
    final_data = read_residual_file(final_file)
    final_sp_data = read_residual_file(final_sp_file)
    
    # 设置绘图样式
    plt.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['Times New Roman'],
        'font.size': 10,
        'axes.labelsize': 10,
        'axes.titlesize': 12,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 8,
        'figure.dpi': 300,
        'figure.figsize': [15, 5],
        'axes.grid': True,
        'grid.alpha': 0.3,
        'axes.axisbelow': True,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
        'xtick.minor.width': 0.5,
        'ytick.minor.width': 0.5,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'axes.spines.top': True,
        'axes.spines.right': True,
        'axes.spines.left': True,
        'axes.spines.bottom': True,
        'figure.facecolor': 'white',
        'axes.facecolor': 'white'
    })
    
    # 创建图形
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    
    # P波残差
    plot_residual_histogram(ax1,
                          process_residuals(initial_data, wave_type='P'),  # 从P波文件中提取P波数据
                          process_residuals(final_data, wave_type='P'),    # 从P波文件中提取P波数据
                          'P-wave Residuals')
    
    # S波残差
    plot_residual_histogram(ax2,
                          process_residuals(initial_data, wave_type='S'),  # 从P波文件中提取S波数据
                          process_residuals(final_data, wave_type='S'),    # 从P波文件中提取S波数据
                          'S-wave Residuals')
    
    # S-P残差
    plot_residual_histogram(ax3,
                          process_residuals(initial_sp_data, wave_type='SP'),  # SP文件中的数据
                          process_residuals(final_sp_data, wave_type='SP'),    # SP文件中的数据
                          'S-P Residuals')
    
    # 调整布局
    plt.tight_layout()
    
    # 保存图片
    plt.savefig(f'{outdir}residuals.pdf', bbox_inches='tight', dpi=300)
    plt.savefig(f'{outdir}residuals.png', bbox_inches='tight', dpi=300)
    plt.close()
    
    print(f"残差分析图已保存到: {outdir}")

if __name__ == "__main__":
    create_residual_plots() 