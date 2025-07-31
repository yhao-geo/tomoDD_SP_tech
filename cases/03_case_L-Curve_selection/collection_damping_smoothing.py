#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
阻尼和平滑参数选择程序 (Damping and Smoothing Parameter Selection)
=======================================================

本程序用于tomoDD的阻尼和平滑参数选择，主要功能包括：
1. 并行计算不同参数组合
2. 自动修改输入文件
3. 收集计算结果
4. 生成L曲线数据

功能特点：
---------
1. 支持多进程并行计算
2. 自动管理计算任务
3. 实时监控计算进度
4. 错误处理和日志记录

输入文件：
--------
- tomoSPD_SP.inp：反演控制文件

输出文件：
--------
- collection_damp.dat：参数选择结果
- term_smooth*_damping*：计算日志文件

参数设置：
--------
- damping_values：阻尼因子列表
- smoothing_values：平滑因子列表
- max_workers：最大并行进程数

作者信息：
--------
作者：杨浩
单位：中国科学技术大学 地球和空间科学学院
邮箱：youngh_geo@mail.ustc.edu.cn

License:
--------
Copyright (c) 2025 Yang Hao. All rights reserved.
"""

import os
import sys
import subprocess
import multiprocessing as mp
from pathlib import Path
import logging
import time
import numpy as np
from typing import List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.ticker import MaxNLocator

class DampingSmoothingSelection:
    def __init__(self, plot_only=False):
        """初始化参数选择程序
        
        Args:
            plot_only (bool): 如果为True，则仅执行绘图功能，不进行计算
        """
        # 运行模式设置
        self.plot_only = plot_only
        
        # 设置参数范围
        self.damping_values = [1, 5, 10, 20, 30, 40, 50, 100]
        self.smoothing_values = [15, 20, 30, 40, 50, 100]
        
        # 程序设置
        self.inp_file = "tomoDD_SP.inp"
        self.program = "./tomoDD_SP"
        self.output_file = "collection_damp.dat"
        self.max_workers = 10  # 最大并行进程数
        
        # 图形输出设置
        self.figure_dir = "figures"
        os.makedirs(self.figure_dir, exist_ok=True)
        
        # 设置matplotlib风格
        plt.rcParams.update({
            'font.family': 'serif',
            'font.serif': ['Times New Roman'],
            'font.size': 12,
            'axes.labelsize': 14,
            'axes.titlesize': 16,
            'xtick.labelsize': 12, 
            'ytick.labelsize': 12,
            'legend.fontsize': 12,
            'axes.linewidth': 1.5,
            'lines.linewidth': 2,
            'lines.markersize': 8,
            'pdf.fonttype': 42,
            'ps.fonttype': 42
        })
        
        # 设置日志
        self.setup_logging()
        
        # 检查文件和程序
        self.check_files()
    
    def setup_logging(self):
        """设置日志系统"""
        self.logger = logging.getLogger("DampingSmoothingSelection")
        self.logger.setLevel(logging.INFO)
        
        # 创建控制台处理器
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        
        # 创建格式化器
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        
        # 添加处理器到日志记录器
        self.logger.addHandler(console_handler)
    
    def check_files(self):
        """检查必要的文件和程序"""
        if not os.path.exists(self.inp_file):
            raise FileNotFoundError(f"输入文件 {self.inp_file} 不存在")
        if not os.path.exists(self.program):
            raise FileNotFoundError(f"程序 {self.program} 不存在")
        
        # 设置程序执行权限
        os.chmod(self.program, 0o755)
        
        # 仅在非plot_only模式下清理旧文件
        if not self.plot_only:
            # 清理旧文件
            if os.path.exists(self.output_file):
                os.remove(self.output_file)
            
            # 清理旧的输出目录
            for f in Path('.').glob('term_smooth*_damping*'):
                f.unlink()
    
    def modify_inp_file(self, damping: float, smoothing: float, temp_id: int) -> str:
        """修改输入文件中的参数
        
        Args:
            damping: 阻尼因子
            smoothing: 平滑因子
            temp_id: 临时文件ID
            
        Returns:
            str: 临时输入文件路径
        """
        # 创建临时文件名
        temp_inp = f"{self.inp_file}_temp_{temp_id}"
        
        # 读取文件内容
        with open(self.inp_file, 'r') as f:
            lines = f.readlines()
        
        # 修改阻尼行（第85行）
        damping_pattern = f"    1     0.01  0.01   -9   -9  0.1  0.08      6  -9  10.0 {damping}   1    0.05    0.05\n"
        lines[84] = damping_pattern
        
        # 修改平滑行（第74行）
        smoothing_pattern = f"{smoothing}  " * 8 + f"{smoothing}\n"
        lines[73] = smoothing_pattern
        
        # 写入临时文件
        with open(temp_inp, 'w') as f:
            f.writelines(lines)
            
        return temp_inp
    
    def run_single_case(self, damping: float, smoothing: float, temp_id: int) -> Tuple[float, float, str]:
        """运行单个参数组合的计算
        
        Args:
            damping: 阻尼因子
            smoothing: 平滑因子
            temp_id: 临时文件ID
            
        Returns:
            Tuple[float, float, str]: 阻尼因子，平滑因子，计算结果
        """
        try:
            # 创建临时输入文件
            temp_inp = self.modify_inp_file(damping, smoothing, temp_id)
            
            # 设置输出文件名
            output = f"term_smooth{smoothing}_damping{damping}"
            
            # 运行程序
            result = subprocess.run(
                [self.program, temp_inp],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            
            # 保存输出
            with open(output, 'w') as f:
                f.write(result.stdout)
            
            # 删除临时输入文件
            os.remove(temp_inp)
            
            # 提取结果
            for line in result.stdout.split('\n'):
                if "smooth,damp,xnorm" in line:
                    values = line.split()[4:10]  # 提取需要的数值
                    return damping, smoothing, " ".join(values)
            
            return damping, smoothing, "ERROR: No results found"
            
        except Exception as e:
            self.logger.error(f"计算出错 (damping={damping}, smoothing={smoothing}): {str(e)}")
            if os.path.exists(temp_inp):
                os.remove(temp_inp)
            return damping, smoothing, f"ERROR: {str(e)}"
    
    def run_parallel(self):
        """并行运行所有参数组合"""
        # 创建参数组合
        combinations = [(d, s) for d in self.damping_values for s in self.smoothing_values]
        total_cases = len(combinations)
        
        self.logger.info(f"开始计算 {total_cases} 个参数组合")
        start_time = time.time()
        
        # 创建进程池
        with mp.Pool(processes=self.max_workers) as pool:
            # 提交所有任务
            results = []
            for i, (damping, smoothing) in enumerate(combinations):
                result = pool.apply_async(self.run_single_case, (damping, smoothing, i))
                results.append(result)
            
            # 收集结果
            with open(self.output_file, 'w') as f:
                completed = 0
                for result in results:
                    damping, smoothing, values = result.get()
                    completed += 1
                    self.logger.info(f"完成进度: {completed}/{total_cases} "
                                   f"(damping={damping}, smoothing={smoothing})")
                    f.write(f"{smoothing} {damping} {values}\n")
        
        end_time = time.time()
        self.logger.info(f"计算完成，用时 {end_time - start_time:.2f} 秒")

    def plot_lcurve(self):
        """绘制L曲线"""
        try:
            # 读取数据
            data = pd.read_csv(self.output_file, sep=' ', header=None,
                             names=['smooth', 'damp', 'xnorm', 'xnorm_vel', 'rnorm', 'rnorm_wt'])
            
            # 全局归一化
            xnorm_global = (data['xnorm'] - data['xnorm'].min()) / (data['xnorm'].max() - data['xnorm'].min())
            rnorm_global = (data['rnorm'] - data['rnorm'].min()) / (data['rnorm'].max() - data['rnorm'].min())
            rnorm_wt_global = (data['rnorm_wt'] - data['rnorm_wt'].min()) / (data['rnorm_wt'].max() - data['rnorm_wt'].min())
            
            # 为每个平滑参数绘制一条L曲线
            unique_smooths = sorted(data['smooth'].unique())
            colors = plt.cm.viridis(np.linspace(0, 1, len(unique_smooths)))
            
            # 创建图形
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
            
            # 绘制解范数 vs 绝对残差范数
            for smooth, color in zip(unique_smooths, colors):
                subset = data[data['smooth'] == smooth]
                subset_indices = subset.index
                
                # 使用全局归一化的数据
                x_norm = xnorm_global[subset_indices]
                r_norm = rnorm_global[subset_indices]
                
                # 按照damping值排序
                sort_idx = subset['damp'].argsort()
                x_norm = x_norm.iloc[sort_idx]
                r_norm = r_norm.iloc[sort_idx]
                damping_values = subset['damp'].iloc[sort_idx]
                
                # 绘制L曲线
                line = ax1.plot(x_norm, r_norm, '-o', color=color, 
                              label=f'Smooth={smooth}', zorder=2)
                
                # 添加damping标注
                for x, y, damp in zip(x_norm, r_norm, damping_values):
                    ax1.annotate(f'{damp}', (x, y), xytext=(5, 5),
                               textcoords='offset points', fontsize=8)
            
            ax1.set_xlabel('Normalized Solution Norm', fontweight='bold')
            ax1.set_ylabel('Normalized Absolute Residual Norm', fontweight='bold')
            ax1.set_title('L-curve (Absolute Residual)', fontweight='bold')
            ax1.grid(True, linestyle='--', alpha=0.3)
            ax1.set_aspect('equal')
            ax1.set_xlim(-0.05, 1.05)
            ax1.set_ylim(-0.05, 1.05)
            ax1.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
            
            # 绘制解范数 vs 加权残差范数
            for smooth, color in zip(unique_smooths, colors):
                subset = data[data['smooth'] == smooth]
                subset_indices = subset.index
                
                # 使用全局归一化的数据
                x_norm = xnorm_global[subset_indices]
                r_norm_wt = rnorm_wt_global[subset_indices]
                
                # 按照damping值排序
                sort_idx = subset['damp'].argsort()
                x_norm = x_norm.iloc[sort_idx]
                r_norm_wt = r_norm_wt.iloc[sort_idx]
                damping_values = subset['damp'].iloc[sort_idx]
                
                # 绘制L曲线
                line = ax2.plot(x_norm, r_norm_wt, '-o', color=color,
                              label=f'Smooth={smooth}', zorder=2)
                
                # 添加damping标注
                for x, y, damp in zip(x_norm, r_norm_wt, damping_values):
                    ax2.annotate(f'{damp}', (x, y), xytext=(5, 5),
                               textcoords='offset points', fontsize=8)
            
            ax2.set_xlabel('Normalized Solution Norm', fontweight='bold')
            ax2.set_ylabel('Normalized Weighted Residual Norm', fontweight='bold')
            ax2.set_title('L-curve (Weighted Residual)', fontweight='bold')
            ax2.grid(True, linestyle='--', alpha=0.3)
            ax2.set_aspect('equal')
            ax2.set_xlim(-0.05, 1.05)
            ax2.set_ylim(-0.05, 1.05)
            ax2.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
            
            # 调整布局并保存
            plt.tight_layout()
            fig.savefig(os.path.join(self.figure_dir, 'lcurve_all.pdf'), dpi=300, bbox_inches='tight')
            fig.savefig(os.path.join(self.figure_dir, 'lcurve_all.png'), dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info("L曲线绘制完成")
            
        except Exception as e:
            self.logger.error(f"绘制L曲线时出错: {str(e)}")
            raise

    def run(self):
        """运行参数选择程序"""
        try:
            if not self.plot_only:
                # 运行并行计算
                self.run_parallel()
            
            # 检查结果文件是否存在
            if not os.path.exists(self.output_file):
                self.logger.error(f"结果文件 {self.output_file} 不存在")
                return 1
            
            # 绘制L曲线
            self.plot_lcurve()
            
            if not self.plot_only:
                os.system('rm -rf term*')
            
            return 0
            
        except Exception as e:
            self.logger.error(f"程序运行出错: {str(e)}")
            return 1

def main():
    """主函数"""
    # 从命令行参数获取运行模式
    plot_only = "--plot-only" in sys.argv
    selection = DampingSmoothingSelection(plot_only=plot_only)
    return selection.run()

if __name__ == "__main__":
    sys.exit(main()) 
