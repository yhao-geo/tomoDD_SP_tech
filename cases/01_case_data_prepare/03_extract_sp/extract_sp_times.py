#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
提取S-P到时差程序 (Extract S-P Time Differences)
==============================================

本程序用于从三种不同类型的数据文件中提取S-P波到时差：
1. absolute.dat - 绝对到时数据
2. dt.ct - 走时差数据（catalog）
3. dt.cc - 走时差数据（cross-correlation，可选）

功能特点：
---------
1. 支持多种数据格式的处理
2. 自动提取并匹配P波和S波数据
3. 生成统一格式的输出文件
4. 提供详细的统计分析和可视化
5. 保持4位小数精度

输入文件：
--------
- absolute.dat：绝对到时数据文件
- dt.ct：走时差数据文件（catalog）
- dt.cc：走时差数据文件（cross-correlation，可选）

输出文件：
--------
- absolute_sp.dat：绝对S-P到时差
- dt_sp.ct：走时差S-P到时差（catalog）
- dt_sp.cc：走时差S-P到时差（cross-correlation，如果有输入）
- sp_time_statistics.pdf：统计分析图（矢量格式）
- sp_time_statistics.png：统计分析图（位图格式）

作者信息：
--------
作者：杨浩
单位：中国科学技术大学 地球和空间科学学院
邮箱：youngh_geo@mail.ustc.edu.cn

License:
--------
Copyright (c) 2025 Yang Hao. All rights reserved.
"""

import argparse
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from scipy import stats

class SPTimeExtractor:
    def __init__(self, output_dir: Path):
        """初始化提取器
        
        Args:
            output_dir: 输出目录路径
        """
        self.output_dir = output_dir
        self.output_dir.mkdir(exist_ok=True)
        
        # 存储S-P时差数据
        self.sp_diffs = {
            'absolute': [],
            'dtct': [],
            'dtcc': []
        }
        
        # 设置日志
        self.setup_logging()
        
        # 设置绘图样式
        self.setup_plot_style()
    
    def setup_logging(self):
        """设置日志系统"""
        self.logger = logging.getLogger("SPTimeExtractor")
        self.logger.setLevel(logging.INFO)
        
        # 创建控制台处理器
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        
        # 创建格式化器
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        
        # 添加处理器到日志记录器
        self.logger.addHandler(console_handler)
    
    def setup_plot_style(self):
        """设置matplotlib的绘图样式"""
        plt.rcParams.update({
            # 字体设置
            'font.family': 'sans-serif',
            'font.sans-serif': ['Arial', 'DejaVu Sans'],
            'font.size': 10,
            'axes.labelsize': 10,
            'axes.titlesize': 11,
            'xtick.labelsize': 9,
            'ytick.labelsize': 9,
            'legend.fontsize': 9,
            
            # 图形大小和分辨率
            'figure.dpi': 300,
            'figure.figsize': [7.2, 4.45],  # Nature single column width
            'savefig.dpi': 300,
            'savefig.bbox': 'tight',
            'savefig.pad_inches': 0.1,
            
            # 网格设置
            'axes.grid': True,
            'grid.alpha': 0.3,
            'grid.color': '#b0b0b0',
            'grid.linestyle': '--',
            'grid.linewidth': 0.5,
            
            # 轴线设置
            'axes.linewidth': 0.8,
            'axes.edgecolor': 'black',
            'axes.axisbelow': True,
            
            # 刻度设置
            'xtick.direction': 'out',
            'ytick.direction': 'out',
            'xtick.major.width': 0.8,
            'ytick.major.width': 0.8,
            'xtick.minor.width': 0.5,
            'ytick.minor.width': 0.5,
            
            # 图例设置
            'legend.frameon': True,
            'legend.framealpha': 0.8,
            'legend.edgecolor': 'black',
            'legend.fancybox': False,
            
            # 其他设置
            'axes.unicode_minus': False,
            'pdf.fonttype': 42,
            'ps.fonttype': 42
        })
    
    def process_absolute_data(self, input_file: Path) -> None:
        """处理absolute.dat文件，提取S-P到时差
        
        Args:
            input_file: absolute.dat文件路径
        """
        self.logger.info(f"处理absolute.dat文件: {input_file}")
        output_file = self.output_dir / "absolute_sp.dat"
        
        try:
            current_event = None
            p_data: Dict[str, Tuple[float, float]] = {}
            s_data: Dict[str, Tuple[float, float]] = {}
            
            with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
                for line in fin:
                    line = line.strip()
                    if not line:
                        continue
                        
                    if line.startswith('#'):
                        # 处理前一个事件的数据
                        if p_data and s_data:
                            self._write_sp_differences(fout, current_event, p_data, s_data, 'absolute')
                        
                        # 开始新的事件
                        current_event = line
                        p_data.clear()
                        s_data.clear()
                    else:
                        parts = line.split()
                        if len(parts) >= 4:
                            station = parts[0]
                            travel_time = float(parts[1])
                            weight = float(parts[2])
                            phase = parts[3]
                            
                            if phase == 'P':
                                p_data[station] = (travel_time, weight)
                            elif phase == 'S':
                                s_data[station] = (travel_time, weight)
                
                # 处理最后一个事件的数据
                if p_data and s_data:
                    self._write_sp_differences(fout, current_event, p_data, s_data, 'absolute')
                    
            self.logger.info(f"完成absolute.dat处理，输出文件: {output_file}")
            
        except Exception as e:
            self.logger.error(f"处理absolute.dat文件时出错: {str(e)}")
            raise
    
    def process_dtct_data(self, input_file: Path) -> None:
        """处理dt.ct文件，提取S-P到时差
        
        Args:
            input_file: dt.ct文件路径
        """
        self.logger.info(f"处理dt.ct文件: {input_file}")
        output_file = self.output_dir / "dt_sp.ct"
        
        try:
            current_event = None
            p_data: Dict[str, Tuple[float, float, float]] = {}  # station -> (time1, time2, weight)
            s_data: Dict[str, Tuple[float, float, float]] = {}  # station -> (time1, time2, weight)
            
            with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
                for line in fin:
                    line = line.strip()
                    if not line:
                        continue
                        
                    if line.startswith('#'):
                        # 处理前一个事件的数据
                        if p_data and s_data:
                            self._write_dtct_differences(fout, current_event, p_data, s_data)
                        
                        # 开始新的事件
                        current_event = line
                        p_data.clear()
                        s_data.clear()
                    else:
                        parts = line.split()
                        if len(parts) >= 5:
                            station = parts[0]
                            time1 = float(parts[1])
                            time2 = float(parts[2])
                            weight = float(parts[3])
                            phase = parts[4]
                            
                            if phase == 'P':
                                p_data[station] = (time1, time2, weight)
                            elif phase == 'S':
                                s_data[station] = (time1, time2, weight)
                
                # 处理最后一个事件的数据
                if p_data and s_data:
                    self._write_dtct_differences(fout, current_event, p_data, s_data)
                    
            self.logger.info(f"完成dt.ct处理，输出文件: {output_file}")
            
        except Exception as e:
            self.logger.error(f"处理dt.ct文件时出错: {str(e)}")
            raise
    
    def process_dtcc_data(self, input_file: Path) -> None:
        """处理dt.cc文件，提取S-P到时差
        
        Args:
            input_file: dt.cc文件路径
        """
        self.logger.info(f"处理dt.cc文件: {input_file}")
        output_file = self.output_dir / "dt_sp.cc"
        
        try:
            current_event = None
            p_data: Dict[str, Tuple[float, float]] = {}  # station -> (travel_time, weight)
            s_data: Dict[str, Tuple[float, float]] = {}  # station -> (travel_time, weight)
            
            with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
                for line in fin:
                    line = line.strip()
                    if not line:
                        continue
                        
                    if line.startswith('#'):
                        # 处理前一个事件的数据
                        if p_data and s_data:
                            self._write_sp_differences(fout, current_event, p_data, s_data, 'dtcc')
                        
                        # 开始新的事件
                        current_event = line
                        p_data.clear()
                        s_data.clear()
                    else:
                        parts = line.split()
                        if len(parts) >= 4:
                            station = parts[0]
                            travel_time = float(parts[1])
                            weight = float(parts[2])
                            phase = parts[3]
                            
                            if phase == 'P':
                                p_data[station] = (travel_time, weight)
                            elif phase == 'S':
                                s_data[station] = (travel_time, weight)
                
                # 处理最后一个事件的数据
                if p_data and s_data:
                    self._write_sp_differences(fout, current_event, p_data, s_data, 'dtcc')
                    
            self.logger.info(f"完成dt.cc处理，输出文件: {output_file}")
            
        except Exception as e:
            self.logger.error(f"处理dt.cc文件时出错: {str(e)}")
            raise
    
    def _write_sp_differences(self, fout, event: str, 
                            p_data: Dict[str, Tuple[float, float]], 
                            s_data: Dict[str, Tuple[float, float]], 
                            data_type: str = 'absolute') -> None:
        """写入S-P到时差（用于absolute.dat和dt.cc）
        
        Args:
            fout: 输出文件对象
            event: 事件信息行
            p_data: P波数据字典
            s_data: S波数据字典
            data_type: 数据类型，'absolute'或'dtcc'
        """
        if not s_data:  # 如果没有S波数据，跳过该事件
            return
            
        # 写入事件信息
        fout.write(f"{event}\n")
        
        # 对于每个同时有P波和S波数据的台站，计算并写入时差
        for station in s_data:
            if station in p_data:
                s_time, s_weight = s_data[station]
                p_time, p_weight = p_data[station]
                time_diff = s_time - p_time
                avg_weight = (s_weight + p_weight) / 2
                fout.write(f"{station} {time_diff:.4f} {avg_weight}\n")
                
                # 存储时差数据用于统计
                self.sp_diffs[data_type].append(time_diff)
    
    def _write_dtct_differences(self, fout, event: str,
                              p_data: Dict[str, Tuple[float, float, float]],
                              s_data: Dict[str, Tuple[float, float, float]]) -> None:
        """写入S-P到时差（用于dt.ct）
        
        Args:
            fout: 输出文件对象
            event: 事件信息行
            p_data: P波数据字典
            s_data: S波数据字典
        """
        if not s_data:  # 如果没有S波数据，跳过该事件
            return
            
        # 写入事件信息
        fout.write(f"{event}\n")
        
        # 对于每个同时有P波和S波数据的台站，计算并写入时差
        for station in s_data:
            if station in p_data:
                s_time1, s_time2, s_weight = s_data[station]
                p_time1, p_time2, p_weight = p_data[station]
                dt = (s_time1 - s_time2) - (p_time1 - p_time2)
                avg_weight = (s_weight + p_weight) / 2
                fout.write(f"{station} {dt:.4f} {avg_weight}\n")
                
                # 存储时差数据用于统计
                self.sp_diffs['dtct'].append(dt)

    def plot_statistics(self):
        """绘制S-P时差的统计图"""
        self.logger.info("开始绘制统计图...")
        
        # 创建图形
        fig = plt.figure(figsize=(15, 5))
        
        # 为每种数据类型创建子图
        for i, (data_type, data) in enumerate(self.sp_diffs.items(), 1):
            if not data:  # 跳过空数据
                continue
                
            ax = fig.add_subplot(1, 3, i)
            
            # 计算统计量
            mean_val = np.mean(data)
            median_val = np.median(data)
            std_val = np.std(data)
            
            # 创建柱状图
            n, bins, patches = ax.hist(data, bins='auto', density=True,
                                     alpha=0.7, color='#2196F3',
                                     edgecolor='black', linewidth=0.5)
            
            # 添加核密度估计曲线
            kde = stats.gaussian_kde(data)
            x_range = np.linspace(min(data), max(data), 200)
            ax.plot(x_range, kde(x_range), 'r-', linewidth=1.5, 
                   label='KDE', color='#E91E63')
            
            # 添加均值和中位数线
            ax.axvline(mean_val, color='#4CAF50', linestyle='--', linewidth=1,
                      label=f'Mean: {mean_val:.2f}')
            ax.axvline(median_val, color='#FFC107', linestyle=':', linewidth=1,
                      label=f'Median: {median_val:.2f}')
            
            # 设置标题和标签
            title = {
                'absolute': 'Absolute S-P Time',
                'dtct': 'Double-Difference S-P Time (CT)',
                'dtcc': 'Double-Difference S-P Time (CC)'
            }
            ax.set_title(title[data_type])
            ax.set_xlabel('Time Difference (s)')
            ax.set_ylabel('Density')
            
            # 添加统计信息文本框
            stats_text = (f'N = {len(data)}\n'
                         f'Mean = {mean_val:.2f}\n'
                         f'Median = {median_val:.2f}\n'
                         f'Std = {std_val:.2f}')
            ax.text(0.95, 0.95, stats_text,
                   transform=ax.transAxes,
                   verticalalignment='top',
                   horizontalalignment='right',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
            
            # 优化刻度
            ax.xaxis.set_major_locator(MaxNLocator(5))
            ax.yaxis.set_major_locator(MaxNLocator(5))
            
            # 添加图例
            ax.legend(loc='upper left', frameon=True, fancybox=False, edgecolor='black')
        
        # 调整布局
        plt.tight_layout()
        
        # 保存图形
        plot_file = self.output_dir / 'sp_time_statistics.pdf'
        plt.savefig(plot_file, bbox_inches='tight', dpi=300)
        self.logger.info(f"统计图已保存到: {plot_file}")
        
        # 同时保存PNG格式用于快速预览
        plot_file_png = self.output_dir / 'sp_time_statistics.png'
        plt.savefig(plot_file_png, bbox_inches='tight', dpi=300)
        
        plt.close()

def main():
    """主函数：直接设置输入输出文件路径"""
    # 设置输入文件路径
    absolute_file = Path("./input/absolute.dat")
    dtct_file = Path("./input/dt.ct")
    dtcc_file = Path("./input/dt.cc")  # 可选文件
    
    # 设置输出目录
    output_dir = Path("./output")  # 移除末尾的斜杠，使用Path对象
    
    # 创建输入目录（如果不存在）
    absolute_file.parent.mkdir(exist_ok=True)
    
    # 创建处理器
    extractor = SPTimeExtractor(output_dir)
    
    try:
        # 处理absolute.dat
        if absolute_file.exists():
            extractor.process_absolute_data(absolute_file)
        else:
            extractor.logger.error(f"文件不存在: {absolute_file}")
            return 1
        
        # 处理dt.ct
        if dtct_file.exists():
            extractor.process_dtct_data(dtct_file)
        else:
            extractor.logger.error(f"文件不存在: {dtct_file}")
            return 1
        
        # 如果dt.cc文件存在，则处理它
        if dtcc_file.exists():
            extractor.process_dtcc_data(dtcc_file)
        else:
            extractor.logger.info(f"dt.cc文件不存在，跳过处理: {dtcc_file}")
        
        # 绘制统计图
        extractor.plot_statistics()
        
    except Exception as e:
        extractor.logger.error(f"处理过程中出错: {str(e)}")
        return 1
    
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main()) 