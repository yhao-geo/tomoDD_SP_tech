#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
震相筛选程序 (Phase Selection Program)
================================

本程序用于对地震数据进行震相筛选，基于走时-距离关系进行数据质量控制。

功能特点：
---------
1. 自动计算震源-台站距离
2. 绘制走时-距离关系图
3. 支持自动/手动拟合走时曲线
4. 多轮sigma剔除的稳健拟合
5. 自动识别和剔除异常数据
6. 支持从缓存加载数据以提高效率
7. 提供详细的日志记录

输入文件：
--------
- station.dat：台站位置文件
  格式：台站名 纬度 经度 高程
- phase.dat：震相数据文件
  格式：事件信息（以#开头）和震相数据

输出文件：
--------
- phase.dat_selection：筛选后的震相数据
- t_dist.dat：时间-距离数据（缓存文件）
- t_dist.png：走时-距离关系图
- processing.log：处理日志
- filtering_stats.txt：数据筛选统计信息

参数设置：
--------
- 拟合方式：'auto'（自动）或'manual'（手动）
- P波边界：b1_P（上边界），b2_P（下边界）
- S波边界：b1_S（上边界），b2_S（下边界）
- 默认边界值：P波±4秒，S波±5秒

使用方法：
--------
1. 准备输入文件：station.dat和phase.dat
2. 设置参数（可选）：
   - fitting_method：'auto'或'manual'
   - boundary_params：边界参数
3. 运行程序：python pha_selection.py

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
import matplotlib
matplotlib.use('Agg')
from pathlib import Path
import math
import logging
from tqdm import tqdm
import sys
import time
import random
from matplotlib.patches import Rectangle

class SeismicDataProcessor:
    def __init__(self, station_file, phase_file, output_dir="./output", log_level=logging.INFO, load_from_cache=False):
        """初始化处理器
        Args:
            station_file (str): 台站文件路径
            phase_file (str): 相位文件路径
            output_dir (str): 输出目录
            log_level (int): 日志级别
            load_from_cache (bool): 是否从缓存文件加载时间-距离数据
        """
        self.station_file = station_file
        self.phase_file = phase_file
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.load_from_cache = load_from_cache
        self.cache_file = self.output_dir / 't_dist.dat'
        
        # 设置日志
        self.setup_logging(log_level)
        
        # 数据
        self.stations = None
        self.phases = None
        self.time_dist_data = None
        
        # 拟合参数（包含斜率、截距和边界参数）
        self.params = {
            'slope_P': None, 'b_P': None,
            'slope_S': None, 'b_S': None,
            'b1_P': 1.5, 'b2_P': 1.5,  # 默认边界参数
            'b1_S': 1.5, 'b2_S': 1.5   # 默认边界参数
        }
        
        self.logger.info(f"初始化完成: 台站文件={station_file}, 相位文件={phase_file}, 输出目录={output_dir}, 从缓存加载={load_from_cache}")
    
    def setup_logging(self, log_level):
        """设置日志系统"""
        self.logger = logging.getLogger("SeismicDataProcessor")
        self.logger.setLevel(log_level)
        
        # 创建控制台处理器
        console_handler = logging.StreamHandler()
        console_handler.setLevel(log_level)
        
        # 创建文件处理器
        file_handler = logging.FileHandler(self.output_dir / "processing.log")
        file_handler.setLevel(log_level)
        
        # 创建格式化器
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        file_handler.setFormatter(formatter)
        
        # 添加处理器到日志记录器
        self.logger.addHandler(console_handler)
        self.logger.addHandler(file_handler)
    
    def load_station_data(self):
        """加载台站数据"""
        self.logger.info("开始加载台站数据")
        
        try:
            self.stations = pd.read_csv(self.station_file, 
                                       sep=r'\s+', 
                                       names=['station_id', 'latitude', 'longitude', 'elevation'])
            # 转换深度单位
            self.stations['depth'] = -self.stations['elevation'] / 1000
            
            self.logger.info(f"成功加载台站数据: {len(self.stations)}个台站")
            return self.stations
            
        except Exception as e:
            self.logger.error(f"加载台站数据失败: {str(e)}")
            raise
    
    def load_phase_data(self):
        """加载相位数据"""
        self.logger.info("开始加载相位数据")
        
        try:
            # 需要特殊处理包含#开头注释行的数据
            with open(self.phase_file, 'r') as f:
                lines = f.readlines()
            
            self.logger.info(f"读取了{len(lines)}行相位数据文件")
            
            events = []
            phases = []
            current_event = None
            
            # 使用tqdm添加进度条
            for line in tqdm(lines, desc="解析相位数据"):
                if line.startswith('#'):
                    # 解析事件信息
                    parts = line.strip().split()
                    # 从原始行中提取事件ID
                    year = parts[1]
                    month = parts[2].zfill(2)
                    day = parts[3].zfill(2)
                    hour = parts[4].zfill(2)
                    minute = parts[5].zfill(2)
                    second = f"{float(parts[6]):06.2f}".zfill(6)
                    event_id = f"{year}{month}{day}{hour}{minute}{second}"
                    current_event = event_id
                    
                    event_info = {
                        'event_id': event_id,
                        'latitude': float(parts[7]),
                        'longitude': float(parts[8]),
                        'depth': float(parts[9]),
                        'line': line.strip()
                    }
                    events.append(event_info)
                else:
                    # 解析相位信息
                    parts = line.strip().split()
                    if len(parts) >= 4:
                        phase_info = {
                            'station_id': parts[0],
                            'travel_time': float(parts[1]),
                            'phase_type': parts[3],
                            'event_id': current_event,
                            'event_lat': events[-1]['latitude'],
                            'event_lon': events[-1]['longitude'],
                            'event_depth': events[-1]['depth'],
                            'line': line.strip()
                        }
                        phases.append(phase_info)
            
            self.phases = pd.DataFrame(phases)
            self.events = pd.DataFrame(events)
            
            self.logger.info(f"成功加载相位数据: {len(self.events)}个事件, {len(self.phases)}个相位")
            return self.phases
            
        except Exception as e:
            self.logger.error(f"加载相位数据失败: {str(e)}")
            raise
    
    def calculate_distance(self):
        """计算震源到台站的距离"""
        self.logger.info("开始处理震源距离数据")
        
        # 检查是否从缓存加载
        if self.load_from_cache and self.cache_file.exists():
            try:
                self.logger.info(f"从缓存文件加载时间-距离数据: {self.cache_file}")
                self.time_dist_data = pd.read_csv(self.cache_file, sep='\s+', 
                                                names=['travel_time', 'distance', 'phase_type'])
                
                # 确保phases和stations数据已加载
                if self.stations is None:
                    self.load_station_data()
                if self.phases is None:
                    self.load_phase_data()
                
                # 重建phases_with_dist数据
                self.phases_with_dist = pd.merge(
                    self.phases,
                    self.stations[['station_id', 'latitude', 'longitude', 'depth']],
                    on='station_id',
                    how='left',
                    suffixes=('', '_station')
                )
                self.phases_with_dist['distance'] = self.time_dist_data['distance'].values
                
                self.logger.info(f"成功从缓存加载数据: {len(self.time_dist_data)}个数据点")
                return self.time_dist_data
            except Exception as e:
                self.logger.warning(f"从缓存加载失败: {str(e)}，将重新计算距离")
        
        if self.stations is None:
            self.logger.info("台站数据未加载，正在加载...")
            self.load_station_data()
            
        if self.phases is None:
            self.logger.info("相位数据未加载，正在加载...")
            self.load_phase_data()
        
        try:
            # 合并台站信息到相位数据
            self.logger.info("合并台站信息到相位数据...")
            phases_with_stations = pd.merge(
                self.phases, 
                self.stations[['station_id', 'latitude', 'longitude', 'depth']], 
                on='station_id',
                how='left',
                suffixes=('', '_station')
            )
            
            # 调试信息：打印列名和数据样例
            self.logger.info(f"合并后的数据列名: {phases_with_stations.columns.tolist()}")
            self.logger.info("数据样例（前1行）:")
            self.logger.info(phases_with_stations.iloc[0].to_string())
            
            # 检查合并后是否有丢失数据
            missing_stations = phases_with_stations[phases_with_stations['latitude'].isna()]['station_id'].unique()
            if len(missing_stations) > 0:
                self.logger.warning(f"有{len(missing_stations)}个台站在台站文件中找不到: {', '.join(missing_stations[:5])}{'...' if len(missing_stations) > 5 else ''}")
                
            # 过滤掉包含NaN值的数据
            valid_data = phases_with_stations.dropna(subset=['latitude', 'longitude', 'depth'])
            dropped_count = len(phases_with_stations) - len(valid_data)
            self.logger.info(f"过滤掉{dropped_count}条包含NaN值的数据记录")
            phases_with_stations = valid_data
            
            # 计算距离
            self.logger.info("计算震源距离...")
            def calc_hypo_dist(row):
                try:
                    dlat = float(row['event_lat']) - float(row['latitude'])
                    dlon = float(row['event_lon']) - float(row['longitude'])
                    # 使用与awk脚本相同的π值和计算方式
                    PI = 3.1415926
                    dist = math.sqrt(
                        (dlat * 111)**2 + 
                        (dlon * math.cos(float(row['event_lat']) * PI / 180) * 111)**2 + 
                        (float(row['event_depth']) - float(row['depth']))**2
                    )
                    return dist
                except Exception as e:
                    self.logger.error(f"距离计算错误，行数据: {row.to_dict()}")
                    self.logger.error(f"错误信息: {str(e)}")
                    return np.nan
            
            tqdm.pandas(desc="计算距离")
            phases_with_stations['distance'] = phases_with_stations.progress_apply(calc_hypo_dist, axis=1)
            
            # 生成时间-距离数据
            time_dist_data = phases_with_stations[['travel_time', 'distance', 'phase_type']].copy()
            # 将震相类型转换为数字（P=1, S=2）
            time_dist_data.loc[:, 'phase_type'] = time_dist_data['phase_type'].map({'P': 1, 'S': 2})
            
            self.phases_with_dist = phases_with_stations
            self.time_dist_data = time_dist_data
            
            # 保存时间-距离数据到缓存文件
            self.logger.info(f"保存时间-距离数据到缓存文件: {self.cache_file}")
            self.time_dist_data.to_csv(self.cache_file, sep=' ', index=False, header=False)
            
            return self.time_dist_data
            
        except Exception as e:
            self.logger.error(f"处理震源距离数据失败: {str(e)}")
            raise
    
    @staticmethod
    def robust_linear_fit(x, y, n_rounds=3, sigma=1):
        """
        多轮sigma剔除的线性拟合，返回最终斜率、截距和掩码
        """
        x = np.asarray(x)
        y = np.asarray(y)
        
        # 检查数据有效性
        valid_mask = ~(np.isnan(x) | np.isnan(y) | np.isinf(x) | np.isinf(y))
        x = x[valid_mask]
        y = y[valid_mask]
        
        # 确保有足够的数据点进行拟合
        if len(x) < 2:
            raise ValueError("Not enough valid data points for fitting")
            
        mask = np.ones(len(x), dtype=bool)
        for _ in range(n_rounds):
            x_fit = x[mask]
            y_fit = y[mask]
            if len(x_fit) < 2:
                break
            try:
                coef = np.polyfit(x_fit, y_fit, 1)
                y_pred = np.polyval(coef, x_fit)
                residuals = y_fit - y_pred
                std = residuals.std()
                outlier_mask = np.abs(residuals) > sigma * std
                if not np.any(outlier_mask):
                    break
                mask[np.where(mask)[0][outlier_mask]] = False
            except Exception as e:
                raise ValueError(f"Fitting failed: {str(e)}")
        
        # 如果所有点都被剔除了，使用原始数据进行拟合
        if not np.any(mask):
            mask = np.ones(len(x), dtype=bool)
            coef = np.polyfit(x, y, 1)
            
        return coef, mask

    def analyze_time_distance(self):
        """分析时间距离关系并拟合参数（多轮剔除离群点）"""
        self.logger.info("开始分析时间-距离关系")
        if self.time_dist_data is None:
            self.logger.info("时间-距离数据未计算，正在计算...")
            self.calculate_distance()
        try:
            # 分离P波和S波数据
            p_data = self.time_dist_data[self.time_dist_data['phase_type'] == 1][['travel_time', 'distance']]
            s_data = self.time_dist_data[self.time_dist_data['phase_type'] == 2][['travel_time', 'distance']]
            self.logger.info(f"P波数据点: {len(p_data)}个, S波数据点: {len(s_data)}个")
            # 拟合P波，自动剔除离群点
            if len(p_data) > 1:
                coef_p, mask_p = self.robust_linear_fit(p_data['distance'].values, p_data['travel_time'].values, n_rounds=3, sigma=3)
                self.params['slope_P'] = coef_p[0]
                self.params['b_P'] = coef_p[1]
                self.logger.info(f'P波拟合结果: slope_P={self.params["slope_P"]:.6f}, b_P={self.params["b_P"]:.6f}, 剩余点数: {mask_p.sum()}')
            else:
                self.logger.warning("P波数据点不足以进行拟合")
            # 拟合S波，自动剔除离群点
            if len(s_data) > 1:
                coef_s, mask_s = self.robust_linear_fit(s_data['distance'].values, s_data['travel_time'].values, n_rounds=3, sigma=3)
                self.params['slope_S'] = coef_s[0]
                self.params['b_S'] = coef_s[1]
                self.logger.info(f'S波拟合结果: slope_S={self.params["slope_S"]:.6f}, b_S={self.params["b_S"]:.6f}, 剩余点数: {mask_s.sum()}')
            else:
                self.logger.warning("S波数据点不足以进行拟合")
            # 绘图
            self.logger.info("生成时间-距离图...")
            return self.params
        except Exception as e:
            self.logger.error(f"分析时间-距离关系失败: {str(e)}")
            raise
    
    def filter_outliers(self, b1_P=None, b2_P=None, b1_S=None, b2_S=None):
        """根据时间-距离关系过滤离群值"""
        self.logger.info("开始过滤离群值")
        
        # 更新边界参数
        if b1_P is not None: 
            self.params['b1_P'] = b1_P
            self.logger.info(f"更新参数: b1_P = {b1_P}")
        if b2_P is not None: 
            self.params['b2_P'] = b2_P
            self.logger.info(f"更新参数: b2_P = {b2_P}")
        if b1_S is not None: 
            self.params['b1_S'] = b1_S
            self.logger.info(f"更新参数: b1_S = {b1_S}")
        if b2_S is not None: 
            self.params['b2_S'] = b2_S
            self.logger.info(f"更新参数: b2_S = {b2_S}")
        
        if self.params['slope_P'] is None or self.params['slope_S'] is None:
            self.logger.info("拟合参数尚未计算，正在分析时间-距离关系...")
            self.analyze_time_distance()
        
        try:
            # 过滤条件
            def is_valid_p(row):
                if row['phase_type'] != 'P':
                    return False
                expected = self.params['slope_P'] * row['distance'] + self.params['b_P']
                upper_bound = expected + self.params['b1_P']
                lower_bound = expected - self.params['b2_P']
                return (row['travel_time'] <= upper_bound) and (row['travel_time'] >= lower_bound)
            
            def is_valid_s(row):
                if row['phase_type'] != 'S':
                    return False
                expected = self.params['slope_S'] * row['distance'] + self.params['b_S']
                upper_bound = expected + self.params['b1_S']
                lower_bound = expected - self.params['b2_S']
                return (row['travel_time'] <= upper_bound) and (row['travel_time'] >= lower_bound)
            
            # 应用过滤器
            self.logger.info("应用过滤条件...")
            tqdm.pandas(desc="过滤P波")
            self.phases_with_dist['valid_P'] = self.phases_with_dist.progress_apply(is_valid_p, axis=1)
            
            tqdm.pandas(desc="过滤S波")
            self.phases_with_dist['valid_S'] = self.phases_with_dist.progress_apply(is_valid_s, axis=1)
            
            self.phases_with_dist['valid'] = self.phases_with_dist['valid_P'] | self.phases_with_dist['valid_S']
            
            # 保存过滤后的数据
            output_file = self.output_dir / 'phase.dat_selection'
            self.logger.info(f"保存过滤后的数据到 {output_file}...")
            
            # 获取过滤后的数据
            filtered_data = self.phases_with_dist[self.phases_with_dist['valid']]
            
            # 对事件ID进行分组，提前准备好所有事件信息
            event_info = {row['event_id']: row['line'] for _, row in self.events.iterrows()}
            
            # 按事件ID排序，这样相同事件的数据会连续出现
            filtered_data = filtered_data.sort_values('event_id')
            
            # 批量写入文件
            with open(output_file, 'w') as f:
                current_event = None
                lines = []
                
                # 使用更高效的迭代方式
                for event_id, group in filtered_data.groupby('event_id'):
                    # 写入事件信息
                    lines.append(event_info[event_id] + '\n')
                    # 写入该事件的所有相位信息
                    lines.extend(row['line'] + '\n' for _, row in group.iterrows())
                
                # 批量写入所有行
                f.writelines(lines)
            
            # 统计结果
            n_p = self.phases_with_dist['valid_P'].sum()
            n_s = self.phases_with_dist['valid_S'].sum()
            total_p = (self.phases_with_dist['phase_type'] == 'P').sum()
            total_s = (self.phases_with_dist['phase_type'] == 'S').sum()
            
            p_removed = total_p - n_p
            s_removed = total_s - n_s
            
            stats_str = (
                f"过滤前P波数量: {total_p}\n"
                f"过滤后P波数量: {n_p}\n"
                f"删除的P波数量: {p_removed} ({p_removed/total_p*100:.1f}%)\n"
                f"过滤前S波数量: {total_s}\n"
                f"过滤后S波数量: {n_s}\n"
                f"删除的S波数量: {s_removed} ({s_removed/total_s*100:.1f}%)"
            )
            
            print(stats_str)
            self.logger.info(stats_str.replace('\n', ', '))
            
            # 保存统计信息
            with open(self.output_dir / 'filtering_stats.txt', 'w') as f:
                f.write(stats_str)
            
            return filtered_data
        
        except Exception as e:
            self.logger.error(f"过滤离群值失败: {str(e)}")
            raise
    
    def run_workflow(self, b1_P=None, b2_P=None, b1_S=None, b2_S=None):
        """运行完整工作流程"""
        self.logger.info("开始运行完整工作流程")
        start_time = time.time()
        
        try:
            self.load_station_data()
            self.load_phase_data()
            self.calculate_distance()
            self.analyze_time_distance()
            filtered_data = self.filter_outliers(b1_P, b2_P, b1_S, b2_S)
            
            end_time = time.time()
            elapsed_time = end_time - start_time
            self.logger.info(f"工作流程成功完成! 耗时: {elapsed_time:.2f}秒")
            
            return filtered_data
            
        except Exception as e:
            self.logger.error(f"工作流程执行失败: {str(e)}")
            raise

    def plot_time_distance(self, boundary_params=None):
        """绘制时间-距离关系图"""
        self.logger.info("开始绘制时间-距离关系图")
        
        if boundary_params is None:
            boundary_params = {
                'b1_P': 1.5, 'b2_P': 1.5,
                'b1_S': 1.5, 'b2_S': 1.5
            }
        
        if self.time_dist_data is None:
            self.logger.info("时间-距离数据未计算，正在计算...")
            self.calculate_distance()
        
        try:
            # 设置matplotlib参数
            plt.rcParams.update({
                # 使用sans-serif字体族，它会自动选择系统中最接近Helvetica的无衬线字体
                'font.family': 'sans-serif',
                'font.sans-serif': ['DejaVu Sans', 'Helvetica', 'Arial', 'sans-serif'],
                'mathtext.fontset': 'dejavusans',  # 使用DejaVu Sans作为数学字体
                'font.size': 8,
                'axes.labelsize': 8,
                'axes.titlesize': 9,
                'xtick.labelsize': 7,
                'ytick.labelsize': 7,
                'legend.fontsize': 7,
                'figure.dpi': 300,
                'figure.figsize': (7.2, 3.6),  # Nature single column width
                'figure.subplot.wspace': 0.3,
                'lines.linewidth': 1,
                'axes.linewidth': 0.5,
                'grid.linewidth': 0.5,
                'xtick.major.width': 0.5,
                'ytick.major.width': 0.5,
                'xtick.minor.width': 0.5,
                'ytick.minor.width': 0.5,
                'axes.grid': True,
                'grid.alpha': 0.3,
                'grid.linestyle': '--',
                'axes.axisbelow': True,
                'axes.facecolor': 'white',
                'figure.facecolor': 'white',
                'figure.edgecolor': 'white',
                'savefig.facecolor': 'white',
                'savefig.edgecolor': 'white',
                'savefig.bbox': 'tight',
                'savefig.pad_inches': 0.1,
                'pdf.fonttype': 42,  # 使用TrueType字体而不是Type 3
                'ps.fonttype': 42
            })

            # 分离P波和S波数据
            p_data = self.time_dist_data[self.time_dist_data['phase_type'] == 1][['travel_time', 'distance']]
            s_data = self.time_dist_data[self.time_dist_data['phase_type'] == 2][['travel_time', 'distance']]
            
            self.logger.info(f"P波数据点: {len(p_data)}个, S波数据点: {len(s_data)}个")
            
            # 清除之前的图形
            plt.close('all')
            
            # 创建图形
            fig, (ax1, ax2) = plt.subplots(1, 2)
 
             # 添加散点图
            ax1.scatter(p_data['distance'], p_data['travel_time'], 
                       c='blue', s=0.5, alpha=1)  # 更小的点，更透明           
            # P波图 - 使用hexbin显示数据密度
            hb1 = ax1.hexbin(p_data['distance'], p_data['travel_time'], 
                           gridsize=30,  # 调整网格大小
                           cmap='Blues',  # 使用蓝色色谱
                           mincnt=1,     # 最小计数
                           bins='log',   # 使用对数刻度
                           alpha=0.6)    # 增加密度图的不透明度
            
            legend_elements = []
            
            if self.params['slope_P'] is not None:
                x_max = p_data['distance'].max() if len(p_data) > 0 else 800
                x_p = np.linspace(0, x_max * 1.1)
                y_p = self.params['slope_P'] * x_p + self.params['b_P']
                line_p = ax1.plot(x_p, y_p, color='#000000', linestyle='-', label='Fitting line', zorder=3)[0]
                legend_elements.append(line_p)
                
                # 上下边界
                bound_up = ax1.plot(x_p, y_p + boundary_params['b1_P'], color='#4DBBD5', linestyle='--', 
                        label=f'+{boundary_params["b1_P"]} s', zorder=3)[0]
                bound_down = ax1.plot(x_p, y_p - boundary_params['b2_P'], color='#4DBBD5', linestyle='--', 
                        label=f'-{boundary_params["b2_P"]} s', zorder=3)[0]
                legend_elements.extend([bound_up, bound_down])
            
            # 使用自定义图例元素
            if legend_elements:
                ax1.legend(handles=legend_elements, frameon=True, edgecolor='none', fancybox=False)
            
            # 添加颜色条
            cb1 = plt.colorbar(hb1, ax=ax1)
            cb1.set_label('Count (log)', fontsize=7)
            cb1.ax.tick_params(labelsize=6)
            
            # 设置自适应的轴范围
            x_margin = p_data['distance'].std() * 0.1 if len(p_data) > 0 else 80
            y_margin = p_data['travel_time'].std() * 0.1 if len(p_data) > 0 else 20
            ax1.set_xlim(max(0, p_data['distance'].min() - x_margin), 
                        p_data['distance'].max() + x_margin if len(p_data) > 0 else 800)
            ax1.set_ylim(max(0, p_data['travel_time'].min() - y_margin), 
                        p_data['travel_time'].max() + y_margin if len(p_data) > 0 else 200)
            
            ax1.set_title('Time-Distance Curve (P-wave)')
            ax1.set_xlabel('Hypocentral Distance (km)')
            ax1.set_ylabel('Travel Time (s)')
            
            
            # 添加散点图
            ax2.scatter(s_data['distance'], s_data['travel_time'], 
                       c='red', s=0.5, alpha=1)  # 更小的点，更透明
             # S波图 - 使用hexbin显示数据密度
            hb2 = ax2.hexbin(s_data['distance'], s_data['travel_time'], 
                           gridsize=30,  # 调整网格大小
                           cmap='Reds',  # 使用红色色谱
                           mincnt=1,     # 最小计数
                           bins='log',   # 使用对数刻度
                           alpha=0.6)    # 增加密度图的不透明度           
            legend_elements = []
            
            if self.params['slope_S'] is not None:
                x_max = s_data['distance'].max() if len(s_data) > 0 else 800
                x_s = np.linspace(0, x_max * 1.1)
                y_s = self.params['slope_S'] * x_s + self.params['b_S']
                line_s = ax2.plot(x_s, y_s, color='#000000', linestyle='-', label='Fitting line', zorder=3)[0]
                legend_elements.append(line_s)
                
                # 上下边界
                bound_up = ax2.plot(x_s, y_s + boundary_params['b1_S'], color='#4DBBD5', linestyle='--', 
                        label=f'+{boundary_params["b1_S"]} s', zorder=3)[0]
                bound_down = ax2.plot(x_s, y_s - boundary_params['b2_S'], color='#4DBBD5', linestyle='--', 
                        label=f'-{boundary_params["b2_S"]} s', zorder=3)[0]
                legend_elements.extend([bound_up, bound_down])
            
            # 使用自定义图例元素
            if legend_elements:
                ax2.legend(handles=legend_elements, frameon=True, edgecolor='none', fancybox=False)
            
            # 添加颜色条
            cb2 = plt.colorbar(hb2, ax=ax2)
            cb2.set_label('Count (log)', fontsize=7)
            cb2.ax.tick_params(labelsize=6)
            
            # 设置自适应的轴范围
            x_margin = s_data['distance'].std() * 0.1 if len(s_data) > 0 else 80
            y_margin = s_data['travel_time'].std() * 0.1 if len(s_data) > 0 else 20
            ax2.set_xlim(max(0, s_data['distance'].min() - x_margin), 
                        s_data['distance'].max() + x_margin if len(s_data) > 0 else 800)
            ax2.set_ylim(max(0, s_data['travel_time'].min() - y_margin), 
                        s_data['travel_time'].max() + y_margin if len(s_data) > 0 else 200)
            
            ax2.set_title('Time-Distance Curve (S-wave)')
            ax2.set_xlabel('Hypocentral Distance (km)')
            ax2.legend(frameon=True, edgecolor='none', fancybox=False)
            
            # 调整布局
            plt.tight_layout()
            
            # # 保存图片
            # plot_file = self.output_dir / 't_dist.pdf'  # 使用PDF格式以保持矢量图质量
            # plt.savefig(plot_file)
            # self.logger.info(f"时间-距离图已保存到 {plot_file}")
            
            # 同时保存PNG格式用于快速预览
            plot_file_png = self.output_dir / 't_dist.png'
            plt.savefig(plot_file_png, dpi=300)
            
            return fig
            
        except Exception as e:
            self.logger.error(f"绘制时间-距离关系图失败: {str(e)}")
            raise

    def set_fitting_params(self, slope_P, b_P, slope_S, b_S, b1_P=None, b2_P=None, b1_S=None, b2_S=None):
        """手动设置拟合参数并更新图形"""
        self.logger.info("设置拟合参数")
        
        # 更新参数
        self.params['slope_P'] = slope_P
        self.params['b_P'] = b_P
        self.params['slope_S'] = slope_S
        self.params['b_S'] = b_S
        
        # 记录参数
        params_str = (
            f"P波: 斜率和截距 (slope_P, b_P) 为 {self.params['slope_P']:.6f} 和 {self.params['b_P']:.6f}\n"
            f"S波: 斜率和截距 (slope_S, b_S) 为 {self.params['slope_S']:.6f} 和 {self.params['b_S']:.6f}"
        )
        
        print(params_str)
        self.logger.info("已更新拟合参数")
        
        # 保存参数到文件
        with open(self.output_dir / 'fitting_params.txt', 'w') as f:
            f.write(params_str)
        
        return self.params


def main():
    """主函数"""
    # 设置参数
    station_file = './input/station.dat'
    phase_file = './input/phase.dat'
    output_dir = './output/'
    
    # 是否从缓存加载时间-距离数据
    load_from_cache = False  # 设置为True则从缓存加载，False则重新计算
    
    # 创建数据处理器
    processor = SeismicDataProcessor(
        station_file=station_file,
        phase_file=phase_file,
        output_dir=output_dir,
        load_from_cache=load_from_cache
    )
    
    try:
        # 设置拟合方式: 'auto' 为自动线性拟合, 'manual' 为手动设置参数
        fitting_method = 'auto'  # 在这里修改选择: 'auto' 或 'manual'
        
        # 设置边界参数（两种模式都会用到）
        boundary_params = {
            'b1_P': 4,  # P波上边界
            'b2_P': 4,  # P波下边界
            'b1_S': 5,  # S波上边界
            'b2_S': 5   # S波下边界
        }
        
        if fitting_method == 'auto':
            # 使用自动线性拟合
            processor.logger.info("使用自动线性拟合...")
            processor.analyze_time_distance()
            # 保存拟合参数
            processor.set_fitting_params(
                slope_P=processor.params['slope_P'],
                b_P=processor.params['b_P'],
                slope_S=processor.params['slope_S'],
                b_S=processor.params['b_S']
            )
            # 绘制图形
            processor.plot_time_distance(boundary_params)
            
        elif fitting_method == 'manual':
            # 先绘制原始数据图
            processor.logger.info("使用手动参数设置...")
            processor.calculate_distance()
            
            # 手动设置参数
            processor.set_fitting_params(
                slope_P=0.15,  # P波斜率
                b_P=0.3,      # P波截距
                slope_S=0.28,  # S波斜率
                b_S=0.5      # S波截距
            )
            # 绘制图形
            processor.plot_time_distance(boundary_params)
        else:
            processor.logger.warning("无效的拟合方式设置！将使用自动线性拟合...")
            processor.analyze_time_distance()
            # 绘制图形
            processor.plot_time_distance(boundary_params)
        
        # 进行数据过滤
        filtered_data = processor.filter_outliers(**boundary_params)
        processor.logger.info(f"处理完成! 过滤后数据已保存到 {Path(output_dir) / 'phase.dat_selection'}")
        
        return filtered_data
        
    except Exception as e:
        processor.logger.error(f"处理过程中发生错误: {str(e)}")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
