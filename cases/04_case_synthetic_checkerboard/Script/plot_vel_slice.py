#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
速度模型切片绘图程序 (Velocity Model Slice Plotting)
============================================

本程序用于绘制三维速度模型的水平切片棋盘测试结果，包括：
1. P波速度（Vp）切片
2. S波速度（Vs）切片
3. 速度比（Vp/Vs）切片

功能特点：
---------
1. 支持多深度层面的切片绘制
2. 自动插值到均匀网格
3. 支持线性和三次插值方法
4. 叠加地震事件分布
5. 自适应图形尺寸和配色方案
6. 输出高质量矢量图和位图

输入文件：
--------
- MOD：模型网格文件
- Vp_model.dat：P波速度模型
- Vs_model.dat：S波速度模型
- VpVs_model.dat：速度比模型
- tomoFDD.reloc：地震重定位结果

输出文件：
--------
- Vp_slice_*.pdf/png：P波速度切片图
- Vs_slice_*.pdf/png：S波速度切片图
- VpVs_slice_*.pdf/png：速度比切片图
（*表示深度值）

参数设置：
--------
- slice_dep：切片深度列表（km）
- event_dep_range_km：地震事件筛选范围（km）
- nx_uni, ny_uni：插值网格分辨率
- interp_method：插值方法（'linear'或'cubic'）

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
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import datetime
from scipy.stats import gaussian_kde
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
matplotlib.use('Agg')
from scipy.interpolate import interp2d  # 替换griddata
import os
import sys
sys.path.append('../../../tools/')
import ModelingTools.MODTools as modt

# ========== 参数设置 ==========
slice_dep = [5, 10.0, 20, 30.0]  # 需要绘制的深度（单位：km）
cmap_vp = 'jet_r'
cmap_vp = 'jet_r'
cmap_vpvs = 'jet'
outdir = './figure/vel_slices/'
os.makedirs(outdir, exist_ok=True)
event_dep_range_km = 5.0  # 地震事件筛选深度范围

mod_file = "../MOD"
modtrue_file = "../Syn/MOD"
vp_file = "../Vel/Vp_model.dat"
vs_file = "../Vel/Vs_model.dat"
vpvs_file = "../Vel/VpVs_model.dat"
event_file = "../Vel/tomoFDD.reloc"

n, nx, ny, nz, X_all, Y_all, Z_all, vp_ini, vs_ini = modt.return_mod(mod_file)
n_true, nx_true, ny_true, nz_true, X_all_true, Y_all_true, Z_all_true, vp_ini_true, vs_ini_true = modt.return_mod(modtrue_file)
vpvs_ini = vp_ini/vs_ini
X = X_all[1:-1]
Y = Y_all[1:-1]
Z = Z_all[1:-1]
vp_true = vp_ini_true[1:-1,1:-1,1:-1]
vs_true = vs_ini_true[1:-1,1:-1,1:-1]
vpvs_true = vp_true/vs_true

vp_data = 100*(np.loadtxt(vp_file).reshape(nz, ny, nx)[1:-1,1:-1,1:-1]-vp_ini[1:-1,1:-1,1:-1])/vp_ini[1:-1,1:-1,1:-1]
vs_data = 100*(np.loadtxt(vs_file).reshape(nz, ny, nx)[1:-1,1:-1,1:-1]-vs_ini[1:-1,1:-1,1:-1])/vs_ini[1:-1,1:-1,1:-1]
vpvs_data = 100*(np.loadtxt(vpvs_file).reshape(nz, ny, nx)[1:-1,1:-1,1:-1]-vpvs_ini[1:-1,1:-1,1:-1])/vpvs_ini[1:-1,1:-1,1:-1]

vp_true_data = 100*(vp_true-vp_ini[1:-1,1:-1,1:-1])/vp_ini[1:-1,1:-1,1:-1]
vs_true_data = 100*(vs_true-vs_ini[1:-1,1:-1,1:-1])/vs_ini[1:-1,1:-1,1:-1]
vpvs_true_data = 100*(vpvs_true-vpvs_ini[1:-1,1:-1,1:-1])/vpvs_ini[1:-1,1:-1,1:-1]

event_data = np.loadtxt(event_file)

# 地震事件信息
event_lat = event_data[:, 1]
event_lon = event_data[:, 2]
event_dep = event_data[:, 3]

# ========== 均匀网格 ==========
x_min, x_max = X.min(), X.max()
y_min, y_max = Y.min(), Y.max()
nx_uni, ny_uni = 10*5, 12*5  # 可调分辨率
x_uni = np.linspace(x_min, x_max, nx_uni)
y_uni = np.linspace(y_min, y_max, ny_uni)
X_uni, Y_uni = np.meshgrid(x_uni, y_uni)

# ========== 绘图函数 ==========
def plot_slice(data3d, Z, X, Y, slice_dep, X_uni, Y_uni, cmap, vlabel, out_prefix, interp_method='cubic',
               event_x=None, event_y=None, event_dep=None, event_dep_range_km=1.0, vmin=None, vmax=None):
    for dep in slice_dep:
        # 找到最接近的深度层
        iz = np.argmin(np.abs(Z - dep))
        dep_real = Z[iz]
        data_slice = data3d[iz, :, :]  # shape: (ny, nx)
        
        # 使用interp2d进行插值
        f = interp2d(X, Y, data_slice, kind=interp_method)
        grid_z = f(x_uni, y_uni)
        
        # 绘图
        plt.rcParams.update({
            'font.family': 'serif',
            'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif', 'serif'],
            'font.size': 12,
            'axes.labelsize': 14,
            'axes.titlesize': 16,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'axes.linewidth': 1.2,
            'pdf.fonttype': 42,
            'ps.fonttype': 42
        })
        
        # 自动设置图片宽高比与地理范围一致
        aspect_ratio = (x_max - x_min) / (y_max - y_min)
        base_height = 6  # 英寸
        fig_width = base_height * aspect_ratio
        fig, ax = plt.subplots(figsize=(fig_width, base_height))
        
        im = ax.imshow(grid_z, extent=[x_min, x_max, y_min, y_max], origin='lower',
                       cmap=cmap, interpolation='bicubic', aspect='equal', vmin=vmin, vmax=vmax)
        
        cbar = fig.colorbar(im, ax=ax, pad=0.05, fraction=0.05)
        cbar.set_label(vlabel, fontsize=13, fontweight='bold')
        ax.set_xlabel('X [km]', fontsize=14, fontweight='bold')
        ax.set_ylabel('Y [km]', fontsize=14, fontweight='bold')
        ax.set_ylim(y_max,y_min)
        ax.set_xlim(x_min,x_max)
        ax.set_title(f'{vlabel} Slice at {dep_real:.2f} km', fontsize=16, fontweight='bold')
        ax.grid(True, which='major', linestyle='--', linewidth=0.8, color='gray', alpha=0.3)
        
        # 叠加地震事件（学术风格优化）
        if (event_x is not None) and (event_y is not None) and (event_dep is not None):
            mask = np.abs(event_dep - dep_real) <= event_dep_range_km
            if np.any(mask):
                ax.scatter(event_x[mask], event_y[mask], s=3, c='k', alpha=0.7, marker='o',
                           edgecolors='w', linewidths=0.7, zorder=10, label='Events')
                ax.legend(fontsize=11, loc='best', frameon=True)
                
        ax.set_ylim(y_min,y_max)
        ax.set_xlim(x_min,x_max)    
        plt.tight_layout()
        out_png = f'{out_prefix}_{dep_real:.2f}km.png'
        out_pdf = f'{out_prefix}_{dep_real:.2f}km.pdf'
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.savefig(out_pdf, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"切片 {dep_real:.2f} km 已保存为 {out_png} / {out_pdf}")

# ========== 主流程 ==========
if __name__ == "__main__":
    for interp_method in ['cubic']:
        # 绘制反演结果
        plot_slice(vp_data, Z, X, Y, slice_dep, X_uni, Y_uni, cmap_vp, 'Vp (%)', f'{outdir}Vp_slice', interp_method=interp_method,
                   event_x=event_lon, event_y=event_lat, event_dep=event_dep, event_dep_range_km=event_dep_range_km, vmin=-5, vmax=5)
        plot_slice(vpvs_data, Z, X, Y, slice_dep, X_uni, Y_uni, cmap_vpvs, 'Vp/Vs (%)', f'{outdir}VpVs_slice', interp_method=interp_method,
                   event_x=event_lon, event_y=event_lat, event_dep=event_dep, event_dep_range_km=event_dep_range_km, vmin=-9.5, vmax=10.5)
        plot_slice(vs_data, Z, X, Y, slice_dep, X_uni, Y_uni, cmap_vpvs, 'Vs (%)', f'{outdir}Vs_slice', interp_method=interp_method,
                   event_x=event_lon, event_y=event_lat, event_dep=event_dep, event_dep_range_km=event_dep_range_km, vmin=-5, vmax=5)

        # 绘制真实模型
        plot_slice(vp_true_data, Z, X, Y, slice_dep, X_uni, Y_uni, cmap_vp, 'Vp (%)', f'{outdir}Vp_slice_true', interp_method=interp_method,
                   event_x=event_lon, event_y=event_lat, event_dep=event_dep, event_dep_range_km=event_dep_range_km, vmin=-5, vmax=5)
        plot_slice(vpvs_true_data, Z, X, Y, slice_dep, X_uni, Y_uni, cmap_vpvs, 'Vp/Vs (%)', f'{outdir}VpVs_slice_true', interp_method=interp_method,
                   event_x=event_lon, event_y=event_lat, event_dep=event_dep, event_dep_range_km=event_dep_range_km, vmin=-9.5, vmax=10.5)
        plot_slice(vs_true_data, Z, X, Y, slice_dep, X_uni, Y_uni, cmap_vpvs, 'Vs (%)', f'{outdir}Vs_slice_true', interp_method=interp_method,
                   event_x=event_lon, event_y=event_lat, event_dep=event_dep, event_dep_range_km=event_dep_range_km, vmin=-5, vmax=5)

