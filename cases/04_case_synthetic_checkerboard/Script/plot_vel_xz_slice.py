#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
速度模型经度-深度剖面绘图程序 (Velocity Model X-Z Cross-section Plotting)
=======================================================

本程序用于绘制三维速度模型的经度-深度剖面棋盘测试结果，包括：
1. P波速度（Vp）剖面
2. S波速度（Vs）剖面
3. 速度比（Vp/Vs）剖面

功能特点：
---------
1. 支持多个纬度位置的剖面绘制
2. 自动插值到均匀网格
3. 支持线性和三次插值方法
4. 叠加地震事件投影
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
- Vp_xz_slice_*.pdf/png：P波速度剖面图
- Vs_xz_slice_*.pdf/png：S波速度剖面图
- VpVs_xz_slice_*.pdf/png：速度比剖面图
（*表示纬度值）

参数设置：
--------
- slice_lat：剖面纬度列表（度）
- event_lat_range：地震事件筛选范围（度）
- nx_uni, nz_uni：插值网格分辨率
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
import matplotlib
matplotlib.use('Agg')
from scipy.interpolate import interp2d
import os
import sys
sys.path.append('../../../tools/')
import ModelingTools.MODTools as modt

# ========== 参数设置 ==========
slice_lat = [28, 29, 30, 31.0]  # 需要绘制的纬度（单位：度）
cmap_vp = 'jet_r'
cmap_vs = 'jet_r'
cmap_vpvs = 'jet'
outdir = './figure/vel_xz_slices/'
os.makedirs(outdir, exist_ok=True)
event_lat_range = 0.5  # 地震事件筛选范围（度）

mod_file = "../MOD"
modtrue_file = "../Syn/MOD"
vp_file = "../Vel/Vp_model.dat"
vs_file = "../Vel/Vs_model.dat"
vpvs_file = "../Vel/VpVs_model.dat"
event_file = "../Vel/tomoFDD.reloc"

# ========== 读取数据 ==========
n, nx, ny, nz, X_all, Y_all, Z_all, vp_ini, vs_ini = modt.return_mod(mod_file)
n_true, nx_true, ny_true, nz_true, X_all_true, Y_all_true, Z_all_true, vp_ini_true, vs_ini_true = modt.return_mod(modtrue_file)
vpvs_ini = vp_ini/vs_ini
X = X_all[1:-1]
Y = Y_all[1:-1]
Z = Z_all[1:-4]
vp_true = vp_ini_true[1:-4,1:-1,1:-1]
vs_true = vs_ini_true[1:-4,1:-1,1:-1]
vpvs_true = vp_true/vs_true

# 计算速度异常百分比
vp_data = 100*(np.loadtxt(vp_file).reshape(nz, ny, nx)[1:-4,1:-1,1:-1]-vp_ini[1:-4,1:-1,1:-1])/vp_ini[1:-4,1:-1,1:-1]
vs_data = 100*(np.loadtxt(vs_file).reshape(nz, ny, nx)[1:-4,1:-1,1:-1]-vs_ini[1:-4,1:-1,1:-1])/vs_ini[1:-4,1:-1,1:-1]
vpvs_data = 100*(np.loadtxt(vpvs_file).reshape(nz, ny, nx)[1:-4,1:-1,1:-1]-vpvs_ini[1:-4,1:-1,1:-1])/vpvs_ini[1:-4,1:-1,1:-1]

# 计算真实模型异常百分比
vp_true_data = 100*(vp_true-vp_ini[1:-4,1:-1,1:-1])/vp_ini[1:-4,1:-1,1:-1]
vs_true_data = 100*(vs_true-vs_ini[1:-4,1:-1,1:-1])/vs_ini[1:-4,1:-1,1:-1]
vpvs_true_data = 100*(vpvs_true-vpvs_ini[1:-4,1:-1,1:-1])/vpvs_ini[1:-4,1:-1,1:-1]

event_data = np.loadtxt(event_file)

# 地震事件信息
event_lat = event_data[:, 1]
event_lon = event_data[:, 2]
event_dep = event_data[:, 3]

# ========== 均匀网格 ==========
x_min, x_max = X.min(), X.max()
z_min, z_max = Z.min(), Z.max()
nx_uni, nz_uni = 50, 30  # 可调分辨率
x_uni = np.linspace(x_min, x_max, nx_uni)
z_uni = np.linspace(z_min, z_max, nz_uni)
X_uni, Z_uni = np.meshgrid(x_uni, z_uni)

# ========== 绘图函数 ==========
def plot_xz_slice(data3d, Z, X, Y, slice_lat, X_uni, Z_uni, cmap, vlabel, out_prefix, interp_method='cubic',
                  event_x=None, event_y=None, event_z=None, event_lat_range=0.1, vmin=None, vmax=None):
    """绘制经度-深度剖面
    
    Args:
        data3d: 3D速度数据 (nz, ny, nx)
        Z, X, Y: 网格坐标
        slice_lat: 剖面纬度列表
        X_uni, Z_uni: 插值用的均匀网格
        cmap: 颜色方案
        vlabel: 速度标签
        out_prefix: 输出文件前缀
        interp_method: 插值方法
        event_x, event_y, event_z: 地震事件坐标
        event_lat_range: 地震事件筛选范围（度）
        vmin, vmax: 颜色范围
    """
    for lat in slice_lat:
        # 找到最接近的纬度层
        iy = np.argmin(np.abs(Y - lat))
        lat_real = Y[iy]
        data_slice = data3d[:, iy, :]  # shape: (nz, nx)
        
        # 使用interp2d进行插值
        f = interp2d(X, Z, data_slice, kind=interp_method)
        grid_v = f(x_uni, z_uni)
        
        # 绘图设置
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
        
        # 自动设置图片宽高比
        aspect_ratio = (x_max - x_min) / (z_max - z_min)*12
        base_height = 4  # 英寸
        fig_width = base_height * aspect_ratio
        fig, ax = plt.subplots(figsize=(fig_width, base_height))
        
        # 绘制速度场
        im = ax.imshow(grid_v, extent=[x_min, x_max, z_max, z_min], 
                      origin='upper', cmap=cmap, interpolation='bicubic', 
                      aspect='auto', vmin=vmin, vmax=vmax)
        
        # 添加颜色条
        cbar = fig.colorbar(im, ax=ax, pad=0.05, fraction=0.05)
        cbar.set_label(vlabel, fontsize=13, fontweight='bold')
        
        # 设置轴标签
        ax.set_xlabel('Longitude (km)', fontsize=14, fontweight='bold')
        ax.set_ylabel('Depth (km)', fontsize=14, fontweight='bold')
        ax.set_title(f'{vlabel} X-Z Section at {lat_real:.2f}°N', fontsize=16, fontweight='bold')
        
        # 添加网格
        ax.grid(True, which='major', linestyle='--', linewidth=0.8, color='gray', alpha=0.3)
        
        # 叠加地震事件
        if all(v is not None for v in [event_x, event_y, event_z]):
            mask = np.abs(event_y - lat_real) <= event_lat_range
            if np.any(mask):
                ax.scatter(event_x[mask], event_z[mask], s=3, c='k', alpha=0.7,
                          marker='o', edgecolors='w', linewidths=0.7, zorder=10,
                          label='Events')
                ax.legend(fontsize=11, loc='upper right', frameon=True)
        
        ax.set_ylim(z_max,z_min)
        ax.set_xlim(x_min,x_max)
        
        # 保存图片
        plt.tight_layout()
        out_png = f'{out_prefix}_{lat_real:.2f}N.png'
        out_pdf = f'{out_prefix}_{lat_real:.2f}N.pdf'
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.savefig(out_pdf, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"剖面 {lat_real:.2f}°N 已保存为 {out_png} / {out_pdf}")

# ========== 主流程 ==========
if __name__ == "__main__":
    for interp_method in ['cubic']:
        # 绘制反演结果
        plot_xz_slice(vp_data, Z, X, Y, slice_lat, X_uni, Z_uni, cmap_vp, 'Vp (%)', 
                     f'{outdir}Vp_xz_slice', interp_method=interp_method,
                     event_x=event_lon, event_y=event_lat, event_z=event_dep, 
                     event_lat_range=event_lat_range, vmin=-5, vmax=5)
        
        plot_xz_slice(vs_data, Z, X, Y, slice_lat, X_uni, Z_uni, cmap_vs, 'Vs (%)', 
                     f'{outdir}Vs_xz_slice', interp_method=interp_method,
                     event_x=event_lon, event_y=event_lat, event_z=event_dep, 
                     event_lat_range=event_lat_range, vmin=-5, vmax=5)
        
        plot_xz_slice(vpvs_data, Z, X, Y, slice_lat, X_uni, Z_uni, cmap_vpvs, 'Vp/Vs (%)', 
                     f'{outdir}VpVs_xz_slice', interp_method=interp_method,
                     event_x=event_lon, event_y=event_lat, event_z=event_dep, 
                     event_lat_range=event_lat_range, vmin=-9.5, vmax=10.5)
        
        # 绘制真实模型
        plot_xz_slice(vp_true_data, Z, X, Y, slice_lat, X_uni, Z_uni, cmap_vp, 'Vp (%)', 
                     f'{outdir}Vp_xz_slice_true', interp_method=interp_method,
                     event_x=event_lon, event_y=event_lat, event_z=event_dep, 
                     event_lat_range=event_lat_range, vmin=-5, vmax=5)
        
        plot_xz_slice(vs_true_data, Z, X, Y, slice_lat, X_uni, Z_uni, cmap_vs, 'Vs (%)', 
                     f'{outdir}Vs_xz_slice_true', interp_method=interp_method,
                     event_x=event_lon, event_y=event_lat, event_z=event_dep, 
                     event_lat_range=event_lat_range, vmin=-5, vmax=5)
        
        plot_xz_slice(vpvs_true_data, Z, X, Y, slice_lat, X_uni, Z_uni, cmap_vpvs, 'Vp/Vs (%)', 
                     f'{outdir}VpVs_xz_slice_true', interp_method=interp_method,
                     event_x=event_lon, event_y=event_lat, event_z=event_dep, 
                     event_lat_range=event_lat_range, vmin=-9.5, vmax=10.5) 
