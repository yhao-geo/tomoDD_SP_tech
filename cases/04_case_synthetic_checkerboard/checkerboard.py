#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
棋盘格测试控制脚本 (Checkerboard Test Control Script)
==============================================

本程序用于tomoDD的棋盘格测试，主要功能包括：
1. 生成棋盘格速度模型
2. 生成合成数据
3. 反演恢复棋盘格

功能特点：
---------
1. 支持P波和S波速度异常
2. 自动生成合成数据
3. 自动进行反演恢复
4. 完整的错误处理和日志记录

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
import shutil
import numpy as np
from pathlib import Path
import subprocess
import logging

# 添加tools路径
sys.path.append('../../../../tools/')
import ModelingTools.MODTools as modt

class CheckerboardTest:
    def __init__(self):
        """初始化棋盘格测试参数"""
        # 程序设置
        self.tomoDD = "tomoDD_SP"
        self.tomoDD_syn = "tomoDD_SP_syn"
        self.tomoDD_inp = "tomoDD_SP.inp"
        self.tomoDD_syn_inp = "tomoDD_SP_syn.inp"
        
        # 棋盘格参数
        self.block_x = 1  # x方向块大小
        self.block_y = 1  # y方向块大小
        self.block_z = 1  # z方向块大小
        self.vp_perturbation = 0.05  # P波速度异常振幅 (5%)
        self.vs_perturbation = -0.05  # S波速度异常振幅 (-5%)
        
        # 目录和文件设置
        self.syn_dir = Path("Syn")
        self.vel_dir = Path("Vel")
        self.mod_file = Path("MOD")  # 原始MOD文件
        self.vp_model = "Vp_model.dat"
        self.vs_model = "Vs_model.dat"
        
        # 设置日志
        self.setup_logging()
        
        # 读取原始MOD文件
        self.logger.info("读取原始MOD文件")
        n, nx, ny, nz, X, Y, Z, vp, vs = modt.return_mod(str(self.mod_file))
        
        # 保存完整的网格信息
        self.n = n
        self.nx = nx
        self.ny = ny
        self.nz = nz
        self.X = X
        self.Y = Y
        self.Z = Z
        self.vp = vp
        self.vs = vs
        self.logger.info(f"模型尺寸: {nx}x{ny}x{nz}")
    
    def setup_logging(self):
        """设置日志系统"""
        self.logger = logging.getLogger("CheckerboardTest")
        self.logger.setLevel(logging.INFO)
        
        # 创建控制台处理器
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        
        # 创建格式化器
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        console_handler.setFormatter(formatter)
        
        # 添加处理器到日志记录器
        self.logger.addHandler(console_handler)
    
    def create_directories(self):
        """创建必要的目录"""
        self.syn_dir.mkdir(exist_ok=True)
        self.vel_dir.mkdir(exist_ok=True)
        shutil.copy("MOD", self.vel_dir / "MOD")
    
    def read_mod_file(self, mod_file: Path):
        """读取MOD文件
        
        Args:
            mod_file: MOD文件路径
        """
        try:
            # 使用MODTools读取模型
            n, nx, ny, nz, X, Y, Z, vp, vs = modt.return_mod(str(mod_file))
            
            # 保存完整的网格信息
            self.n=n
            self.nx = nx
            self.ny = ny
            self.nz = nz
            self.X = X
            self.Y = Y
            self.Z = Z
            
            return vp, vs
                
        except Exception as e:
            self.logger.error(f"读取MOD文件时出错: {str(e)}")
            raise
    
    def create_checkerboard(self, vel_data: np.ndarray, perturbation: float) -> np.ndarray:
        """生成棋盘格速度模型
        
        Args:
            vel_data: 原始速度模型
            perturbation: 速度异常振幅
            
        Returns:
            np.ndarray: 棋盘格速度模型
        """
        # 复制原始数据
        vel_checker = vel_data.copy()
        
        # 只对内部区域进行棋盘格扰动
        vel_checker[1:-1, 1:-1, 1:-1] = modt.checkerboard(
            vel_data[1:-1, 1:-1, 1:-1].copy(),
            self.block_z,
            self.block_y,
            self.block_x,
            perturbation
        )
        
        return vel_checker
    
    def write_mod_file(self, output_file: Path, vp_data: np.ndarray, vs_data: np.ndarray):
        """写入MOD文件
        
        Args:
            output_file: 输出文件路径
            vp_data: P波速度模型
            vs_data: S波速度模型
        """
        try:
            # 使用MODTools的write_mod函数写入文件
            modt.write_mod(
                str(output_file),
                self.n,  # n 参数，通常为1
                self.nx,
                self.ny,
                self.nz,
                self.X,
                self.Y,
                self.Z,
                vp_data,
                vs_data
            )
                    
        except Exception as e:
            self.logger.error(f"写入MOD文件时出错: {str(e)}")
            raise
    
    def dc2dc(self):
        """处理交叉相关走时数据
        将dtcc.syn转换为syn.dt.cc格式
        """
        try:
            if not os.path.exists("dtcc.syn"):
                self.logger.info("未找到dtcc.syn文件，跳过cc数据处理")
                return
                
            self.logger.info("开始处理cc数据...")
            with open("dtcc.syn", 'r') as f_in, open("syn.dt.cc", 'w') as f_out:
                ev1_old = None
                ev2_old = None
                for line in f_in:
                    parts = line.strip().split()
                    ev1, ev2, sta, t1, pha = parts
                    
                    if ev1_old is None:
                        ev1_old = ev1
                        ev2_old = ev2
                        f_out.write(f"# {ev1} {ev2} 0.0\n")
                    
                    if ev1 != ev1_old or ev2 != ev2_old:
                        ev1_old = ev1
                        ev2_old = ev2
                        f_out.write(f"# {ev1} {ev2} 0.0\n")
                    
                    f_out.write(f"{sta} {t1} {pha}\n")
                    
            self.logger.info("cc数据处理完成")
            
        except Exception as e:
            self.logger.error(f"处理cc数据时出错: {str(e)}")
            raise

    def process_synthetic_data(self):
        """处理合成数据"""
        # 复制必要的文件
        files_to_copy = [
            self.tomoDD_syn_inp,
            "ak135.15.SKS",
            "layer-16.dat",
            "extract_absolute_SP.awk",
            "extract_dtct_SP.awk",
            "extract_dtcc_SP.awk",  # 添加cc数据处理脚本
            self.tomoDD_syn
        ]
        
        for file in files_to_copy:
            shutil.copy(file, self.syn_dir)
        
        # 切换到合成数据目录
        os.chdir(self.syn_dir)
        
        # 设置执行权限并运行程序
        os.chmod(self.tomoDD_syn, 0o777)
        subprocess.run([f"./{self.tomoDD_syn}", self.tomoDD_syn_inp])
        
        # 处理数据文件
        self.abs2abs()
        self.dt2dt()
        self.dc2dc()  # 添加cc数据处理
        
        # 提取数据
        subprocess.run(["awk", "-f", "extract_absolute_SP.awk", "syn.absolute.dat"])
        subprocess.run(["awk", "-f", "extract_dtct_SP.awk", "syn.dt.ct"])
        
        # 如果存在cc数据，则处理
        if os.path.exists("syn.dt.cc"):
            self.logger.info("检测到syn.dt.cc文件，开始提取S-P数据...")
            subprocess.run(["awk", "-f", "extract_dtcc_SP.awk", "syn.dt.cc"])
            self.logger.info("S-P数据提取完成")
        
        os.chdir("..")
    
    def process_inversion(self):
        """进行反演"""
        # 复制必要的文件
        files_to_copy = [
            self.tomoDD_inp,
            "ak135.15.SKS",
            "layer-16.dat",
            self.tomoDD
        ]
        
        for file in files_to_copy:
            shutil.copy(file, self.vel_dir)
        
        # 切换到反演目录
        os.chdir(self.vel_dir)
        
        # 设置执行权限并运行程序
        os.chmod(self.tomoDD, 0o777)
        subprocess.run([f"./{self.tomoDD}", self.tomoDD_inp])
        
        os.chdir("..")
    
    def abs2abs(self):
        """处理绝对走时数据"""
        with open("absolute.syn", 'r') as f_in, open("syn.absolute.dat", 'w') as f_out:
            event_old = None
            for line in f_in:
                event, sta, time, weight, phase = line.strip().split()
                if event_old is None:
                    event_old = event
                    f_out.write(f"#  {event}\n")
                
                if event != event_old:
                    f_out.write(f"#  {event}\n")
                    event_old = event
                
                f_out.write(f"{sta}  {time}  {weight}  {phase}\n")
    
    def dt2dt(self):
        """处理相对走时数据"""
        with open("dt.syn", 'r') as f_in, open("syn.dt.ct", 'w') as f_out:
            ev1_old = None
            ev2_old = None
            for line in f_in:
                ev1, ev2, sta, t1, t2, qual, pha = line.strip().split()
                if ev1_old is None:
                    ev1_old = ev1
                    ev2_old = ev2
                    f_out.write(f"# {ev1} {ev2}\n")
                
                if ev1 != ev1_old or ev2 != ev2_old:
                    ev1_old = ev1
                    ev2_old = ev2
                    f_out.write(f"# {ev1} {ev2}\n")
                
                f_out.write(f"{sta} {t1} {t2} {qual} {pha}\n")
    
    def run(self):
        """运行棋盘格测试"""
        try:
            self.logger.info("开始棋盘格测试")
            
            # 创建目录
            self.create_directories()
            
            # 创建棋盘格模型
            vp_checker = self.create_checkerboard(self.vp, self.vp_perturbation)
            vs_checker = self.create_checkerboard(self.vs, self.vs_perturbation)
            self.logger.info("成功生成棋盘格模型")
            
            # 写入合成数据目录的MOD文件
            self.write_mod_file(self.syn_dir / "MOD", vp_checker, vs_checker)
            
            # 生成合成数据
            self.logger.info("开始生成合成数据")
            self.process_synthetic_data()
            self.logger.info("完成合成数据生成")
            
            # 进行反演
            self.logger.info("开始反演恢复")
            self.process_inversion()
            self.logger.info("完成反演恢复")
            
            self.logger.info("棋盘格测试完成")
            return 0
            
        except Exception as e:
            self.logger.error(f"棋盘格测试失败: {str(e)}")
            return 1

def main():
    """主函数"""
    checker = CheckerboardTest()
    return checker.run()

if __name__ == "__main__":
    sys.exit(main()) 