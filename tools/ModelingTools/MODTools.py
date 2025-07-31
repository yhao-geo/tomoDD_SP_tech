# 关于作者和文件的说明
# 文件: MODTools.py
# 作者: 杨浩
# 描述: 本文件包含用于处理模型数据的多个函数，主要包括读取、写入、替换数组行、棋盘格扰动、以及数据平滑等功能。

# 导入所需的库
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.interpolate import RegularGridInterpolator
import numpy as np

def return_mod(filename):
    """
    从指定文件中读取模型数据。
    
    参数:
        filename: str，模型数据文件的路径。
    
    返回:
        n, nx, ny, nz, X, Y, Z, vp, vs：从文件中提取的各个参数和数组。
    """
    filename_MOD = filename
    with open(filename_MOD,'r') as file:
        n, nx, ny, nz = file.readline().split()
        nx = int(nx); ny = int(ny); nz = int(nz)
        X = np.array(file.readline().split()).astype(float)
        Y = np.array(file.readline().split()).astype(float)
        Z = np.array(file.readline().split()).astype(float)
        temp_p = np.zeros((nz, ny, nx))
        temp_ps = np.zeros((nz, ny, nx))
        for k in range(nz):
            for j in range(ny):
                temp_p[k,j,:] = file.readline().split()
        for k in range(nz):
            for j in range(ny):
                temp_ps[k,j,:] = file.readline().split()
    
    vp = np.array(temp_p).astype(float)
    vs = np.array(temp_p).astype(float) / np.array(temp_ps).astype(float)

    return n, nx, ny, nz, X, Y, Z, vp, vs

def write_mod(filename, n, nx, ny, nz, X, Y, Z, vp, vs):
    """
    将模型数据写入指定文件。
    
    参数:
        filename: str，目标文件的路径。
        n, nx, ny, nz: 模型的维度。
        X, Y, Z: 模型的坐标数组。
        vp, vs: 模型的速度数据。
    """
    with open(filename, 'w') as file:
        file.write(f"{n} {nx} {ny} {nz}\n")
        file.write(' '.join(map(str, X)) + '\n')
        file.write(' '.join(map(str, Y)) + '\n')
        file.write(' '.join(map(str, Z)) + '\n')
        for z in range(nz):
            for y in range(ny):
                for x in range(nx):
                    file.write(f'{vp[z, y, x]:.3f} ')
                file.write('\n')
        for z in range(nz):
            for y in range(ny):
                for x in range(nx):
                    file.write(f'{vp[z, y, x] / vs[z, y, x]:.3f} ')
                file.write('\n')

def checkerboard(v, nz, ny, nx, perturbation):
    """
    基于指定的块大小对三维numpy数组进行棋盘格模式的扰动。
    
    参数:
    v: 三维numpy数组，待扰动的数组。
    nz, ny, nx: 整数，表示z, y, x方向上的块大小。
    perturbation: 浮点数，扰动因子。
    
    返回:
    修改后的三维numpy数组。
    """
    # 计算扰动值
    value_1 = 1 - perturbation
    value_2 = 1 + perturbation
    
    # 获取数组的维度
    z_dim, y_dim, x_dim = v.shape
    
    # 遍历数组元素并应用棋盘格模式的扰动
    for k in range(z_dim):
        for j in range(y_dim):
            for i in range(x_dim):
                # 确定当前元素属于哪种类型的块
                if ((k // nz) + (j // ny) + (i // nx)) % 2 == 0:
                    v[k, j, i] *= value_1  # 应用第一个扰动值
                else:
                    v[k, j, i] *= value_2  # 应用第二个扰动值

    return v

def replace_rows(arr1, arr2):
    """
    根据第二个数组的第一列元素替换第一个数组中对应行的特定列。
    
    参数:
        arr1: np.ndarray，第一个数组，包含待替换行的原始数据。
        arr2: np.ndarray，第二个数组，包含用于替换的数据，第一列作为匹配依据。
    
    返回:
        np.ndarray，替换后的数组。
    """
    # 创建一个字典，以第二个数组的第一列元素为键，对应行为值
    dict_arr2 = {str(int(row[0])): row[1:4] for row in arr2}

    # 创建一个新的列表，用于存储需要保留的行
    retained_rows = []

    # 遍历第一个数组，替换或删除行
    for row in arr1:
        key = str(int(row[9]))  # 用作匹配的键
        if key in dict_arr2:
            # 如果匹配，将相关列替换
            row[2:5] = dict_arr2[key]
            retained_rows.append(row)  # 将匹配的行加入新的列表
        else:
            # 如果不匹配，不将其添加到 retained_rows 中
            continue

    # 将列表转换为 numpy 数组
    arr1 = np.array(retained_rows)

    return arr1



def interpolate_velocity_3d(initial_grid, target_grid, values, method='cubic'):
    """
    对初始网格的数据进行插值，生成目标网格的对应值。
    
    参数:
        initial_grid: tuple of np.ndarray，初始网格的坐标数组，分别对应各个维度（z, y, x）。
        values: np.ndarray，初始网格上的数据值。
        target_grid: tuple of np.ndarray，目标网格的坐标数组。
    
    返回:
        np.ndarray，目标网格上的插值结果。
    """
    # 创建插值器
    interpolator = RegularGridInterpolator(initial_grid, values, method=method, bounds_error=False, fill_value=None)
    
    # 对目标网格进行插值
    target_values = np.zeros((len(target_grid[0]), len(target_grid[1]), len(target_grid[2])))

    # 进行插值
    for k in range(len(target_grid[0])):
        for j in range(len(target_grid[1])):
            for i in range(len(target_grid[2])):
                target_values[k, j, i] = interpolator(np.array([target_grid[0][k], target_grid[1][j], target_grid[2][i]])).item()
    
    return target_values

def smooth(vel, initial_grid, level=1, method='cubic'):
    """
    对输入的速度数组进行插值到均匀网格后，再进行平滑处理，最后插值回初始网格。
    
    参数:
        vel: np.ndarray，待平滑的速度数组。
        initial_grid: tuple of np.ndarray，初始网格的坐标数组，分别对应各个维度（z, y, x）。
        level: 整数，平滑的程度，默认为1。
        method: 字符串，插值的方法，默认为 'cubic'。
    
    返回:
        np.ndarray，插值和平滑后的速度数组，大小与initial_grid一致。
    """
    # Step 1: 根据初始网格生成均匀网格
    z_grid = np.linspace(initial_grid[0][0], initial_grid[0][-1], num=initial_grid[0].size * 2)  # 目标网格 Z 轴
    y_grid = np.linspace(initial_grid[1][0], initial_grid[1][-1], num=initial_grid[1].size * 2)  # 目标网格 Y 轴
    x_grid = np.linspace(initial_grid[2][0], initial_grid[2][-1], num=initial_grid[2].size * 2)  # 目标网格 X 轴
    
    # Step 2: 将原始数据插值到均匀网格
    vel_interpolated = interpolate_velocity_3d(initial_grid, (z_grid, y_grid, x_grid), vel, method=method)
    
    # Step 3: 对插值后的数据进行平滑
    vel_smooth = gaussian_filter(vel_interpolated, level)
    
    # Step 4: 将平滑后的数据重新插值回初始网格
    vel_final = interpolate_velocity_3d((z_grid, y_grid, x_grid), initial_grid, vel_smooth, method=method)
    
    return vel_final
