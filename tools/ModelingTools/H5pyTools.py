# 关于作者和文件的说明
# 文件: H5pyTools.py
# 作者: 杨浩
# 描述: 本文件包含对地震数据进行处理的函数，主要包括读取HDF5文件、去除离群点、以及根据样本数据生成模型等功能。

# 导入所需的库
import numpy as np
import h5py
from sklearn.ensemble import IsolationForest

def read_samples(fn_h5py, nv, n1=0, n2=0, filter_outliers=False, z_threshold=5):
    """
    读取指定的HDF5文件并处理样本数据。
    
    参数:
        fn_h5py: str，HDF5文件的路径。
        nv: int，模型的参数大小。
        n1: int，样本起始位置，默认为0。
        n2: int，样本结束位置，默认为0。
        filter_outliers: bool，是否过滤离群点，默认为False。
        z_threshold: float，Z-score阈值，默认为5，用于过滤离群点。
    
    返回:
        mean: np.ndarray，数据均值。
        std: np.ndarray，数据标准差。
        mean2: np.ndarray，比例数据均值。
        std2: np.ndarray，比例数据标准差。
        last: np.ndarray，最后一个样本数据。
    """
    
    # 读取 HDF5 文件
    with h5py.File(fn_h5py, 'r+') as f:
        samples = np.array(f['samples']).astype(float)

    # 选取指定范围的样本
    data = samples[n1:n2, :, :].reshape((-1, samples.shape[2]))
    ratio_data = np.divide(samples[n1:n2, :, :nv], samples[n1:n2, :, nv:2*nv], 
                          out=np.zeros_like(samples[n1:n2, :, :nv]), where=samples[n1:n2, :, nv:2*nv] != 0).reshape((-1, nv))

    # 过滤数据以去除离群点
    if filter_outliers:
        data = filter_outliers_zscore(data, z_threshold)
        ratio_data = filter_outliers_zscore(ratio_data, z_threshold)

    # 计算均值和标准差
    mean = np.mean(data, axis=0)
    std = np.std(data, axis=0)
    mean2 = np.mean(ratio_data, axis=0)
    std2 = np.std(ratio_data, axis=0)

    # 取最后一个样本
    last = samples[-1, :, :]

    return mean, std, mean2, std2, last

def return_distribution(fn_h5py, param_index, start=0, end=-1, filter_outliers=False, mad_threshold=2):
    """
    提取指定参数索引位置的所有样本数据，支持异常值过滤。

    参数:
        fn_h5py: str，HDF5文件路径。
        param_index: list[int]，要提取的参数索引列表。
        start: int，迭代次数的起始位置，默认从0开始。
        end: int，迭代次数的结束位置，默认到最后一项。
        filter_outliers: bool，是否过滤异常值，默认False。
        mad_threshold: float，MAD 阈值，仅当filter_outliers为True时生效，默认值为2。

    返回:
        samples_1d: np.ndarray，二维数组，形状为 [iterations * particles, len(param_index)]。
                    每列对应一个参数的数据，每行是所有迭代和粒子的样本数据。
    """
    with h5py.File(fn_h5py, 'r') as f:
        # 获取样本数据
        samples_data = f['samples'][start:end, :, :]
        
        # 创建一个空的列表来存储每个参数的展平数据
        all_samples = []
        
        # 遍历每个参数索引，提取并展平对应的样本数据
        for i in param_index:
            # 提取当前参数的数据，形状为 [iterations, particles]
            current_param_samples = samples_data[:, :, i]
            
            # 将当前参数的样本数据展平为一维数组，形状为 [iterations * particles]
            current_param_samples_1d = current_param_samples.reshape(-1)
            
            # 将展平后的数据添加到列表中
            all_samples.append(current_param_samples_1d)
        
        # 将所有参数的展平数据按列拼接成一个二维数组，形状为 [iterations * particles, len(param_index)]
        samples_1d = np.column_stack(all_samples)
        
        # 检查是否需要过滤异常值
        if filter_outliers:
            # 对每个参数的样本分别过滤异常值，返回长度不同的列表
            filtered_samples = [filter_outliers_mad(samples_1d[:, i], mad_threshold) for i in range(samples_1d.shape[1])]
            
            # 找到过滤后最长的列长度
            max_len = max(len(col) for col in filtered_samples)
            
            # 将所有过滤后的样本填充为同样长度，使用 np.nan 进行填充
            filtered_samples_padded = np.column_stack([np.pad(col, (0, max_len - len(col)), constant_values=np.nan) for col in filtered_samples])
            
            return filtered_samples_padded
        else:
            return samples_1d



def h5py2model(fn_h5py, n1, n2, nz, ny, nx):
    """
    根据样本数据生成模型。

    参数:
        fn_h5py: str，HDF5文件路径。
        n1: int，样本起始位置。
        n2: int，样本结束位置。
        nz, ny, nx: int，模型的维度。

    返回:
        vp_mean, vs_mean, event_mean, vp_std, vs_std, event_std: np.ndarray，分别表示Vp、Vs模型均值和标准差，以及事件模型均值和标准差。
    """
    
    mean, std, last = read_samples(fn_h5py, n1, n2)
    nv = nz * ny * nx
    events_num = int((len(mean) - 2 * nv) / 3)
    print('Event numbers : ', events_num)

    # 分离各个参数
    vp_mean = mean[0:nv].reshape((nz, ny, nx))
    vs_mean = mean[nv:2 * nv].reshape((nz, ny, nx))
    vp_std = std[0:nv].reshape((nz, ny, nx))
    vs_std = std[nv:2 * nv].reshape((nz, ny, nx))

    # 事件数据
    eventx_mean = mean[2 * nv:2 * nv + events_num]
    eventy_mean = mean[2 * nv + events_num:2 * nv + 2 * events_num]
    eventz_mean = mean[2 * nv + 2 * events_num:2 * nv + 3 * events_num]
    eventx_std = std[2 * nv:2 * nv + events_num]
    eventy_std = std[2 * nv + events_num:2 * nv + 2 * events_num]
    eventz_std = std[2 * nv + 2 * events_num:2 * nv + 3 * events_num]

    event_mean = np.column_stack((eventx_mean, eventy_mean, eventz_mean))
    event_std = np.column_stack((eventx_std, eventy_std, eventz_std))

    return vp_mean, vs_mean, event_mean, vp_std, vs_std, event_std

# 使用 Z-score 去除离群点
def filter_outliers_zscore(data, z_threshold=2):
    """
    使用 Z-score 方法过滤离群点。

    参数:
        data: np.ndarray，输入数据。
        z_threshold: float，Z-score 阈值。

    返回:
        np.ndarray，过滤后的数据。
    """
    # 计算均值和标准差
    mean = np.mean(data, axis=0)
    std_dev = np.std(data, axis=0)

    # 计算 Z-score，处理标准差为 0 的情况
    with np.errstate(divide='ignore', invalid='ignore'):
        z_scores = (data - mean) / std_dev

    # 处理多维数组的情况
    if data.ndim > 1:
        # 标准差为 0 的列的掩码
        std_dev_zero_mask = (std_dev == 0)
        z_scores[:, std_dev_zero_mask] = 0
        # 保留在每列上都满足 z_threshold 的行
        filtered_data = data[(np.abs(z_scores) < z_threshold).all(axis=1)]
    else:
        # 处理一维数组的情况
        if std_dev == 0:
            return data  # 如果所有值相同，不执行过滤
        filtered_data = data[np.abs(z_scores) < z_threshold]
    
    return filtered_data

def filter_outliers_mad(data, mad_threshold=2):
    """
    使用中位数绝对偏差（MAD）方法过滤离群点。

    参数:
        data: np.ndarray，输入数据。
        mad_threshold: float，MAD 阈值。

    返回:
        np.ndarray，过滤后的数据。
    """
    # 计算中位数和MAD
    median = np.median(data, axis=0)
    mad = np.median(np.abs(data - median), axis=0)

    # 计算每个数据点的偏差，处理 MAD 为 0 的情况
    with np.errstate(divide='ignore', invalid='ignore'):
        mad_scores = np.abs(data - median) / mad

    # 处理多维数组的情况
    if data.ndim > 1:
        # 标准差为 0 的列的掩码
        mad_zero_mask = (mad == 0)
        mad_scores[:, mad_zero_mask] = 0
        # 保留在每列上都满足 mad_threshold 的行
        filtered_data = data[(mad_scores < mad_threshold).all(axis=1)]
    else:
        # 处理一维数组的情况
        if mad == 0:
            return data  # 如果所有值相同，不执行过滤
        filtered_data = data[mad_scores < mad_threshold]
    
    return filtered_data



# 使用 IsolationForest 去除离群点
def filter_outliers_isolation_forest(data, contamination=0.1):
    """
    使用IsolationForest方法过滤离群点。
    
    参数:
        data: np.ndarray，输入数据。
        contamination: float，污染因子，默认值为0.1。
    
    返回:
        np.ndarray，过滤后的数据。
    """
    iso_forest = IsolationForest(contamination=contamination)
    outliers = iso_forest.fit_predict(data)
    return data[outliers == 1]


def h5py2model_sp(fn_h5py, n1, n2, nz, ny, nx):
    """
    根据样本数据生成带有Vp、Vs和Vp/Vs模型的扩展模型。

    参数:
        fn_h5py: str，HDF5文件路径。
        n1: int，样本起始位置。
        n2: int，样本结束位置。
        nz, ny, nx: int，模型的维度。

    返回:
        vp_mean, vs_mean, vpvs_mean, event_mean, vp_std, vs_std, vpvs_std, event_std: np.ndarray，分别表示Vp、Vs、Vp/Vs模型的均值和标准差，以及事件模型的均值和标准差。
    """
    
    mean, std, last = read_samples(fn_h5py, n1, n2)
    nv = nz * ny * nx
    events_num = int((len(mean) - 3 * nv) / 3)
    print('Event numbers : ', events_num)

    # 分离各个参数
    vp_mean = mean[0:nv].reshape((nz, ny, nx))
    vs_mean = mean[nv:2 * nv].reshape((nz, ny, nx))
    vpvs_mean = mean[2 * nv:3 * nv].reshape((nz, ny, nx))
    vp_std = std[0:nv].reshape((nz, ny, nx))
    vs_std = std[nv:2 * nv].reshape((nz, ny, nx))
    vpvs_std = std[2 * nv:3 * nv].reshape((nz, ny, nx))

    # 事件数据
    eventx_mean = mean[3 * nv:3 * nv + events_num]
    eventy_mean = mean[3 * nv + events_num:3 * nv + 2 * events_num]
    eventz_mean = mean[3 * nv + 2 * events_num:3 * nv + 3 * events_num]
    eventx_std = std[3 * nv:3 * nv + events_num]
    eventy_std = std[3 * nv + events_num:3 * nv + 2 * events_num]
    eventz_std = std[3 * nv + 2 * events_num:3 * nv + 3 * events_num]

    event_mean = np.column_stack((eventx_mean, eventy_mean, eventz_mean))
    event_std = np.column_stack((eventx_std, eventy_std, eventz_std))

    return vp_mean, vs_mean, vpvs_mean, event_mean, vp_std, vs_std, vpvs_std, event_std

def calculate_index(point, ny, nx):
    """
    根据给定点位置(point)，计算展平后的一维索引。

    参数:
        point: tuple/list，包含三个元素的点位置[x,y,z]。
        ny: int，数组的ny维度。
        nx: int，数组的nx维度。

    返回:
        para1: int，展平后的一维索引。
    """
    # 计算展平后的一维索引
    para_index = point[2] * ny * nx + point[1] * nx + point[0]
    return para_index
