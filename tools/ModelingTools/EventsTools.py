# 关于作者和文件的说明
# 文件: EventsTools.py
# 作者: 杨浩
# 描述: 本文件包含处理地震事件数据，主要包括替换数组行中的特定值等功能。

# 导入所需的库
import numpy as np

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

def filter_earthquakes(event, label, indice, threshold=0.15):
    """
    筛选地震数据，返回符合条件的掩码。

    参数：
    event : list or array
        地震数据的 y 坐标（例如，某一特定属性，如地震震源深度等）
    label : list or array
        与 event 对应的标签值，通常是地震相关的某些特征或位置坐标
    indice : int
        用于索引 label 中的特定值进行筛选
    threshold : float, optional
        筛选范围的阈值，默认为 0.15，表示以 label[indice] 为中心，上下浮动的最大值范围

    返回：
    masks : list
        布尔值列表，表示每个 event 是否在筛选范围内
    """
    # 获取 label 中对应 indice 索引的值作为筛选基准
    value = label[indice]
    
    # 根据阈值和筛选条件生成掩码
    masks = [(event[i] > (value - threshold)) and 
             (event[i] < (value + threshold)) 
             for i in range(len(event))]
    
    return masks

