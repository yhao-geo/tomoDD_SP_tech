import numpy as np
import matplotlib.pyplot as plt
import math

def mygrt(elt, eln, slt, sln):
    PI = 3.141592653
    
    # Convert degrees to radians
    slat = slt * PI / 180.0
    slon = sln * PI / 180.0
    elat = elt * PI / 180.0
    elon = eln * PI / 180.0

    # Correct for ellipticity
    slat = math.atan(0.996647 * math.tan(slat))
    elat = math.atan(0.996647 * math.tan(elat))

    # Go to colatitudes
    slat = PI / 2.0 - slat
    elat = PI / 2.0 - elat

    # Make all longitudes positive
    if slon < 0.0:
        slon = slon + 2.0 * PI
    if elon < 0.0:
        elon = elon + 2.0 * PI

    # Compute direction cosines
    a = math.sin(elat) * math.cos(elon)
    b = math.sin(elat) * math.sin(elon)
    c = math.cos(elat)
    a1 = math.sin(slat) * math.cos(slon)
    b1 = math.sin(slat) * math.sin(slon)
    c1 = math.cos(slat)

    cd = a * a1 + b * b1 + c * c1

    # Make sure acos won't fail
    cd = min(max(cd, -1.0), 1.0)

    info = {}

    # Calculate the delta (angle between two points)
    info['del'] = math.acos(cd) * 180.0 / PI
    info['dist'] = info['del'] * PI * 6371.0 / 180.0

    # Calculate azimuth (az)
    tmp1 = math.cos(elon) * math.cos(slon) + math.sin(elon) * math.sin(slon)
    tmp2a = 1.0 - cd * cd

    if tmp2a <= 0.0:
        tmp2 = 0.0
        tmp3 = 1.0
    else:
        tmp2 = math.sqrt(tmp2a)
        tmp3 = (math.sin(elat) * math.cos(slat) - math.cos(elat) * math.sin(slat) * tmp1) / tmp2

    # Make sure acos won't fail
    tmp3 = min(max(tmp3, -1.0), 1.0)

    z = math.acos(tmp3)

    # Adjust orientation for azimuth
    if (math.sin(slon) * math.cos(elon) - math.cos(slon) * math.sin(elon)) < 0.0:
        z = 2.0 * PI - z

    info['az'] = 180.0 * z / PI

    # Calculate back azimuth (baz)
    tmp1 = math.cos(slon) * math.cos(elon) + math.sin(slon) * math.sin(elon)
    tmp2a = 1.0 - cd * cd

    if tmp2a <= 0.0:
        tmp2 = 0.0
        tmp3 = 1.0
    else:
        tmp2 = math.sqrt(tmp2a)
        tmp3 = (math.sin(slat) * math.cos(elat) - math.cos(slat) * math.sin(elat) * tmp1) / tmp2

    # Make sure acos won't fail
    tmp3 = min(max(tmp3, -1.0), 1.0)

    bz = math.acos(tmp3)

    # Adjust orientation for back azimuth
    if (math.sin(elon) * math.cos(slon) - math.cos(elon) * math.sin(slon)) < 0.0:
        bz = 2.0 * PI - bz

    info['baz'] = 180.0 * bz / PI

    return info

def sph2car(dell, az, dep):
    s_fac = 180.0 / 3.1415926

    d1 = 6371.0
    d3 = d1 - dep
    d2 = math.sqrt(d1 * d1 + d3 * d3 - 2.0 * d1 * d3 * math.cos(dell / s_fac))

    z = 0.5 * (d1 * d1 + d2 * d2 - d3 * d3) / d1

    dist = math.sqrt(d2 * d2 - z * z)

    x = dist * math.sin(az / s_fac)
    y = dist * math.cos(az / s_fac)

    return x, y, z

def rotate(x, y, theta):
    xr = x * math.cos(theta) - y * math.sin(theta)
    yr = x * math.sin(theta) + y * math.cos(theta)
    return xr, yr


def sph2car_ft(lat_ft, lon_ft, dep_ft, wlat_ft, wlon_ft, theta_ft):
    # 转换成弧度
    lat = lat_ft
    lon = lon_ft
    wlat = wlat_ft
    wlon = wlon_ft
    theta = theta_ft

    # 计算球面坐标
    dat = mygrt(wlat, wlon, lat, lon)
    del_angle = dat['del']
    az_angle = dat['az']

    # 转换为笛卡尔坐标系
    x, y, z = sph2car(del_angle, az_angle, dep_ft)

    # 旋转
    xr, yr = rotate(x, y, theta)

    return xr, yr, z


if __name__=="__main__":

    wlat           =  35.97420383
    wlon           = -120.5521413
    rota           = (-0/180)*math.pi

    src_x,src_y,src_z = 35.97420383,-120.5521413,0
    src_x_e,src_y_e,src_z_e = sph2car_ft(src_x+1, src_y, src_z, wlat, wlon, rota)

    print(f"{src_x_e,src_y_e,src_z_e}")

