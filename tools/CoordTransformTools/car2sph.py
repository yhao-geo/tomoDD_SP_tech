import numpy as np
import math

# Function to convert Cartesian to spherical coordinates
def car2sph_ft(x_ft, y_ft, z_ft, wlat_ft, wlon_ft, theta_ft):
    wlat = wlat_ft
    wlon = wlon_ft
    x = x_ft
    y = y_ft
    z = z_ft
    theta = theta_ft
    
    theta = -theta
    xr, yr = rotate(x, y, theta)
    
    del_val, az, dep = car2sph(xr, yr, z)
    lat, lon = fpos(del_val, az, wlat, wlon)
    if lon > 180:
    	lon -= 360
    
    return lat, lon, dep

# Function to rotate points in (x, y) plane
def rotate(x, y, theta):
    xr = x * np.cos(theta) - y * np.sin(theta)
    yr = x * np.sin(theta) + y * np.cos(theta)
    return xr, yr

# Convert Cartesian coordinates (x, y, z) to spherical coordinates (distance, azimuth, elevation)
def car2sph(x, y, z):
    sfac = 180.0 / np.pi
    dx = x
    dy = y
    dz = z

    az = get_az(x, y)

    d1 = 6371.0
    d2 = np.sqrt(dx * dx + dy * dy + dz * dz)
    d3 = np.sqrt(dx * dx + dy * dy + (dz - 6371.0) * (dz - 6371.0))

    # Calculate elevation and distance in spherical coordinates
    del_val = sfac * np.arccos((d2 * d2 - d1 * d1 - d3 * d3) / (-2.0 * d1 * d3))
    dep = d1 - d3
    
    return del_val, az, dep

# Function to calculate azimuth
def get_az(x, y):
    dx = x
    dy = y
    sfac = 180.0 / np.pi

    if abs(dx) < 1.0e-10:
        if dy > 0.0:
            az = 0.0
            return az
        else:
            az = 180.0
            return az

    if dx >= 0.0:
        az = 90.0 - (sfac * np.arctan(dy / dx))
        return az
    else:
        az = 270.0 - (sfac * np.arctan(dy / dx))
        return az

# Calculate latitude and longitude from distance, azimuth, target's latitude, and longitude
def fpos(del_val, az, wlat, wlon):
    f1 = np.pi / 180.0

    lt, ln = locate(f1 * wlat, f1 * wlon, f1 * del_val, f1 * az)
    lat = lt / f1
    lon = ln / f1

    return lat, lon

# Function to locate on a spherical surface
def locate(oldlat, oldlon, angd, az):
    PI = np.pi
    tlat1, tlon1 = rot2(PI / 2.0 - angd, PI - az, PI / 2.0 - oldlat)
    newlat, newlon = rot3(tlat1, tlon1, oldlon)
    return newlat, newlon

# First rotation function
def rot2(oldlat, oldlon, angle):
    x1, x2, x3 = sphcar(oldlat, oldlon)
    temp = x1
    x1 = x1 * np.cos(angle) + x3 * np.sin(angle)
    x3 = x3 * np.cos(angle) - temp * np.sin(angle)
    newlat, newlon = carsph(x1, x2, x3)
    return newlat, newlon

# Second rotation function
def rot3(oldlat, oldlon, angle):
    pi = np.pi
    newlat = oldlat
    newlon = oldlon + angle

    if newlon > 2.0 * pi:
        newlon = newlon - 2.0 * pi
    elif newlon < -2.0 * pi:
        newlon = newlon + 2.0 * pi
    
    return newlat, newlon

def sphcar(lat, lon):
    x1 = np.cos(lat) * np.cos(lon)
    x2 = np.cos(lat) * np.sin(lon)
    x3 = np.sin(lat)
    return x1, x2, x3

# Convert (x1, x2, x3) coordinates in spherical coordinates to latitude and longitude in spherical coordinate system
def carsph(x1, x2, x3):
    PI = np.pi

    if x3 > 1.0:
        x3 = 1.0

    if x3 < -1.0:
        x3 = -1.0

    lat = np.arcsin(x3)

    if x1 == 0.0:
        if x2 > 0.0:
            lon = PI / 2.0
        if x2 < 0.0:
            lon = 3.0 * PI / 2.0
    else:
        lon = np.arctan(x2 / x1)

        if x1 < 0.0:
            lon = lon + PI

        if x1 > 0.0 and x2 < 0.0:
            lon = lon + 2.0 * PI

    return lat, lon


if __name__=="__main__":

    wlat           =  35.97420383
    wlon           = -120.5521413
    rota           = (0/180)*math.pi

    src_x,src_y,src_z = 0,0,12
    src_x_e,src_y_e,src_z_e = car2sph_ft(src_x, src_y, src_z, wlat, wlon, rota)

    print(f"{src_x_e,src_y_e,src_z_e}")
