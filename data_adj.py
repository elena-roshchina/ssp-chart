import copy
import math
import os
import re
from fractions import Fraction

from netCDF4._netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
import datetime as dt
from datetime import timedelta, time
import pandas as pd
from matplotlib import pyplot as plt



STD_HORIZONS = [0.0, 5.0, 10.0, 15.0, 20.0,  25.0, 30.0, 40.0, 50.0, 60.0, 75.0,
                100.0,  125.0,  150.0, 200.0, 250.0, 300.0, 350.0,
                400.0, 450.0, 500.0, 600.0, 800.0,  1000.0,  1200.0,  1500.0, 1800.0, 2000.0]


def density_w(t):
    # sea water density as function of temperature
    a_coeff = np.array([999.842594, 6.793952E-2, -9.095290E-3,
                        1.001685E-4, -1.120083E-6, 6.536332E-9], float)
    t_deg = np.array([t**x for x in range(6)], float)
    pw = np.dot(a_coeff, t_deg)
    return pw


def bulk_modulus(sal, t, press):
    #  K(S, t, p)
    t_degr = np.array([t ** x for x in range(5)], float)
    p_degr = np.array([t ** x for x in range(3)], float)

    #  K_w
    coeff_e = np.array([19652.21, 148.4206, -2.327105, 1.360477E-2, -5.155288E-5], float)
    coeff_f = np.array([54.6746, -0.603459, 1.09987E-2, -6.1670E-5, 0.00000000000], float) * sal
    coeff_g = np.array([7.944E-2, 1.6483E-2, -5.3009E-4, 0.00000000000, 0.000000000], float) * (sal ** Fraction(3,2))

    k_s_t_02 = np.dot(coeff_e + coeff_f + coeff_g, t_degr)

    #  А
    coeff_i = np.array([2.2838E-3, -1.0981E-5, -1.6078E-6, 0.0, 0.0], float) * sal
    coeff_h = np.array([3.239908, 1.43713E-3, 1.16092E-4, -5.77905E-7, 0.0], float)
    coeff_j = 1.91075E-4
    aa = np.dot(coeff_h + coeff_i, t_degr) + coeff_j * (sal ** Fraction(3,2))
    # B

    coeff_k = np.array([8.50935E-5, -6.12293E-6, 5.2787E-8, 0.0, 0.0], float)
    coeff_m = np.array([-9.9348E-7, 2.0816E-8, 9.1697E-10, 0.0, 0.0], float) * sal
    bb = np.dot(coeff_k + coeff_m, t_degr)

    #  K
    coefficients = np.array([k_s_t_02, aa, bb])

    return np.dot(coefficients, p_degr)


def density(sal, temp, pressure):
    # pressure unit = decibar, salinity = psu, temperatures in degree C
    # scale pressure to bars
    pressure /= 10
    t_degr = np.array([temp ** x for x in range(5)], float)
    b = np.array([8.24493E-1, -4.0899E-3, 7.6438E-5, -8.2467E-7, 5.3875E-9]) * sal
    c = np.array([-5.72466E-3, 1.0227E-4, -1.6546E-6, 0.0, 0.0]) * (sal ** Fraction(3,2))
    d0 = 4.8314E-4
    dens_s_t_0 = density_w(t=temp) + np.dot(b + c,t_degr) + d0 * (sal ** 2)
    return dens_s_t_0 / (1 - pressure / bulk_modulus(sal=sal, t=temp, press=pressure))


def vilson_sound_velocity(sal, temp, pressure):
    # Sound velocity general formulae
    # c = c0 + delta_c_temp + delta_c_sal + delta_c_press + delta_c_tsp
    # pressure unit = decibar, salinity = psu, temperatures in degree C
    # scale pressure to bars
    pressure /= 10
    sal0 = 35.0
    sal -= sal0
    c = 1449.14
    c_temp = [4.5721, -4.4532E-2, -2.6045E-4, 7.9851E-6]
    for i in range(len(c_temp)):
        c += c_temp[i] * (temp ** (i+1))

    c_sal = [1.39799, 1.69202E-3]
    for i in range(len(c_sal)):
        c += c_sal[i] * (sal ** (i+1))

    c_press = [1.60272E-1, 1.0268E-5, 3.5216E-9, -3.3603E-12]
    for i in range(len(c_press)):
        c += c_press[i] * (pressure ** (i + 1))

    a = [-1.1244E-2, 7.7711E-7, 7.7016E-5,
             -1.2943E-7, 3.1580E-8, 1.5790E-9,
             -1.8607E-4, 7.4812E-6, 4.5283E-8,
             -2.5294E-7, 1.8563E-9, -1.9646E-10]

    t2 = temp ** 2
    p2 = pressure ** 2
    pt = pressure * temp
    pt2 = pt * temp
    pt3 = pt2 * temp
    p2t = pt * pressure
    p3t = p2t * pressure
    p2t2 = pt * pt

    param = [sal * temp, sal * t2, sal * pressure, sal * p2, sal * pt, sal * pt2,
             pt, pt2, pt3, p2t, p2t2, p3t]

    for i in range(len(a)):
        c += a[i] * param[i]
    return c


def unesco_sound_velosity(sal, temp, pressure):
    # Sound velocity UNESCO formulae
    # Fofonov, 1983
    # sound_vel = c_w + AS + B (S ** 3/2) + D * (S ** 2)
    # pressure unit = decibar, salinity = psu, temperatures in degree C

    # scale pressure to bars
    pressure /= 10
    p2 = pressure ** 2
    p3 = pressure ** 3
    c_coeff = [[1402.388, 5.03711, -5.80852E-2, 3.3420E-4, -1.47800E-6, 3.1464E-9],
         [0.153563, 6.8982E-4, -8.1788E-6, 1.3621E-7, -6.1185E-10],
         [3.1260E-5, -1.7107E-6, 2.5974E-8, -2.5335E-10, 1.0405E-12],
         [-9.7729E-9, 3.8504E-10, -2.3643E-12]]
    c_w = 0.0
    for i in range(len(c_coeff)):
        for j in range(len(c_coeff[i])):
            c_w += c_coeff[i][j] * (temp ** j) * (pressure ** i)
    a_coeff = [[1.389, -1.262E-2, 7.164E-5, 2.006E-6, -3.21E-8],
         [9.4742E-5, -1.2580E-5, -6.4885E-8, 1.0507E-8, -2.0122E-10],
         [-3.9064E-7, 9.1041E-9, -1.6002E-10, 7.988E-12],
         [1.100E-10, 6.649E-12, -3.389E-13]]
    a = 0.0
    for i in range(len(a_coeff)):
        for j in range(len(a_coeff[i])):
            a += a_coeff[i][j] * (temp ** j) * (pressure ** i)
    b_coeff = [[-1.922E-2, -4.42E-5],
               [7.3637E-5, 1.7945E-7]]
    b = 0.0
    for i in range(len(b_coeff)):
        for j in range(len(b_coeff[i])):
            b += b_coeff[i][j] * (temp ** j) * (pressure ** i)
    d = 1.727E-3 - 7.9836E-6 * pressure
    return c_w + a * sal + b * (sal ** Fraction(3,2)) + d * (sal ** 2)


def calculate_press(latitude, horison):
    # p - pressure, 100 дБар
    # horizon - depth, m
    # latitude, deg
    phi = abs(latitude)
    a = -2.204E-2
    b = 0
    coeff = [99.404, 4.983E-4, -2.06E-4, 1.492E-6]
    for i in range(len(coeff)):
        b += coeff[i] * (phi ** i)
    c = -horison
    discriminant = b ** 2 - 4 * a * c
    p = 99999
    if discriminant >= 0:
        p1 = 100 * (math.sqrt(discriminant) - b) / (2 * a)
        p2 = 100 * (-b - math.sqrt(discriminant)) / (2 * a)
        if math.fabs(p1 - horison) < math.fabs(p2 - horison):
            p = p1
        else:
            p = p2
    return float("{0:.4f}".format(p))


def calculate_depth(latitude, pressure):
    # latitude, degrees
    # h - depth, meter
    # pressure, decibar
    # p, 100 дБарecibar
    phi = abs(latitude)
    p = pressure / 100
    a = [99.404, 4.983E-4, -2.06E-4, 1.492E-6]
    h = 0.0
    for i in range(len(a)):
        h += a[i] * (phi ** i)
    h -= 2.204E-2 * p
    return h * p


def linear_reg(one, two, xc, place_of_x, place_of_y):
    x1 = one[place_of_x]
    x2 = two[place_of_x]
    y1 = one[place_of_y]
    y2 = two[place_of_y]
    return y1 + (y2 - y1) / (x2 - x1) * (xc - x1)


def find_horizon(my_list, start, end, step, position, value):
    # my_list - list of lists, on place with index = position
    # should be the value - if it is found
    # the index of this element will be returned
    for k in range(start, end, step):
        if my_list[k][position] == value:
            return k
    return -1


def interpolation3(input_array, measurement_name):
    array = copy.deepcopy(input_array)
    try:
        max_depth = max(array, key=lambda x: x[2])[2]
        min_depth = min(array, key=lambda x: x[2])[2]
    except ValueError:
        print('Value Error', array)
        return array, False
    epsilon = 0.1
    standard = []
    for j in range(len(STD_HORIZONS)):
        try:
            if min_depth < STD_HORIZONS[j] < max_depth:
                standard.append(STD_HORIZONS[j])
        except TypeError:
            print(min_depth)
            print(max_depth)

    for i in range(len(array)):
        for hor in standard:
            if math.fabs(hor - array[i][2]) < epsilon:
                standard.remove(hor)
                array[i][0] = True   # flag of standard horizon
                array[i][2] = hor
                break
    if len(standard) == 0:
        return array, False
    for j in range(len(standard)):
        #                0             1          2                   3
        #             [standard interpolated   horizon             sound_vel ]
        array.append([True, True, float("{0:.2f}".format(standard[j])), 0])
    array = copy.deepcopy(sorted(array, key=lambda x: x[2], reverse=False))
    num = len(array)
    first_observed_index = find_horizon(my_list=array, start=0, end=num,
                                        step=1, position=1, value=False)
    first_standard_index = find_horizon(my_list=array, start=first_observed_index, end=num,
                                        step=1, position=1, value=False)
    for i in range(first_standard_index, num-1):
        if array[i][1]:
            # find previous observed horizon
            prev = find_horizon(my_list=array, start=i, end=-1, step=-1, position=1, value=False)
            # find next observed horizon
            next = find_horizon(my_list=array, start=i, end=num, step=1, position=1, value=False)
            first = array[prev]
            second = array[next]
            # 2 - x: depth, horizon
            # 3 - sound vel
            sv = linear_reg(one=first, two=second, place_of_x=2, place_of_y=3, xc=array[i][2])
            array[i][3] = sv
    output_just_standard = {'name': measurement_name}
    for item in array:
        if item[0]:
            output_just_standard[item[2]] = item[3]
    return output_just_standard, True


def read_datetime_with_fill_value(data, key):
    result = ''
    for item in data.variables[key][:]:
        s = str(item)[2:3]
        result += s
    d = dt.datetime(year=int(result[0:4]),
                          month=int(result[4:6]),
                          day=int(result[6:8]),
                          hour=int(result[8:10]),
                          minute=int(result[10:12]),
                          second=int(result[12:])).replace(tzinfo=None)

    return d


def get_jul_dates(dataset):
    n = dataset.dimensions['N_PROF'].size
    reference_jd_date = read_datetime_with_fill_value(dataset, 'REFERENCE_DATE_TIME')
    juld_location = []
    for i in range(n):
        juld = reference_jd_date + timedelta(seconds=int(dataset.variables['JULD_LOCATION'][i].data * 24 * 3600))
        # juld_location.append(juld.strftime(TIME_TEMPLATE))
        juld.replace(tzinfo=None)
        juld_location.append(juld)
    return juld_location


def prn_vars(dts):
    print("key; .variables[key].name; .dimensions name; [.dimensions size] ")
    for key in dts.variables:
        dim_names = ''
        dim_sizes = ''
        for dims in dts.variables[key].dimensions:
            dim_names += dts.dimensions[dims].name + '; ' + str(dts.dimensions[dims].size) + '; '
        print("{}; {}; {}".format(key, dts.variables[key].name, dim_names, dim_sizes))


def data_validation(sal, temp, pres, sal_qc, temp_qc, pres_qc):
    # check if values of temperature, salinity and pressuare are in reasonable range and of good quality
    try:
        if str(pres_qc) == '--':
            flag_pres_qc = 0
        else:
            flag_pres_qc = int(str(pres_qc)[2])
        if str(sal_qc) == '--':
            flag_sal_qc = 0
        else:
            flag_sal_qc = int(str(sal_qc)[2])
        if str(temp_qc) == '--':
            flag_temp_qc = 0
        else:
            flag_temp_qc = int(str(temp_qc)[2])

        if sal != '--' and temp != '--' and pres != '--':
            # Quality Control checking qc=4 or 9 should be ignored
            pres_qc_valid = flag_pres_qc != 4 and flag_pres_qc != 9
            psal_qc_valid = flag_sal_qc != 4 and flag_sal_qc != 9
            temp_qc_valid = flag_temp_qc != 4 and flag_temp_qc != 9
            # 0 < p < 10000 decibar
            # 0 < s < 42 psu
            # -2 < t < 40
            # Values Control
            values_valid_min = pres > 0.0 and sal > 0.0 and temp > -1.9
            values_valid_max = pres < 10000.0 and sal < 42.0 and temp < 40.0
            values_valid = values_valid_min and values_valid_max
            is_valid = pres_qc_valid and psal_qc_valid and temp_qc_valid and values_valid
        else:
            is_valid = False
    except TypeError:
        flag_sal_qc, flag_temp_qc, flag_pres_qc = 9, 9, 9
        is_valid = False
    except IndexError:
        flag_sal_qc, flag_temp_qc, flag_pres_qc = 9, 9, 9
        is_valid = False
    return is_valid


if __name__ == '__main__':
    nc_path = os.getcwd() + '\\nc\\'
    file_names = []
    if os.path.isfile(nc_path):
        print(nc_path, ' - файл')
    elif os.path.isdir(nc_path):
        for item in os.listdir(nc_path):
            if os.path.isfile(nc_path + item):
                fn_pattern = r"filtered_[\d]{8}_prof_*[\d]*\.nc"
                parce = re.search(fn_pattern, item)
                if parce:
                    file_names.append(nc_path + item)
    profiles = []
    indices = []
    count = 0
    for fn in file_names:
        print(fn)
        dataset = Dataset(fn, "r", format="NETCDF3_CLASSIC")
        n_prof = dataset.dimensions['N_PROF'].size  # profile number
        n_levels = dataset.dimensions['N_LEVELS'].size  # profile number
        juld_location = get_jul_dates(dataset)
        latitude = dataset.variables['LATITUDE']
        salinity = dataset.variables['PSAL_ADJUSTED']
        salinity_qc = dataset.variables['PSAL_ADJUSTED_QC']
        salinity_err = dataset.variables['PSAL_ADJUSTED_ERROR']
        pressure = dataset.variables['PRES_ADJUSTED']
        pressure_qc = dataset.variables['PRES_ADJUSTED_QC']
        pressure_err = dataset.variables['PRES_ADJUSTED_ERROR']
        temperature = dataset.variables['TEMP_ADJUSTED']
        temperature_qc = dataset.variables['TEMP_ADJUSTED_QC']
        temperature_err = dataset.variables['TEMP_ADJUSTED_ERROR']

        for i in range(n_prof):
            svel = []
            depth = []
            psal_adj = salinity[i]
            psal_adj_qc = salinity_qc[i]
            psal_adj_err = salinity_err[i].data
            pres_adj = pressure[i]
            pres_adj_qc = pressure_qc[i]
            pres_adj_err = pressure_err[i].data
            temp_adj = temperature[i]
            temp_adj_qc = temperature_qc[i]
            temp_adj_err = temperature_err[i].data
            svel_arr = []
            salt_arr = []
            temp_arr = []
            for j in range(n_levels):
                if data_validation(sal=psal_adj[j], temp=temp_adj[j], pres=pres_adj[j],
                                  sal_qc=psal_adj_qc[j], temp_qc=temp_adj_qc[j], pres_qc=pres_adj_qc[j]):
                    dep = calculate_depth(latitude=latitude[i], pressure=pres_adj[j])
                    sv = unesco_sound_velosity(sal=psal_adj[j], temp=temp_adj[j], pressure=pres_adj[j])
                    depth.append(dep)
                    svel.append(sv)
                    #  0             1          2            3
                    # [standard interpolated   horizon   sound_vel/salt/temp ]
                    svel_arr.append([False, False, dep, sv])
                    salt_arr.append([False, False, dep, psal_adj[j]])
                    temp_arr.append([False, False, dep, temp_adj[j]])
                else:
                    depth.append(np.nan)
                    svel.append(np.nan)

            if len(svel_arr) > 15:
                svel_standard, success_sv = interpolation3(svel_arr, 'SVEL')
                if success_sv:
                    profiles.append(svel_standard)
                    indices.append(juld_location[i])
                    count += 1
                salt_standard, success_sal = interpolation3(salt_arr, 'SALT')
                if success_sal:
                    profiles.append(salt_standard)
                    indices.append(juld_location[i])
                    count += 1
                temp_standard, success_temp = interpolation3(temp_arr, 'TEMP')
                if success_temp:
                    profiles.append(temp_standard)
                    indices.append(juld_location[i])
                    count += 1

        dataset.close()
    df = pd.DataFrame(profiles, index=indices)
    df.to_csv('nc/labrador_sea_full.csv')




