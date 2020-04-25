import os
import re

from netCDF4._netCDF4 import Dataset
import numpy.ma as ma

NC_PATH = 'C:/argo_storage/2filter/'
SETTINGS = 'settings.txt'
LOG = 'log.txt'

KEY_VAR_TO_FILTER = ['PLATFORM_NUMBER','PROJECT_NAME','PI_NAME','STATION_PARAMETERS','CYCLE_NUMBER','DIRECTION',
                 'DATA_CENTRE','DC_REFERENCE','DATA_STATE_INDICATOR','DATA_MODE','PLATFORM_TYPE','FLOAT_SERIAL_NO',
                 'FIRMWARE_VERSION','WMO_INST_TYPE','JULD','JULD_QC','JULD_LOCATION','LATITUDE','LONGITUDE',
                 'POSITION_QC','POSITIONING_SYSTEM','PROFILE_PRES_QC','PROFILE_TEMP_QC','PROFILE_PSAL_QC',
                 'VERTICAL_SAMPLING_SCHEME','CONFIG_MISSION_NUMBER','PRES','PRES_QC','PRES_ADJUSTED','PRES_ADJUSTED_QC',
                 'PRES_ADJUSTED_ERROR','TEMP','TEMP_QC','TEMP_ADJUSTED','TEMP_ADJUSTED_QC','TEMP_ADJUSTED_ERROR','PSAL',
                 'PSAL_QC','PSAL_ADJUSTED','PSAL_ADJUSTED_QC','PSAL_ADJUSTED_ERROR','PARAMETER',
                 'SCIENTIFIC_CALIB_EQUATION','SCIENTIFIC_CALIB_COEFFICIENT','SCIENTIFIC_CALIB_COMMENT',
                 'SCIENTIFIC_CALIB_DATE']

KEY_DIM_TO_CHANGE = 'N_PROF'

KEY_HISTORY = ['HISTORY_INSTITUTION','HISTORY_STEP','HISTORY_SOFTWARE','HISTORY_SOFTWARE_RELEASE','HISTORY_REFERENCE',
               'HISTORY_DATE','HISTORY_ACTION','HISTORY_PARAMETER','HISTORY_START_PRES','HISTORY_STOP_PRES',
               'HISTORY_PREVIOUS_VALUE','HISTORY_QCTEST']

if __name__ == '__main__':
    nc_path = os.getcwd() + '/'
    print(nc_path)

    if nc_path.find('Django') != -1:
        nc_path = NC_PATH

    select = {}
    with open(nc_path + SETTINGS, 'r') as g:
        keys = g.readline().replace('\n','').split(';')
        values = g.readline().replace('\n','').split(';')
        if len(keys) == len(values):
            for i in range(len(keys)):
                select[keys[i]] = float(values[i])
        else:
            print('Warning: неверные параметры выборки')
            print(keys)
            print(values)
            raise ValueError

    file_names = []
    if os.path.isfile(nc_path):
        print(nc_path, ' - файл')
    elif os.path.isdir(nc_path):
        print(nc_path, ' path found')
        for item in os.listdir(nc_path):
            if os.path.isfile(nc_path + item):
                # D20190104_prof_1.nc
                fn_pattern = r"[\d]{8}_prof_*[\d]*\.nc"  # 05/03/20 поставлен \ перед .docx
                parce = re.search(fn_pattern, item)
                if parce:
                    print(item)
                    file_names.append(item)

    lat_edge = (select.get("latitude_min"), select.get("latitude_max"))
    long_edge = (select.get("longitude_min"), select.get("longitude_max"))
    output_dir_name = 'lat[{},{}],long[{},{}]/'.format(str(lat_edge[0]), str(lat_edge[1]), str(long_edge[0]), str(long_edge[1]))
    save_path = nc_path + output_dir_name
    print(save_path)
    if not os.path.exists(save_path):
        try:
            os.mkdir(save_path, mode=0o777)
            print('directory {} was created'.format(save_path))
        except OSError:
            print('OS error occurs, a directory {} was not created'.format(save_path))
            raise OSError
    with open(save_path + LOG, 'w') as g:
        g.write('Current dir: ' + nc_path + '\n\n')
        g.write('Coordinate range: ' + output_dir_name + '\n\n')

    # If the file is open for write access (mode='w', 'r+' or 'a'), you may write any type of data including new
    # dimensions, groups, variables and attributes.
    # data_model= NETCDF3_CLASSIC
    file_count = 0
    for i in range(len(file_names)):

        prof_numbers_to_save = []
        dataset = Dataset(nc_path + file_names[i], "r", format="NETCDF3_CLASSIC")
        n_prof = dataset.dimensions['N_PROF'].size  # число профилей
        latitude = dataset.variables['LATITUDE']    # массив 1 х n_prof из широты
        longitude = dataset.variables['LONGITUDE']  # массив 1 х n_prof из широты
        with open(save_path + LOG, 'a') as g:
            g.write(file_names[i] + '\n')
            g.write('Coordinate range: Lat: {}, {}; Long: {}, {} - '.format(min(latitude), max(latitude),
                                                                            min(longitude), max(longitude)))

        for j in range(n_prof):
            if (lat_edge[0] <= latitude[j] <= lat_edge[1]) and (long_edge[0] <= longitude[j] <= long_edge[1]):
                prof_numbers_to_save.append(j)
        if len(prof_numbers_to_save) == 0:
            dataset.close()
            with open(save_path + LOG, 'a') as g:
                g.write(' 0 profiles saved\n\n')
        else:
            new_fn = 'filtered_' + file_names[i]
            filtered_set = Dataset(save_path + new_fn, 'w', format="NETCDF3_CLASSIC")
            # To copy the global attributes of the netCDF file
            for attname in dataset.ncattrs():
                setattr(filtered_set, attname, getattr(dataset, attname))
            # To copy the dimension of the netCDF file
            for dimname, dim in dataset.dimensions.items():
                # if you want to make changes in the dimensions of the new file
                # you should add your own conditions here before the creation of the dimension.
                if dimname != 'N_PROF':
                    filtered_set.createDimension(dimname, len(dim))
                else:
                    filtered_set.createDimension(dimname, len(prof_numbers_to_save))
            # To copy the variables of the netCDF file
            for varname, ncvar in dataset.variables.items():
                # if you want to make changes in the variables of the new file
                # you should add your own conditions here before the creation of the variable.
                # переписываем только те профили, которые удовлетворяют условиям выборки
                if varname in KEY_VAR_TO_FILTER:
                    ncvar_list = list(ncvar[:])
                    temp = []
                    for num in range(len(ncvar)):
                        if num in prof_numbers_to_save:
                            temp.append(ncvar_list[num])
                    var_temp = ma.array(temp, mask=ma.nomask, dtype=ncvar.dtype)
                    var = filtered_set.createVariable(varname=varname, datatype=ncvar.dtype,
                                                      dimensions=ncvar.dimensions,
                                                      fill_value=getattr(dataset.variables[varname], '_FillValue'))
                    # Finally copy the variable data to the new created variable
                    var[:] = var_temp[:]
                    # Proceed to copy the variable attributes
                    for attname in ncvar.ncattrs():
                        if attname != '_FillValue':
                            setattr(var, attname, getattr(ncvar, attname))
                    #print(varname, ncvar.dtype, ncvar.dimensions)
                elif varname in KEY_HISTORY:
                    # данные имеют первую размерность ноль их не записываем
                    pass
                else:
                    var = filtered_set.createVariable(varname=varname, datatype=ncvar.dtype,
                                                      dimensions=ncvar.dimensions,
                                                      fill_value=getattr(dataset.variables[varname], '_FillValue'))
                    # Finally copy the variable data to the new created variable
                    var[:] = ncvar[:]
                    # Proceed to copy the variable attributes
                    for attname in ncvar.ncattrs():
                        if attname != '_FillValue':
                            setattr(var, attname, getattr(ncvar, attname))

            n = filtered_set.dimensions['N_PROF'].size
            n_old = dataset.dimensions['N_PROF'].size
            dataset.close()
            filtered_set.close()
            print(new_fn, ": ", n, "profiles saved")
            with open(save_path + LOG, 'a') as g:
                g.write(str(n) + ' profiles saved from ' + str(n_old) + '\n\n')