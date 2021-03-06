#!/usr/local/bin/python
# -*- coding: utf-8 -*-

import sqlite3 
import sys
import time
import logging



def make_db(hr_time_str,
            OBJECT_ID,
            filename,
            object_length,
            object_width,
            length_to_width,
            mean_IVT,
            meridional_IVT,
            object_orientation_direction,
            eccentricity,
            landfalling ,
            landfall_lat,
            landfall_lon,
            wind_dir_mean,
            wind_dir_var,
            wind_speed,  
            end_lat,
            end_lon,
            start_lat,
            start_lon,
            gc_distance,
            AR_BASE_FLAG):

  while True:

    dBase = sqlite3.connect('Atmospheric_River.db',timeout=10)
    cursor = dBase.cursor()

    try:
      cursor.execute('''CREATE TABLE IF NOT EXISTS Rivers(
                              hr_time_str TEXT, 
                              OBJECT_ID TEXT,
                              filename TEXT,
                              object_length REAL,
                              object_width REAL,
                              length_to_width REAL,
                              mean_IVT REAL,
                              meridional_IVT REAL,
                              object_orientation_direction REAL,
                              eccentricity REAL,
                              landfalling TEXT,   
                              landfall_lat REAL,
                              landfall_lon REAL,
                              wind_dir_mean REAl,
                              wind_dir_var REAL, 
                              wind_speed REAL,
                              end_lat REAL,
                              end_lon REAL,
                              start_lat REAL,
                              start_lon REAL,
                              gc_distance REAL,
                              AR_BASE_FLAG TEXT
                              )''')


#
 
      cursor.execute('''INSERT INTO Rivers(
                                hr_time_str,
                                OBJECT_ID,
                                filename,
                                object_length,
                                object_width,
                                length_to_width,
                                mean_IVT,
                                meridional_IVT, 
                                object_orientation_direction,
                                eccentricity,
                                landfalling ,
                                landfall_lat,
                                landfall_lon,                          
                                wind_dir_mean,
                                wind_dir_var,
                                wind_speed,
                                end_lat,
                                end_lon,                        
                                start_lat,
                                start_lon,
                                gc_distance,
                                AR_BASE_FLAG
                                )


                       VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''', 

                       (      hr_time_str,
                              OBJECT_ID,
                              filename,
                              object_length,
                              object_width,
                              length_to_width,
                              mean_IVT,
                              meridional_IVT, 
                              object_orientation_direction,
                              eccentricity,
                              landfalling,
                              landfall_lat,
                              landfall_lon,                          
                              wind_dir_mean,
                              wind_dir_var,
                              wind_speed,                      
                              end_lat,
                              end_lon,
                              start_lat,
                              start_lon,
                              gc_distance,
                              AR_BASE_FLAG
                                ))
      dBase.commit()
      dBase.close()
      break

    except sqlite3.OperationalError as E:
      print E
      print 'database is locked. waiting'


# bax = {'Filename': 'a',
#         'hr_time_str': 'b',
#         'OBJECT_ID': 'c',
#         'object_length': 'd', 
#         'object_width': 'e',
#         'length_to_width': 'f',
#         'mean_IVT': 'g',
#         'object_orientation_direction': 'h',
#         'Path_len': 'i',
#         'eccentricity': None,
#         'Landfalling': 'z', 
#         'Landfall_point': 'x',
#         'End_Point': 'c'}

# make_db(**bax)




