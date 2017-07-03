#!/usr/local/bin/python
# -*- coding: utf-8 -*-

import sqlite3 
import sys
import time
import logging


def log_AR(hr_time_str,AR_Name,fname, object_orientation_direction, object_length, object_width,  mean_IVT, mean_wind_dir, poleward_IVT, length_to_width, Hemispere, Path_len):

 dBase = sqlite3.connect('Atmospheric_River.db')
 cursor = dBase.cursor()



 cursor.execute('''CREATE TABLE Rivers(
 hr_time_str TEXT,
 AR_Name TEXT, 
 fname TEXT,
 object_orientation_direction REAL,
 object_length REAL,
 object_width REAL,
 mean_IVT REAL,
 mean_wind_dir REAL,
 poleward_IVT REAL,
 length_to_width REAL,
 Hemispere TEXT,
 Path_len REAL)'''
 )


 cursor.execute('''INSERT INTO Rivers(
 hr_time_str, 
 AR_Name, 
 fname, 
 object_orientation_direction,
 object_length,
 object_width,
 mean_IVT,
 mean_wind_dir,
 poleward_IVT,
 length_to_width,
 Hemispere,
 Path_len
 )

 VALUES(?,?,?,?,?,?,?,?,?,?,?,?)''', (hr_time_str, AR_Name, fname, object_orientation_direction, object_length, object_width, mean_IVT, mean_wind_dir, poleward_IVT, length_to_width, Hemispere, Path_len))

 dBase.commit()
 dBase.close()
 logging.warning('success')
