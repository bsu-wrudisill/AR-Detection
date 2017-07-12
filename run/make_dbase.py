#!/usr/local/bin/python
# -*- coding: utf-8 -*-

import sqlite3 
import sys
import time
import logging



def make_db(
  hr_time_str, 
  AR_Name,
  fname, 
  object_length, 
  object_width, 
  mean_IVT, 
  mean_wind_dir, 
  poleward_IVT,
  object_orientation_direction_a,
  object_orientation_direction_b, 
  length_to_width, 
  Hemispere, 
  Path_len, 
  ID_Landfalling, 
  Landfalling, 
  Landfall_point
):

  dBase = sqlite3.connect('Atmospheric_River.db')
  cursor = dBase.cursor()
  cursor.execute('''CREATE TABLE IF NOT EXISTS Rivers(
                            hr_time_str TEXT, 
                            AR_Name TEXT,
                            fname TEXT,
                            object_length REAL,
                            object_width REAL,
                            mean_IVT REAL,
                            mean_wind_dir REAL,
                            poleward_IVT REAL,
                            length_to_width REAL,
                            object_orientation_direction_a REAL,
                            object_orientation_direction_b REAL,
                            Hemispere TEXT,
                            Path_len REAL,
                            ID_Landfalling TEXT,
                            Landfalling TEXT,   
                            Landfall_point TEXT
                )''')


#
 
  cursor.execute('''INSERT INTO Rivers(
                                      hr_time_str,
                                      AR_Name, 
                                      fname,
                                      object_length,
                                      object_width,
                                      mean_IVT,
                                      mean_wind_dir,
                                      poleward_IVT,
                                      length_to_width,
                                      object_orientation_direction_a,
                                      object_orientation_direction_b,
                                      Hemispere,
                                      Path_len,
                                      ID_Landfalling,
                                      Landfalling,
                                      Landfall_point) 

                   VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''', 

                   (hr_time_str, 
                    AR_Name,
                    fname, 
                    object_length, 
                    object_width, 
                    mean_IVT, 
                    mean_wind_dir, 
                    poleward_IVT,
                    object_orientation_direction_a,
                    object_orientation_direction_b, 
                    length_to_width, 
                    Hemispere, 
                    Path_len, 
                    ID_Landfalling, 
                    Landfalling, 
                    Landfall_point
                    ))



  dBase.commit()
  dBase.close()




#make_db('foo','foo','foo','foo','foo','foo','foo','foo','foo','foo','foo','foo','foo','foo','foo','foo')



