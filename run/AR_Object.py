import numpy as np


class AR_Object():
	def __init__(self):
		self.hr_time_str= None
		self.OBJECT_ID = None
		self.filename = None
		self.object_length = None
		self.object_width = None
		self.length_to_width = None
		self.mean_IVT = None
		self.object_orientation_direction = None
		self.eccentricity = None
		self.landfalling = None 
		self.landfall_point = None
		self.wind_dir_mean = None
		self.wind_dir_var = None
		self.end_lat = None
		self.end_lon = None
		self.start_lat = None
		self.start_lon = None
		self.AR_FLAG = None

	def Make_Db(self):

		from make_dbase import make_db
		dic = {'hr_time_str': self.hr_time_str,
				'OBJECT_ID': self.OBJECT_ID,
				'filename': self.filename,
				'object_length': self.object_length, 
				'object_width': self.object_width,
				'length_to_width': self.length_to_width,
				'mean_IVT': self.mean_IVT,
				'object_orientation_direction': self.object_orientation_direction,
				'eccentricity': self.eccentricity,
				'landfalling': self.landfalling, 
				'landfall_point': self.landfall_point,
				'wind_dir_mean': self.wind_dir_mean,
				'wind_dir_var': self.wind_dir_var,
				'wind_speed': self.wind_speed,			
				'end_lat': self.end_lat,
				'end_lon': self.end_lon,
				'start_lat':self.start_lat,
				'start_lon':self.start_lon,
				'AR_FLAG':self.AR_FLAG}

		make_db(**dic)

#	def Test_Flags(self):
		


	def Save_File(self, label_indices, ivt, path, v_wgt_mn, u_wgt_mn):
		# 1) Object label indices
		# 2) IVT grid 
		# 3) Maybe more later .... 

		save_me = np.zeros((4, 361, 720))
		save_me[0][label_indices] = ivt[label_indices]
		save_me[1][label_indices] = v_wgt_mn[label_indices] 
		save_me[2][label_indices] = u_wgt_mn[label_indices] 
		save_me[3] = path  # path is the entire grid so we don't need indices
		directory  = '../data/outfiles/'
		save_me_name = directory + self.hr_time_str+'_OBJECT_'+str(self.OBJECT_ID)
		np.save(save_me_name, save_me)


	def __call__(self):
		pass


