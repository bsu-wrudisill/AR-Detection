import gc 
from skimage.measure import regionprops
import numpy as np
from centerline import Centerline
from atmDetect import *


np.seterr(all='print')

def Blob_Tester(blob_num,label,label_array,u_i, v_i, wnd_360, ivt,land_mask,lats, lons,lats_mesh, lons_mesh, hr_time_str, fname):

    ##INPUTS##
    # blob_num
    # Label
    # Label_array
    # u_i, v_i
    # ivt
    # land_mask
    # lats, lons
    # lats_mesh, lons_mesh
    # hr_time_str

 #   gc.enable()
    cell_to_km = 50
    
    #-------------------------------------------------------------------------#
    #  Test Flags; set to true if passing 
    #-------------------------------------------------------------------------#
    TC_a = False      # Mean IVT 
    TC_b = False      # Coherency in IVT Direction  (variance)
    TC_c = False      # Object Mean Meridonal IVT (flowing towards a pole)
    TC_d = False      # Object/IVT direction Consistency
    TC_e = False      # Length/Width Ratio
    TC_f = False      # Interior USA Landfalling 
    TC_g = False      # NOT crossing equator (we don't want AR crossing eq)
    #-------------------------------------------------------------------------#

    # Indicies of blob of interest
    label_indices     = np.where(label_array == label)
    label_indices_alt = np.argwhere(label_array == label)  # indices in [(a,b), ....] fmt

    # Remove all elements except blob of interest
    sub_array = np.where(label_array == label, 1, 0)

    #-------------------------------------------------------------------------#
    # Region Properties of blob of interest 
    #    
    #    From scipy.morphology.regionprops
    #    Lenght/Width Criteria                
    #    TC_e: Length must be gt. 2000km and l:w gt 2. 
    #-------------------------------------------------------------------------#

    blob          = regionprops(sub_array)[0]                    # set to 0; only 1 region        
    blob_length   = blob.major_axis_length*cell_to_km
    blob_width    = blob.minor_axis_length*cell_to_km
    blob_size     = blob.filled_area*cell_to_km

    if blob_width == 0: 
        return 

    blob_length_width_ratio = blob_length/blob_width

    #-------------------------------------------------------------------------#
    # Blob Direction Calulations;
    # We must flip blob dir (*-1) since the array is flipped w.r.t lat
    #-------------------------------------------------------------------------#
    blob_dir_raw         = blob.orientation / np.pi * 180
    blob_dir_corrected   = blob_dir_raw * -1 

    # Blob Orientation; we must deal with ambiguity of blob orientation w.r.t to wind dir
    # So we give it two possible directions
    # Recalls: orientation is between -90 and 90. So we convert that to [0,359.9]

    if  blob_dir_corrected >= 0:                                        # Positive blob dir
        blob_dir = (blob_dir_corrected, 180+blob_dir_corrected)
    else:                                                               # Negative blob dir
        blob_dir = (360 - abs(blob_dir_corrected), 360 - abs(blob_dir_corrected) - 180)


    #-------------------------------------------------------------------------#
    # Blob Width/Length 
    #-------------------------------------------------------------------------#
    
    # Break loop if there is no width... ???

    if blob_width == 0: 
        return 


    
    #-------------------------------------------------------------------------#
    #    Mean IVT Calculation
    #    
    #   TC_a: Mean IVT must be greater than 250 kg/m/s
    #-------------------------------------------------------------------------#

    mean_ivt = np.mean(ivt[label_indices])              # Mean IVT of blob


    #-------------------------------------------------------------------------#
    #    Wind Direction calculations
    #
    #    TC_b:  No more than 1/2 of gridcells in blob should deviate more than 45* 
    #           from objects mean IVT direction
    #    TC_c:  Mean wind dir and Object Orientation (blob_dir) Should not deviate by more than 45*
    #    TC_d:  Mean wind dir must have 'significant' poleward (meridonal) component;    
    #    
    #    ----------Mean Wind direction-----------
    #    Calculated by finding the mean of the u and v components respectively, then 
    #    Finding the angle between them and converting to degree coordinates
    #    And finally converting to the range [0, 359.0]
    #--------------------------------------------------------------------
    
    mean_u_i        = np.mean(u_i[label_indices])
    mean_v_i        = np.mean(v_i[label_indices])
    mean_wind_dir_  = np.arctan2(mean_v_i,mean_u_i)/np.pi * 180
    mean_wind_dir   = np.where(mean_wind_dir_< 0, 360 + mean_wind_dir_ , mean_wind_dir_)
    
    #    ----------Angular difference between object and mean wind dir-------
    angle_diff = map(lambda X: smallest_angle(X, mean_wind_dir), wnd_360[label_indices])
    
    # Counts angles gt 45 deg from mean dir
    angle_gt_mean = map(lambda x: x>45, angle_diff).count(True) 

    #    ----------Poleward IVT----------
    poleward_IVT = mean_ivt*np.sin(mean_wind_dir/180*np.pi)


    #--------------------------------------------------------------------
    # Hemispere
    #--------------------------------------------------------------------
    if any(n < 0 for n in lats[label_indices[0]]) == True:
        Hemisphere = 'Southern'

    else:
        Hemisphere = 'Northern'

        #----------------------------------------------------------------------------------#
# Landfalling 
        #----------------------------------------------------------------------------------#

    if any( n == 1 for n in land_mask[label_indices]) == True:
        Landfalling = True 
    else:
        Landfalling = False
    

    #------------------------------------------------------------------
    # Does AR make it into interior West? Between 120 W and 100 W 
    #--------------------------------------------------------------------

    if any( (n > 240) and (n < 260 ) for n in lons_mesh[label_indices]) == True:
        if any( (n > 40) and (n < 50) for n in lats_mesh[label_indices]) == True:
            Interior_Landfalling = True
        else:
            Interior_Landfalling = False
    else:
        Interior_Landfalling = False




    #-------------------------------------------------------------------------#
    # Blog Length Calc -- route through array method
    # And Landfalling Location 
    #-------------------------------------------------------------------------#

    # For description see Centerline object 
    # Centerline doesn't ingest lats_mesh or lons_mesh
    # So we do that here... for now ..... 

    try:
        center   = Centerline(label_indices_alt, label_indices, ivt, land_mask)
        
        if center.landfall_location:
            landfall_point = []
            for i in center.landfall_location:
                landfall_point.append((lats_mesh[i],lons_mesh[i]))

        else:
            landfall_point = None
    
    except Exception as e:

        Logger().failure(fname + '__' + hr_time_str)
        raise e
        return 


    #----------------------------------------------------------------------------------#
    # Object Testing 
    #----------------------------------------------------------------------------------#

    if mean_ivt > 350.0:
        TC_a = True
    else:
        return 


        #----------------------------------------------------------------------------------#
    if angle_gt_mean < len(angle_diff)/2:
        TC_b = True
    else:
        return 
                    
    #----------------------------------------------------------------------------------#
    if smallest_angle(mean_wind_dir, blob_dir[0]) < 50.0 or smallest_angle(mean_wind_dir, blob_dir[1]) < 50.0:
        TC_c = True

    else:
        return

    #----------------------------------------------------------------------------------#
    if Hemisphere == 'Northern':
        if poleward_IVT > 50.0:     
            TC_d = True
    elif poleward_IVT <  -50:
        TC_d = True 

    else:
        return

    #----------------------------------------------------------------------------------#
    if (blob_length > 2000.0) and (blob_length_width_ratio > 2.0): # Later add a width component...maybe this does not matter
        TC_e = True

    else:
        return

    #----------------------------------------------------------------------------------#
    TC_f = True         
    # Add landfalling later


    #---------------------------------------------------------------------------------#
    # Write output to dictionary or do Nothing
    #---------------------------------------------------------------------------------#

    if TC_a + TC_b + TC_c + TC_d + TC_e + TC_f == 6:  # In Python, True == 1 

        AR_EXISTS = True

        #-------------------------------------------------------------------------#
        # Write output to dictionary 
        # optional: save output plot
        #-------------------------------------------------------------------------#

        # ID of AR blob feature
        blob_num = blob_num + 1                                

        # Name of AR 
        AR_Name = 'AR_' + str(blob_num)


        # Dictionary Entry 


        info =  {
            'hr_time_str'                      : hr_time_str,
            'AR_Name'                          : AR_Name,
            'fname'                            : fname,
            'object_length'                    : float(blob_length),
            'object_width'                     : float(blob_width),
            'mean_IVT'                         : float(mean_ivt),
            'mean_wind_dir'                    : float(mean_wind_dir),
            'poleward_IVT'                     : float(poleward_IVT),
            'object_orientation_direction_a'   : float(blob_dir[0]),
            'object_orientation_direction_b'   : float(blob_dir[1]),
            'length_to_width'                  : float(blob_length_width_ratio),
            'Hemispere'                        : str(Hemisphere),       
            'Path_len'                         : float(center.path_len),
            'ID_Landfalling'                   : str(Interior_Landfalling),
            'Landfalling'                      : str(Landfalling),
            'Landfall_point'                   : str(landfall_point)
        }

        #-----------------------------------------------------------------------# 
        # import make_dbase function;
        # log AR info to sqli
        #-----------------------------------------------------------------------# 
        import make_dbase
        make_dbase.make_db(**info)
     

    else: #This is probably not needed... 
        return 

    gc.collect()



#----------------------------------------------------------------------------------------#
# END Blob-Loop
#---------------------------------------------------------------------------------#
