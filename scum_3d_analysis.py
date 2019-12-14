import numpy as np 
import matplotlib.pyplot as plt 
import sys 
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import cv2 as opencv

VICON_SCALE = 1

class Lighthouse:
	def __init__(self,x=0,y=0,z=0,pose=np.eye(3)):

		self.pose = pose
		         

		self.translation = np.array([[1,0,0,-x],
		                         [0,1,0,-y],
		                         [0,0,1,-z]])
		self.K = np.array([[1,0,0],
		               [0,1,0],
		               [0,0,1]])

		self.P = np.matmul(np.matmul(self.K,self.pose),self.translation)

		self.cvRot = None
		self.cvTrans = None
		self.cvCam = None


	def setOpenCVParams(self, tvec, rvec, K):
		self.cvRvec = rvec
		self.cvTvec = tvec
		self.cvK = K


		rot, jacobian = opencv.Rodrigues(rvec)
		self.pose = np.squeeze(np.array(rot)) 

		modified_t =  -np.linalg.inv(self.pose) @ np.squeeze(tvec)

		self.translation = np.array([[1,0,0,-modified_t[0]], [0,1,0,-modified_t[1]], [0,0,1,-modified_t[2]]]) 

		print("Pose")
		print(self.pose)
		self.K = np.squeeze(np.array(K))
		self.P = self.pose @ self.translation

#this class generates xyz coordinates from lighthouse azimuths
class Triangulator:  
    def __init__(self,lighthouse1, lighthouse2):
        self.P1 = lighthouse1.P
        self.P2 = lighthouse2.P
        
        self.xcam1 = None
        self.ycam1 = None

        '''
        self.pose1 = np.array([[1,0,0], 
                     [0,1, 0],
                     [0,0, 1]])

        self.pose2 = np.array([[ 0.60669553,  0.12838786,  0.78449798],
                             [-0.07049024,  0.99167271, -0.10777928],
                                [-0.79180279,  0.01008975,  0.61069349]])

        self.P1 = np.matmul(self.pose1,lighthouse1.translation)
        self.P2 = np.matmul(self.pose2,lighthouse2.translation)
        '''

    def triangulate(self,az1, el1, az2, el2):
        #calculate triangulation (remember a on lighthouse is actually b
        #here

        xcam1 = np.tan(az1 - np.pi/2) #focal distance of 1 (virtual)
        ycam1 = np.tan((el1 - np.pi/2))
        
        self.xcam1 = xcam1
        self.ycam1 = ycam1

        xcam2 = np.tan(az2 - np.pi/2); #focal distance of 1 (virtual)
        ycam2 = np.tan((el2 - np.pi/2))
        
        #stack A matrix (entry for each lighthouse is X x (PX) = 0) direct linear transform
        A = [xcam1 * self.P1[2,:] - self.P1[0,:],
             ycam1 * self.P1[2,:] - self.P1[1,:],
             xcam2 * self.P2[2,:] - self.P2[0,:],
             ycam2 * self.P2[2,:] - self.P2[1,:]]
        
        #find SVD 
        u,s,vt = np.linalg.svd(A)
        v= np.transpose(vt)
        #solution is final column of V
        xhat = v[:,3]/v[3,3]
        return (np.squeeze(xhat)[0],np.squeeze(xhat)[1],np.squeeze(xhat)[2])    
        

def plot_raw_data(azimuth1_gnd, azimuth2_gnd, elevation1_gnd, elevation2_gnd, lighthouse_df, scum_gnd_df):
	#################################
	plt.figure()
	plt.subplot(2,2,1)
	plt.plot(-azimuth1_gnd*180/3.14159 +90)
	plt.ylim([50,150])
	plt.title('Az 1')

	plt.subplot(2,2,3)
	#print(lighthouse_df['azA'])

	plt.plot(lighthouse_df[ lighthouse_df['azA'] > -360 ]['azA'])
	plt.ylim([50,150])
	plt.title('Az A')

	plt.subplot(2,2,2)
	plt.plot(-azimuth2_gnd*180/3.14159 +90)
	plt.ylim([50,150])
	plt.title('Az 2')

	plt.subplot(2,2,4)
	#print(lighthouse_df['azA'])
	plt.plot(lighthouse_df[ lighthouse_df['azB'] > -360 ]['azB'])
	plt.ylim([50,150])
	plt.title('Az B')

#############
	plt.figure()
	plt.subplot(2,2,1)
	plt.plot(-elevation1_gnd*180/3.14159 +90)
	plt.ylim([50,150])
	plt.title('el 1')

	plt.subplot(2,2,3)
	#print(lighthouse_df['azA'])
	plt.plot(lighthouse_df[ lighthouse_df['elA'] > -360 ]['elA'])
	plt.ylim([50,150])
	plt.title('el A')

	plt.subplot(2,2,2)
	plt.plot(scum_gnd_df['Time','Time'],-elevation2_gnd*180/3.14159 +90)
	plt.ylim([50,150])
	plt.title('el 2')

	plt.subplot(2,2,4)
	#print(lighthouse_df['azA'])
	plt.plot(lighthouse_df[ lighthouse_df['elB'] > -360 ]['timestamp_optitrack'], lighthouse_df[ lighthouse_df['elB'] > -360 ]['elB'])
	plt.ylim([50,150])
	plt.title('el b')

#returns processed and cleaned data dict from lighthouse and ground truth data
def generate_data_dict(azimuth1_gnd, azimuth2_gnd, elevation1_gnd, elevation2_gnd, lighthouse_df, scum_gnd_df):
	lighthouse_processed_dict = {}

	lighthouse_processed_dict['azA'] = {}
	lighthouse_processed_dict['azA']['time'] = lighthouse_df[ lighthouse_df['azA'] > -360 ]['timestamp_optitrack']
	lighthouse_processed_dict['azA']['lh'] = lighthouse_df[ lighthouse_df['azA'] > -360 ]['azA']
	lighthouse_processed_dict['azA']['truth']  = np.interp(lighthouse_processed_dict['azA']['time'], scum_gnd_df['Time','Time'],azimuth1_gnd*180/3.14159 +90) 

	lighthouse_processed_dict['elA'] = {}
	lighthouse_processed_dict['elA']['time'] = lighthouse_df[ lighthouse_df['elA'] > -360 ]['timestamp_optitrack']
	lighthouse_processed_dict['elA']['lh'] = lighthouse_df[ lighthouse_df['elA'] > -360 ]['elA']
	lighthouse_processed_dict['elA']['truth']  = np.interp(lighthouse_processed_dict['elA']['time'], scum_gnd_df['Time','Time'],elevation1_gnd*180/3.14159 +90) 

	lighthouse_processed_dict['azB'] = {}
	lighthouse_processed_dict['azB']['time'] = lighthouse_df[ lighthouse_df['azB'] > -360 ]['timestamp_optitrack']
	lighthouse_processed_dict['azB']['lh'] = lighthouse_df[ lighthouse_df['azB'] > -360 ]['azB']
	lighthouse_processed_dict['azB']['truth']  = np.interp(lighthouse_processed_dict['azB']['time'], scum_gnd_df['Time','Time'],azimuth2_gnd*180/3.14159 +90) 

	lighthouse_processed_dict['elB'] = {}
	lighthouse_processed_dict['elB']['time'] = lighthouse_df[ lighthouse_df['elB'] > -360 ]['timestamp_optitrack']
	lighthouse_processed_dict['elB']['lh'] =  lighthouse_df[ lighthouse_df['elB'] > -360 ]['elB']
	lighthouse_processed_dict['elB']['truth']  = np.interp(lighthouse_processed_dict['elB']['time'], scum_gnd_df['Time','Time'],elevation2_gnd*180/3.14159 +90) 

	return lighthouse_processed_dict


def plot_truth_traj(scum_df,lh1, lh2):
	plt.figure()
	plt.subplot(3,1,1)
	plt.plot(scum_df['Time', 'Time'],scum_df['Position','X']/VICON_SCALE)

	plt.subplot(3,1,2)
	plt.plot(scum_df['Time', 'Time'],scum_df['Position','Y']/VICON_SCALE)

	plt.subplot(3,1,3)
	plt.plot(scum_df['Time', 'Time'],scum_df['Position','Z']/VICON_SCALE)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot(xs = scum_df['Position','X']/VICON_SCALE, ys = scum_df['Position','Y']/VICON_SCALE, zs = scum_df['Position','Z']/VICON_SCALE)
	ax.scatter(xs = [lh1[1][0], lh2[1][0]], 
				ys = [lh1[1][1], lh2[1][1]], 
				zs = [lh1[1][2], lh2[1][2]])
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	ax.axis('scaled')

#this function will project the ground truth motions to azimuth and elevation angles relative to the lighthouse. 
def project_truth_to_az_el(scum_df,lh1,lh2,lh1_to_gnd,lh2_to_gnd):

	scum_pos = scum_df['Position']
	#rotate ground truth xyz to lh1 coordinates

	traj_lh1 = np.dot(np.linalg.inv(lh1_to_gnd),(scum_pos[['X','Y','Z']] - lh1).T)

	print('######## Trajectory in LH 1 Coordinate Frame')
	print(traj_lh1)
	#find azimuth and elevation using atan2
	elevation1 = np.arctan2(traj_lh1[1,:],traj_lh1[2,:])
	azimuth1 = np.arctan2(traj_lh1[0,:],traj_lh1[2,:])


	#rotate ground truth xyz to lh1 coordinates
	traj_lh2 = np.dot(np.linalg.inv(lh2_to_gnd),(scum_pos[['X','Y','Z']] - lh2).T)

	#find azimuth and elevation using atan2
	elevation2 = np.arctan2(traj_lh2[1,:],traj_lh2[2,:])
	azimuth2 = np.arctan2(traj_lh2[0,:],traj_lh2[2,:])


	return azimuth1, elevation1, azimuth2, elevation2 

def process_lh_locs(lh_locs,board_to_gnd):
	#translations are from board to lighthouse in lighthouse coordinates
	#rotations are lighthouse to obard
	lh1_row = lh_locs.loc['Lighthouse1']
	R1 = np.array([ [lh1_row['R00'], lh1_row['R01'], lh1_row['R02']],
				    [lh1_row['R10'], lh1_row['R11'], lh1_row['R12']],
				    [lh1_row['R20'], lh1_row['R21'], lh1_row['R22']] ])

	lh2_row = lh_locs.loc['Lighthouse2']
	R2 = np.array([ [lh2_row['R00'], lh2_row['R01'], lh2_row['R02']],
				    [lh2_row['R10'], lh2_row['R11'], lh2_row['R12']],
				    [lh2_row['R20'], lh2_row['R21'], lh2_row['R22']]])


	#negate translation
	lh1_trans_lh1 = -1*np.array(lh1_row[['X','Y','Z']]).T
	lh2_trans_lh2 = -1*np.array(lh2_row[['X','Y','Z']]).T

	lh1_to_board = np.linalg.inv(R1)
	lh2_to_board = np.linalg.inv(R2)

	#rotate to board coordinates
	lh1_trans_board = np.dot(lh1_to_board, lh1_trans_lh1)
	lh2_trans_board = np.dot(lh2_to_board, lh2_trans_lh2)

	lh1_trans_ground = np.dot(board_to_gnd, lh1_trans_board)
	lh2_trans_ground = np.dot(board_to_gnd, lh2_trans_board)

	#rotate to ground coordinates
	lh1_to_gnd = np.dot(board_to_gnd, lh1_to_board)
	lh2_to_gnd = np.dot(board_to_gnd, lh2_to_board)

	lh1_info = (lh1_to_gnd,lh1_trans_ground)
	lh2_info = (lh2_to_gnd,lh2_trans_ground)

	return lh1_info, lh2_info

#this function will align the lighthouse data to the vicon data. It will
#return lighthouse data timestamp column that is temporally aligned with respect to the opititrack data
def align_data(lh_data,optitrack_data = None):

	#time offsets (lh - optitrack)
	offsets = {'take1' : 284.6-85.4039,
				'take2' : 257.359 - 2.94922,
				'take3' : 58.08 - 5.671,
				'take4' : 150.123 - 96.4074,
				'take5' : 116.671-84.6218}

	aligned_timestamps = lh_data['timestamp_seconds'] - offsets[sys.argv[1]] 
	return aligned_timestamps

def triangulate_scum(lh_data,lighthouse1, lighthouse2):

	#create triangulator
	triangulator = Triangulator(lighthouse1,lighthouse2)

	#we triangulate every time we get a new measurement, which means we need to use the most recent 
	#for the calculation values of the other angles
	
	current_indicies = [0,0,0,0]

	azA_len = len(lh_data['azA']['lh'])
	elA_len = len(lh_data['elA']['lh'])
	azB_len = len(lh_data['azB']['lh'])
	elB_len = len(lh_data['elA']['lh'])

	dict_keys = ['azA', 'elA', 'azB', 'elB']
	
	last_values = [-1000,-1000,-1000,-1000]
	traj = []

	#iterate through every timestep
	for i in range(0,azA_len+elA_len+azB_len+elB_len):
		#find smallest current time
		current_times = [lh_data['azA']['time'].iloc[current_indicies[0]], 
						 lh_data['elA']['time'].iloc[current_indicies[1]], 
						 lh_data['azB']['time'].iloc[current_indicies[2]], 
						 lh_data['elB']['time'].iloc[current_indicies[3]]]
		#print(current_times)
		#replace negative or zero current times with positive inf
		if len([ time for time in current_times if time <= 0]  ) > 0:
			for index in range(0,4):
				if current_times[index] <= 0:
					current_times[index] = float('inf')
					current_indicies[index] += 1

		#find current time to use
		min_time_idx = np.argmin(current_times)
		current_time = current_times[min_time_idx]
		current_measure = lh_data[dict_keys[min_time_idx]]['lh'].iloc[current_indicies[min_time_idx]]
		#print(min_time_idx)
		#print( lh_data['azA']['truth'])
		last_values[min_time_idx] = current_measure/180.0*np.pi
		#print(last_values)
		if -1000 not in last_values and not np.isnan(last_values).any():
			point = triangulator.triangulate(az1=last_values[0], el1 = last_values[1], az2 = last_values[2], el2 = last_values[3])
			traj.append([current_time,point[0],point[1],point[2]])
			#print([current_time,current_measure,point[0],point[1],point[2]])
		#print('time: ', current_times[min_time_idx], 'Measurement: ',dict_keys[min_time_idx], ', ', current_measure , current_times)
		current_indicies[min_time_idx] += 1
	#print(traj)

	return np.squeeze(np.array(traj))

def triangulate_gnd(gnd_data,lighthouse1, lighthouse2):

	#create triangulator
	triangulator = Triangulator(lighthouse1,lighthouse2)

	#we triangulate every time we get a new measurement, which means we need to use the most recent 
	#for the calculation values of the other angles
	
	current_indicies = [0,0,0,0]
	df = pd.DataFrame(gnd_data)

	traj = []

	#iterate through every timestep
	for i in range(0,len(df)):
		azA = df['az1'].iloc[i]
		elA = df['el1'].iloc[i]
		azB = df['az2'].iloc[i]
		elB = df['el2'].iloc[i]
		if not np.isnan([azA,elA,azB,elB]).any():
			point = triangulator.triangulate(az1 = -azA, el1 = -elA, az2 = -azB, el2 = -elB)
			traj.append([df['time'].iloc[i], point[0],point[1],point[2]])
	

	return np.squeeze(np.array(traj))

#lh_data is the raw lighthouse dataframe
def triangulate_scum_nodict(lh_data,lighthouse1, lighthouse2):

	#create triangulator
	triangulator = Triangulator(lighthouse1,lighthouse2)

	#we triangulate every time we get a new measurement, which means we need to use the most recent 
	#for the calculation values of the other angles
	
	azA = [np.nan]
	azB = [np.nan]
	elA = [np.nan]
	elB = [np.nan]

	print(lh_data)
	xcam1 = [np.nan]
	ycam1 = [np.nan] 
	xcam2 = [np.nan]
	ycam2 = [np.nan]

	traj = []
	cam_points = []
	#iterate through every timestep
	for i in range(0,len(lh_data)):

		if abs(lh_data['azA'].iloc[i]) < 180:
			azA.append(lh_data['azA'].iloc[i] * np.pi/180.0)
			xcam1.append( np.tan(azA[-1] - np.pi/2) )#focal distance of 1 (virtual)
	
		if abs(lh_data['azB'].iloc[i]) < 180:
			azB.append(lh_data['azB'].iloc[i] * np.pi/180.0)
			xcam2.append( np.tan(azB[-1]- np.pi/2) )#focal distance of 1 (virtual)

		if abs(lh_data['elA'].iloc[i]) < 180:
			elA.append(lh_data['elA'].iloc[i] * np.pi/180.0)
			ycam1.append(np.tan((elA[-1] - np.pi/2)))

		if abs(lh_data['elB'].iloc[i]) < 180:
			elB.append(lh_data['elB'].iloc[i] * np.pi/180.0)
			ycam2.append( np.tan((elB[-1] - np.pi/2)))

		if not np.isnan([azA[-1],elA[-1],azB[-1],elB[-1]]).any():
			point = triangulator.triangulate(az1 = azA[-1], el1 = elA[-1], az2 = azB[-1], el2 = elB[-1])
			traj.append([lh_data['timestamp_optitrack'].iloc[i], point[0],point[1],point[2]])
			cam_points.append([lh_data['timestamp_optitrack'].iloc[i], xcam1[-1],ycam1[-1],xcam2[-1],ycam2[-1]])

	return np.squeeze(np.array(traj)), np.squeeze(np.array(cam_points))



def plot_az_el(lighthouse_processed_dict):
	plt.figure()
	plt.subplot(2,2,1)
	plt.scatter(lighthouse_processed_dict['azA']['time'],lighthouse_processed_dict['azA']['lh'],s = 1)
	plt.scatter(lighthouse_processed_dict['azA']['time'],lighthouse_processed_dict['azA']['truth'], s = 1)
	#plt.scatter(scum_gnd_df['Time','Time'], - azimuth1_gnd*180/3.14159 +90, s = 1)
	plt.ylim([0,200])
	plt.xlim([0,lighthouse_processed_dict['azA']['time'].iloc[len(lighthouse_processed_dict['azA']['time'])-1]])
	plt.title('Azimuth 1')
	plt.legend(['Lighthouse Measurement', 'Ground Truth'])
	plt.xlabel('Time (s)')
	plt.ylabel('Degrees')


	plt.subplot(2,2,2)
	plt.scatter(lighthouse_processed_dict['elA']['time'],lighthouse_processed_dict['elA']['lh'],s = 1)
	plt.scatter(lighthouse_processed_dict['elA']['time'],lighthouse_processed_dict['elA']['truth'], s = 1)
	#plt.scatter(scum_gnd_df['Time','Time'], - elevation1_gnd*180/3.14159 +90, s = 1)
	plt.ylim([0,200])
	plt.xlim([0,lighthouse_processed_dict['elA']['time'].iloc[len(lighthouse_processed_dict['elA']['time'])-1]])
	plt.title('Elevation 1')
	plt.legend(['Lighthouse Measurement', 'Ground Truth'])
	plt.xlabel('Time (s)')
	plt.ylabel('Degrees')


	plt.subplot(2,2,3)
	plt.scatter(lighthouse_processed_dict['azB']['time'],lighthouse_processed_dict['azB']['lh'],s = 1)
	plt.scatter(lighthouse_processed_dict['azB']['time'],lighthouse_processed_dict['azB']['truth'], s = 1)
	#plt.scatter(scum_gnd_df['Time','Time'], - azimuth2_gnd*180/3.14159 +90, s = 1)
	plt.ylim([0,200])
	plt.xlim([0,lighthouse_processed_dict['azB']['time'].iloc[len(lighthouse_processed_dict['azB']['time'])-1]])
	plt.title('Azimuth 2')
	plt.legend(['Lighthouse Measurement', 'Ground Truth'])
	plt.xlabel('Time (s)')
	plt.ylabel('Degrees')

	plt.subplot(2,2,4)
	plt.scatter(lighthouse_processed_dict['elB']['time'],lighthouse_processed_dict['elB']['lh'],s = 1)
	plt.scatter(lighthouse_processed_dict['elB']['time'],lighthouse_processed_dict['elB']['truth'], s = 1)
	#plt.scatter(scum_gnd_df['Time','Time'], - elevation2_gnd*180/3.14159 +90, s = 1)
	plt.ylim([0,200])
	plt.xlim([0,lighthouse_processed_dict['elB']['time'].iloc[len(lighthouse_processed_dict['elB']['time'])-1]])
	plt.title('Elevation 2')
	plt.legend(['Lighthouse Measurement', 'Ground Truth'])
	plt.xlabel('Time (s)')
	plt.ylabel('Degrees')

	#error##############################################################################
	plt.figure()
	plt.subplot(2,2,1)
	plt.scatter(lighthouse_processed_dict['azA']['time'],lighthouse_processed_dict['azA']['lh']-lighthouse_processed_dict['azA']['truth'],s = 1)
	#plt.ylim([0,200])
	plt.xlim([0,lighthouse_processed_dict['azA']['time'].iloc[len(lighthouse_processed_dict['azA']['time'])-1]])
	plt.title('Azimuth 1')
	plt.ylabel('Error (deg)')

	plt.subplot(2,2,2)
	plt.scatter(lighthouse_processed_dict['elA']['time'],lighthouse_processed_dict['elA']['lh'] - lighthouse_processed_dict['elA']['truth'],s = 1)
	#plt.ylim([0,200])
	plt.xlim([0,lighthouse_processed_dict['elA']['time'].iloc[len(lighthouse_processed_dict['elA']['time'])-1]])
	plt.title('Elevation 1')
	plt.ylabel('Error (deg)')


	plt.subplot(2,2,3)
	plt.scatter(lighthouse_processed_dict['azB']['time'],lighthouse_processed_dict['azB']['lh']-lighthouse_processed_dict['azB']['truth'],s = 1)
	#plt.ylim([0,200])
	plt.xlim([0,lighthouse_processed_dict['azB']['time'].iloc[len(lighthouse_processed_dict['azB']['time'])-1]])
	plt.title('Azimuth 2')
	plt.ylabel('Error (deg)')

	plt.subplot(2,2,4)
	plt.scatter(lighthouse_processed_dict['elB']['time'],lighthouse_processed_dict['elB']['lh'] -lighthouse_processed_dict['elB']['truth'],s = 1)
	#plt.ylim([0,200])
	plt.xlim([0,lighthouse_processed_dict['elB']['time'].iloc[len(lighthouse_processed_dict['elB']['time'])-1]])
	plt.title('Elevation 2')
	plt.ylabel('Error (deg)')

	error_azA = (lighthouse_processed_dict['azA']['lh']-lighthouse_processed_dict['azA']['truth']).values
	error_azB = (lighthouse_processed_dict['azB']['lh']-lighthouse_processed_dict['azB']['truth']).values
	error_elA = (lighthouse_processed_dict['elA']['lh']-lighthouse_processed_dict['elA']['truth']).values
	error_elB = (lighthouse_processed_dict['elB']['lh']-lighthouse_processed_dict['elB']['truth']).values

	error_azA = error_azA[~np.isnan(error_azA)]
	error_azB = error_azB[~np.isnan(error_azB)]
	error_elA = error_elA[~np.isnan(error_elA)]
	error_elB = error_elB[~np.isnan(error_elB)]

	#histogram of error
	
	plt.figure()
	plt.subplot(2,2,1)
	plt.hist(error_azA,bins = 50)
	plt.title('Azimuth 1')
	plt.xlabel('Error (deg)')

	plt.subplot(2,2,2)
	plt.hist(error_elA,bins = 50)
	plt.title('Elevation 1')
	plt.xlabel('Error (deg)')

	plt.subplot(2,2,3)
	plt.hist(error_azB,bins = 50)
	plt.title('Azimuth 2')
	plt.xlabel('Error (deg)')

	plt.subplot(2,2,4)
	plt.hist(error_elB,bins = 50)
	plt.title('Elevation 2')
	plt.xlabel('Error (deg)')	

#take in azimuth and elevation in degrees 
#this solves the equation X_gnd * P = Y_lh, solving for the P that minimizes error
def calibrate(scum_gnd_df,initial_cam1, initial_cam2,lighthouse_cams):



	begin = 2000
	end = 3500

	#loop through each point 
	
	scum_gnd_df.set_index(scum_gnd_df['Time','Time'].values, inplace = True)
	#print(scum_gnd_df['Time','Time'])
	#print(scum_gnd_df.index)
	Pcam = initial_cam1.copy()
	
	gnd_cam = []
	lh_cam_1 = [ list([row[1],row[2]]) for row in lighthouse_cams[100:400] ]
	lh_cam_1 = []
	
	#print(lh_cam_1)
	for i in range(begin,end):
		current_point_time = lighthouse_cams[i,0]
		gnd_index = scum_gnd_df.index.get_loc(current_point_time,method = 'nearest')
		gnd_row = scum_gnd_df.iloc[gnd_index]
		gnd_x = np.interp(current_point_time, scum_gnd_df['Time','Time'],scum_gnd_df['Position','X']) 
		gnd_y = np.interp(current_point_time, scum_gnd_df['Time','Time'],scum_gnd_df['Position','Y']) 
		gnd_z = np.interp(current_point_time, scum_gnd_df['Time','Time'],scum_gnd_df['Position','Z']) 
		
		gnd_point = list(gnd_row['Position'].values[0:3].astype('float32'))
		if (not np.isnan(gnd_point).any()) and (not (abs(lighthouse_cams[i,1:3]) > 5).any())  and (not (abs(lighthouse_cams[i,0]) > 900)) and (not (lighthouse_cams[i,0] < 0)):
			lh_cam_1.append([lighthouse_cams[i][1], lighthouse_cams[i][2]])
			gnd_cam.append(gnd_point)

	mat = np.array([[1,0,0,0],
				[0,1,0,0],
				[0,0,1,0],
				[0,0,0,0]])

	#print(np.array([gnd_cam]))
	#print(np.array([lh_cam_1]).astype('float32'))
	K = np.array([[1,0,0],
				 [0,1,0],
				 [0,0,1] ],dtype='float32')

	retval, rvec, tvec = opencv.solvePnP(np.array([gnd_cam]),np.array([lh_cam_1]).astype('float32'),K,None)
	print("pnp")
	print(retval)
	print(rvec)
	print(tvec)

	print("lh1 rotation mat")
	rot, jacobian = opencv.Rodrigues(rvec)
	print(rot)
	#print(np.linalg.inv(rot) @ tvec)
	print(tvec)
	plt.figure()
	plt.subplot(2,1,1)
	plt.plot(gnd_cam)
	plt.subplot(2,1,2)
	plt.plot(lh_cam_1)

	lh1_obj = Lighthouse()
	lh1_obj.setOpenCVParams(tvec = tvec,rvec = rvec,K = K)

	imagePoints, jacob = opencv.projectPoints(np.array([gnd_cam]),lh1_obj.cvRvec,lh1_obj.cvTvec,lh1_obj.cvK,None)



	

	cvproj1 = np.squeeze(np.array(imagePoints))

	rot, rando = opencv.Rodrigues(rvec)
	xyz_camera = rot @ (np.array(gnd_cam).T + np.linalg.inv(rot) @ lh1_obj.cvTvec)

	plt.figure()
	plt.subplot(2,1,1)
	plt.plot(xyz_camera[0,:]/xyz_camera[2,:])
	plt.plot(np.array(lh_cam_1)[:,0])
	plt.title("Lighthouse 1 projecttion check")
	plt.plot(cvproj1[:,0])
	plt.legend(['Lighthouse','OpenCV Projectect Truth'])

	plt.subplot(2,1,2)
	plt.plot(np.array(lh_cam_1)[:,1])
	plt.plot(cvproj1[:,1])
	plt.legend(['Lighthouse','OpenCV Projectect Truth'])



###################################################################
	gnd_cam = []
	lh_cam_2 = []
	
	#print(lh_cam_1)
	for i in range(begin,end):
		current_point_time = lighthouse_cams[i,0]
		gnd_index = scum_gnd_df.index.get_loc(current_point_time,method = 'nearest')
		gnd_row = scum_gnd_df.iloc[gnd_index]
		gnd_x = np.interp(current_point_time, scum_gnd_df['Time','Time'],scum_gnd_df['Position','X']) 
		gnd_y = np.interp(current_point_time, scum_gnd_df['Time','Time'],scum_gnd_df['Position','Y']) 
		gnd_z = np.interp(current_point_time, scum_gnd_df['Time','Time'],scum_gnd_df['Position','Z']) 
		
		gnd_point = list(gnd_row['Position'].values[0:3].astype('float32'))
		if (not np.isnan(gnd_point).any()) and (not (abs(lighthouse_cams[i,3:5]) > 5).any()) and (not (abs(lighthouse_cams[i,0]) > 900)) and (not (lighthouse_cams[i,0] < 0)):
			lh_cam_2.append([lighthouse_cams[i][3], lighthouse_cams[i][4]])
			gnd_cam.append(gnd_point)

	mat = np.array([[1,0,0,0],
				[0,1,0,0],
				[0,0,1,0],
				[0,0,0,0]])

	#print(np.array([gnd_cam]))
	#print(np.array([lh_cam_1]).astype('float32'))
	K = np.array([[1,0,0],
				 [0,1,0],
				 [0,0,1] ],dtype='float32')

	retval, rvec, tvec = opencv.solvePnP(np.array([gnd_cam]),np.array([lh_cam_2]).astype('float32'),K,None)
	print("pnp lh 2")
	print(retval)
	print(rvec)
	print(tvec)

	print("lh2 rotation mat")
	rot, jacobian = opencv.Rodrigues(rvec)
	print(rot)
	print(np.linalg.inv(rot) @ tvec)
	plt.figure()
	plt.subplot(2,1,1)
	plt.plot(gnd_cam)
	plt.subplot(2,1,2)
	plt.plot(lh_cam_2)


	lh2_obj = Lighthouse()
	lh2_obj.setOpenCVParams(tvec = tvec,rvec = rvec,K = K)

	imagePoints, jacob = opencv.projectPoints(np.array([gnd_cam]),lh2_obj.cvRvec,lh2_obj.cvTvec,lh2_obj.cvK,None)



	cvproj2 = np.squeeze(np.array(imagePoints))

	plt.figure()
	plt.subplot(2,1,1)
	plt.plot(np.array(lh_cam_2)[:,0])
	plt.title("Lighthouse 2 projection check")
	plt.plot(cvproj2[:,0])
	plt.legend(['Lighthouse','OpenCV Projectect Truth'])

	plt.subplot(2,1,2)
	plt.plot(np.array(lh_cam_2)[:,1])
	plt.plot(cvproj2[:,1])
	plt.legend(['Lighthouse','OpenCV Projectect Truth'])



	return lh1_obj, lh2_obj, cvproj1, cvproj2

	

if __name__ == "__main__":

	
	print("Loading Experiment: ", sys.argv[1])

	optitrack_filename = sys.argv[1]+"_optitrack.csv"
	lighthouse_filename = sys.argv[1] + "_lighthouse.csv"
	lh_location_filename = sys.argv[1] + "_lighthouse_locations.txt"
	#load ground truth file
	with open(optitrack_filename) as file:
		time_name = ['Unnamed: 1_level_0', 'Unnamed: 1_level_1', 'Time (Seconds)']

		ground_truth_df = pd.read_csv(file, header = [2,4,5])

		scum_gnd_df = ground_truth_df['scum']
		gnd_time_df = ground_truth_df['Unnamed: 1_level_0', 'Unnamed: 1_level_1', 'Time (Seconds)']
		scum_gnd_df['Time','Time'] = gnd_time_df
		#print(scum_gnd_df)


	#load lighthouse data file
	with open(lighthouse_filename) as file:
		if sys.argv[1] == 'take1' or sys.argv[1] == 'take5':
			lighthouse_df = pd.read_csv(file, header=None, names = ['azA', 'elA', 'azB', 'elB', 'timestamp_10.82 Hz'], dtype = float)
		else: 
			lighthouse_df = pd.read_csv(file, header=None,lineterminator = ']', sep = ',', names = ['azA', 'elA', 'azB', 'elB', 'timestamp_10.82 Hz'],dtype = float)

		lighthouse_df['timestamp_seconds'] = lighthouse_df['timestamp_10.82 Hz'] / 10.82e6
		lighthouse_df['timestamp_optitrack'] = align_data(lighthouse_df)

	with open(lh_location_filename) as file:
		lighthouse_locations = pd.read_csv(file, index_col = "Base Station",)

		board_rotation = np.array([ [-1, 0,0],
									[0, 0, 1], 
									[0, 1, 0] ])
		#multiply by calibration factor
		theta = 20 / 180 * np.pi
		board_rotation = np.dot(np.array([[np.cos(theta),0,np.sin(theta)],
											[0, 1 ,0],
											[-np.sin(theta), 0, np.cos(theta)]]),
											board_rotation)

		lh1,lh2 = process_lh_locs(lighthouse_locations,board_rotation)
		#print(lh1)
		#print(lh1[0][0],lh1[0][1],lh1[0][2],lh1[1])
		lh1_obj = Lighthouse(lh1[1][0],lh1[1][1],lh1[1][2],np.linalg.inv(lh1[0]))
		lh2_obj = Lighthouse(lh2[1][0],lh2[1][1],lh2[1][2],np.linalg.inv(lh2[0]))

		#print(lh1)
		#print(lh2)

	#plot_truth_traj(scum_gnd_df,lh1,lh2)

	print(lh1[1])
	azimuth1_gnd, elevation1_gnd, azimuth2_gnd, elevation2_gnd = project_truth_to_az_el(scum_df = scum_gnd_df,
																						lh1 = lh1[1],
																						lh2 = lh2[1],
																						lh1_to_gnd = lh1[0], 
																						lh2_to_gnd = lh2[0])
	
	plot_raw_data(azimuth1_gnd, azimuth2_gnd, elevation1_gnd, elevation2_gnd, lighthouse_df, scum_gnd_df)
	#plt.show()

	#interpolate lighthouse data times to get ground truth at those data times 
	lighthouse_processed_dict = generate_data_dict(azimuth1_gnd, azimuth2_gnd, elevation1_gnd, elevation2_gnd, lighthouse_df, scum_gnd_df)
	
	trajectory,cam_points = triangulate_scum_nodict(lighthouse_df, lh1_obj, lh2_obj)
	lh1_calibrated, lh2_calibrated, lh1_gnd_proj, lh2_gnd_proj = calibrate(scum_gnd_df,lh1_obj.P,lh2_obj.P,cam_points)

	##############################
	print(np.squeeze(lh1_calibrated.cvTvec))
	azimuth1_gnd, elevation1_gnd, azimuth2_gnd, elevation2_gnd = project_truth_to_az_el(scum_df = scum_gnd_df,
																						lh1 = -np.linalg.inv(lh1_calibrated.pose) @ np.squeeze(lh1_calibrated.cvTvec),
																						lh2 = -np.linalg.inv(lh2_calibrated.pose) @ np.squeeze(lh2_calibrated.cvTvec),
																						lh1_to_gnd = np.linalg.inv(lh1_calibrated.pose), 
																						lh2_to_gnd = np.linalg.inv(lh2_calibrated.pose))


	#interpolate lighthouse data times to get ground truth at those data times 
	lighthouse_processed_dict = generate_data_dict(azimuth1_gnd, azimuth2_gnd, elevation1_gnd, elevation2_gnd, lighthouse_df, scum_gnd_df)
	##############################


	truth = {'az1' : azimuth1_gnd + np.pi/2,
			 'el1' : elevation1_gnd + np.pi/2,
			 'az2' : azimuth2_gnd + np.pi/2,
			 'el2' : elevation2_gnd + np.pi/2,
			 'time' : scum_gnd_df['Time','Time']}

	gnd_traj = triangulate_gnd(truth,lh1_obj, lh2_obj)
	#print(gnd_traj)
	#truth_processed_dict = generate_data_dict(azimuth1_gnd, azimuth2_gnd, elevation1_gnd, elevation2_gnd, truth, scum_gnd_df)
	plt.figure()
	for i in range(1,4):
		
		plt.subplot(3,1,i)
		plt.scatter(gnd_traj[:,0],gnd_traj[:,i],s = 1)
		plt.title('Truth: '+str(i))
		#plt.ylim([0,5])


	#trajectory,cam_points = triangulate_scum_nodict(lighthouse_df, lh1_obj, lh2_obj)
	#lh1_calibrated, lh2_calibrated, lh1_gnd_proj, lh2_gnd_proj = calibrate(scum_gnd_df,lh1_obj.P,lh2_obj.P,cam_points)

	trajectory,cam_points = triangulate_scum_nodict(lighthouse_df,lh1_calibrated, lh2_calibrated)


	plt.figure()
	xyz_map = ['X','Y','Z']
	for i in range(1,4):
		plt.subplot(3,1,i)
		plt.scatter(trajectory[:,0],trajectory[:,i],s = 1)
		plt.scatter(scum_gnd_df['Time', 'Time'],scum_gnd_df['Position',xyz_map[i-1]],s=1)
		plt.xlim([0,100])
		plt.ylabel(xyz_map[i-1] + ' (m)')
		plt.legend(['Scum Lighthouse', 'Ground Truth'])

	plt.xlabel('Time (s)')


	plt.figure()
	plt.plot(truth['az1'])
	plt.title('truth')

	plot_az_el(lighthouse_processed_dict)
	plt.plot(lh1_gnd_proj[:,0])
	plt.figure()
	plt.plot(np.tan(lh1_gnd_proj[:,0])*180/np.pi + 90)

	plt.figure()
	plt.plot(np.tan(lh1_gnd_proj[:,1])*180/np.pi + 90)

	plt.figure()

	plt.plot(np.tan(lh2_gnd_proj[:,0])*180/np.pi + 90)

	plt.figure()
	plt.plot(np.tan(lh2_gnd_proj[:,1])*180/np.pi + 90)

	print(lh1_obj.pose)
	print(lh1_obj.translation)
	print(lh1_calibrated.pose)
	print(lh1_calibrated.translation)

	plt.show()	
