% test the GPS deduplicate script

ingps_x=[1,2,3,4,3];
ingps_y=[0,0,0,0,0];
ingps_vx=[1,1,1,1,2];
ingps_vy=[1,1,1,1,1];
ingps_sx=[1,1,1,1,.1];
ingps_sy=[1,1,1,1,1];

ingps=[ingps_x', ingps_y', ingps_vx', ingps_vy', ingps_sx', ingps_sy']

outgps=deuplicate_GSRM_gps_data(ingps)
