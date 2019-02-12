#!/anaconda3/bin/pythonw

######################################
## Radial Distribution Function ##
## Naeyma Islam ##
## UNO, 2nd Mar 2018 ##
######################################

#----- Importing the libraries -----#

import numpy as np
import pandas as pd
import math as mth
import matplotlib.pyplot as plt
from scipy.interpolate import spline

input_file_name = "pmaaNa.lammpstrj"  # .lammpstrj file name along with location path
out_file_name = "rdf.out"   # output file name along with path to store rgyr values

print ("Start processing Trajectory file....\n")

#----- Reading .lammpstrj file (from lammps output) -----#

cols_read = ['TIMESTEP', 'ATOM_ID', 'ATOM_TYPE', 'x', 'y', 'z', 'xbox', 'ybox', 'zbox']
input_data = pd.DataFrame(columns=cols_read)    # Dataframe for the Timestep data
input_data_new = pd.DataFrame(columns=cols_read)
input_one_row = []

atom_type_o6 = 6  # Oxygen
atom_type_o7 = 7  # Oxygen
atom_type_na = 8  # Na

dr = 0.1

xbox = 0
ybox = 0
zbox = 0
vol_box_sum = 0
count_o6 = 0
count_o7 = 0
count_na = 0
time_step = -1
frame_count = 0

no_bins = 2000                   # assuming total number of bins
g_o6na = []
g_o7na = []
for i in range (0, no_bins):     # initializing bins with '0' value
    g_o6na.insert(i, 0)
    g_o7na.insert(i, 0)

#cols_output = ['R', 'RDF_O6', 'SUM_O6', 'RDF_O7', 'SUM_O7s']
cols_output = ['R', 'AVG_RDF', 'AVG_SUM']
output_data = pd.DataFrame(columns=cols_output)  # Dataframe to store the calculated Rgyr values along with Time
new_output = []

try:
    input_file = open(input_file_name, "r")  #opens the .data file in the working directory in 'reading mode'
    out_file = open (out_file_name, "w")   # opening the out file to write rgyr values
    
    no_of_atoms = 0

    for line in input_file:          
        #print (line)                          # line points to the first line of the FRAME ('ITEM: TIMESTEP')
        next_line = input_file.readline()      # pointing to the 'ITEM: TIMESTEP' line 
        #print (next_line)
        line_as_list = next_line.split()
        time_step = int(line_as_list[0])       # value of the TIMESTEP/FRAME
        print ("->Processing Frame: %d" % time_step)
        frame_count = frame_count + 1          # Increasing frame counter
        
        text_line = input_file.readline()      # pointing to the 'ITEM: NUMBER OF ATOMS' line
        next_line = input_file.readline()      # pointing to the next line of 'ITEM: NUMBER OF ATOMS'
        #print (next_line)
        line_as_list = next_line.split()
        no_of_atoms = int(line_as_list[0])     # number of atoms per frame
        
        text_line = input_file.readline()      # pointing to the 'ITEM: BOX BOUNDS' line
       # print (text_line)
        
        next_line = input_file.readline()      # pointing to the line having X low and high for the box
        #print (next_line)
        line_as_list = next_line.split()
        xbox = float(line_as_list[1]) - float(line_as_list[0])  # X distance for Box
        
        next_line = input_file.readline()      # pointing to the line having Y low and high for the box
        line_as_list = next_line.split()
        ybox = float(line_as_list[1]) - float(line_as_list[0])  # Y distance for Box
        
        next_line = input_file.readline()      # pointing to the line having Z low and high for the box
        line_as_list = next_line.split()
        zbox = float(line_as_list[1]) - float(line_as_list[0])  # Z distance for Box
        
        vol_box_sum = vol_box_sum + (xbox * ybox * zbox)  # adding-up Volume of the Box
        
        text_line = input_file.readline()      # pointing to the 'ITEM: ATOMS' line
            # reading all atoms in the frame
        for i in range (0, no_of_atoms):
            next_line = input_file.readline()  # pointing to the line of each Atoms
            #print (next_line)
            line_as_list = next_line.split()
            if line_as_list and line_as_list[0].isdigit():   # check whether the line is not empty and the first word is a number
                atom_type = int(line_as_list[1])
                if (atom_type == atom_type_o6) or (atom_type == atom_type_o7 or (atom_type == atom_type_na)):      # Process all atoms except Atom_Type_to_exclude (water)

                    x = float(line_as_list[2])   
                    y = float(line_as_list[3])  
                    z = float(line_as_list[4])  
                    
                    input_one_row = [time_step, int(line_as_list[0]), int(line_as_list[1]), x, y, z, xbox, ybox, zbox]   # loading one row data
                    new_row_df = pd.DataFrame ([input_one_row], columns=cols_read, index=[int(line_as_list[0])])          # making dataframe for loaded row
                    input_data = input_data.append(new_row_df)                           # add new row into the Dataframe
                    input_one_row = []                                                 # resetting row placeholder
        # reading of atoms completed
          
            # Check for Boundary condition and update the coordinates accordingly            
        #x_new = []
        #y_new = []
        #z_new = []
        #x_new.insert(0, input_data.iloc[0,3])
        #y_new.insert(0, input_data.iloc[0,4])
        #z_new.insert(0, input_data.iloc[0,5])

        #new_row = [input_data.iloc[0,0], input_data.iloc[0,1], input_data.iloc[0,2], x_new[0], y_new[0], z_new[0], input_data.iloc[0,6], input_data.iloc[0,7], input_data.iloc[0,8]]   # loading one row data
        #new_df = pd.DataFrame ([new_row], columns=cols_read)     # making dataframe for loaded row
        #input_data_new = input_data_new.append(new_df)           # add new row into the Dataframe

        #j = 1    # starting from second atom ('0' is the first one)
        #for i in range(int(input_data.shape[0] - 1)):

            #xj = input_data.iloc[j,3]
            #new_xj = xj - xbox * round ((xj - x_new[j-1])/xbox)
            #x_new.insert(j, new_xj)

            #yj = input_data.iloc[j,4]
            #new_yj = yj - ybox * round ((yj - y_new[j-1])/ybox)
            #y_new.insert(j, new_yj)

            #zj = input_data.iloc[j,5]
            #new_zj = zj - zbox * round ((zj - z_new[j-1])/zbox)
            #z_new.insert(j, new_zj)

            #new_row = [input_data.iloc[j,0], input_data.iloc[j,1], input_data.iloc[j,2], x_new[j], y_new[j], z_new[j], input_data.iloc[j,6], input_data.iloc[j,7], input_data.iloc[j,8]]   # loading one row data
            #new_df = pd.DataFrame ([new_row], columns=cols_read)      # making dataframe for loaded row
            #input_data_new = input_data_new.append(new_df)            # add new row into the Dataframe

            #j = j + 1

        # Boundary checking and coordinate update completed here

            # Creating sub dataframes based on Atom-Type
        df_o6 = input_data.loc[input_data['ATOM_TYPE'] == atom_type_o6]
        df_o7 = input_data.loc[input_data['ATOM_TYPE'] == atom_type_o7]
        df_na = input_data.loc[input_data['ATOM_TYPE'] == atom_type_na]

        count_o6 = count_o6 + int(df_o6.shape[0])       # number of Oxygens of Type-6
        count_o7 = count_o7 + int(df_o7.shape[0])       # number of Oxygens of Type-7
        count_na = count_na + int(df_na.shape[0])       # number of Na of Type-8
            
        for i in range(int(df_o6.shape[0] - 1)):          # for each Oxygen of Type-6
            r_dist = 0
            nbin = 0
            for j in range(int(df_na.shape[0] - 1)):  # for each Na
                
                dx = df_o6.iloc[i,3] - df_na.iloc[j,3]
                dy = df_o6.iloc[i,4] - df_na.iloc[j,4]
                dz = df_o6.iloc[i,5] - df_na.iloc[j,5]
                dx = dx - (xbox * round (dx/xbox))
                dy = dy - (ybox * round (dy/ybox))
                dz = dz - (zbox * round (dz/zbox))
                
                r_dist = mth.sqrt (mth.pow(dx,2) + mth.pow(dy,2) + mth.pow(dz,2))
                nbin = int(r_dist/dr)
                if (nbin >= 0) and (nbin <=  no_bins):
                    g_o6na[nbin] =  g_o6na[nbin] + 1


                #x_dist2 = mth.pow((df_o6.iloc[i,3] - df_na.iloc[j,3]),2)
                #y_dist2 = mth.pow((df_o6.iloc[i,4] - df_na.iloc[j,4]),2)
                #z_dist2 = mth.pow((df_o6.iloc[i,5] - df_na.iloc[j,5]),2)

                #r = mth.sqrt (x_dist2 + y_dist2 + z_dist2)
                #nbin = int(r/dr)
                #if (nbin >= 0) and (nbin <=  no_bins):
                    #g_o6na[nbin] =  g_o6na[nbin] + 1
    
        for i in range(int(df_o7.shape[0] - 1)):          # for each Oxygen of Type-7
            r_dist = 0
            nbin = 0
            for j in range(int(df_na.shape[0] - 1)):  # for each Na
                
                dx = df_o7.iloc[i,3] - df_na.iloc[j,3]
                dy = df_o7.iloc[i,4] - df_na.iloc[j,4]
                dz = df_o7.iloc[i,5] - df_na.iloc[j,5]
                dx = dx - (xbox * round (dx/xbox))
                dy = dy - (ybox * round (dy/ybox))
                dz = dz - (zbox * round (dz/zbox))
                
                r_dist = mth.sqrt (mth.pow(dx,2) + mth.pow(dy,2) + mth.pow(dz,2))
                #print(r_dist)
                nbin = int(r_dist/dr)
                if (nbin >= 0) and (nbin <=  no_bins):
                    g_o7na[nbin] =  g_o7na[nbin] + 1

#r = 0
            #nbin = 0
            #for j in range(int(df_na.shape[0] - 1)):      # for each Na
                #x_dist2 = mth.pow((df_o7.iloc[i,3] - df_na.iloc[j,3]),2)
                #y_dist2 = mth.pow((df_o7.iloc[i,4] - df_na.iloc[j,4]),2)
                #z_dist2 = mth.pow((df_o7.iloc[i,5] - df_na.iloc[j,5]),2)
                
                #r = mth.sqrt (x_dist2 + y_dist2 + z_dist2)
                #nbin = int(r/dr)
                #if (nbin >= 0) and (nbin <=  no_bins):
                    #g_o7na[nbin] =  g_o7na[nbin] + 1
        
        # reset all variables and get ready to process next Frame/Timestep
        no_of_atoms = 0
        xbox = 0
        ybox = 0
        zbox = 0
        input_data = input_data.iloc[0:0]
        input_data_new = input_data_new.iloc[0:0]
        df_o6 = df_o6.iloc[0:0]
        df_o7 = df_o7.iloc[0:0]
        df_na = df_na.iloc[0:0]
            
    # File reading completed here (main for loop ends here)
    print ("\nProcessing of Trajectory file completed.\n")
    
    print ("\nCalculating  th Radial Distribution Function.\n")
        # Now, calculating RDF
    vol_box_avg =  vol_box_sum / frame_count    # average volume of the box
    avg_o6 = count_o6/frame_count
    avg_o7 = count_o7/frame_count
    avg_na = count_na/frame_count
    sum_o6 = 0
    sum_o7 = 0
    #print(count_o6,frame_count,avg_o6,count_na,avg_na)
    norm_o6 = vol_box_avg/(avg_o6 * avg_na)
    norm_o7 = vol_box_avg/(avg_o7 * avg_na)
    
    for i in range (1, 200):        
        vol_temp = 4 * 3.14 * mth.pow(i,2) * mth.pow(dr, 3)
        rdf_o6 = norm_o6 * (g_o6na[i]/(vol_temp * frame_count))
        rdf_o7 = norm_o7 * (g_o7na[i]/(vol_temp * frame_count))
        sum_o6 = sum_o6 + (rdf_o6 * vol_temp/norm_o6/avg_o6)
        sum_o7 = sum_o7 + (rdf_o7 * vol_temp/norm_o7/avg_o7)
        avg_rdf = (rdf_o6 + rdf_o7)/2
        avg_sum = (sum_o6 + sum_o7)/2
            # Writing to file
        #out_file.write(str(round((i * dr),3)) + "\t" +  str(round(rdf_o6,3)) + "\t" + str(round(sum_o6,3)) + "\t" +  str(round(rdf_o7,3)) + "\t" + str(round(sum_o7,3)) + "\n")      # saving rgyr value into the file

        out_file.write(str(round((i * dr),3)) + "\t" +  str(round(avg_rdf,3)) + "\t" + str(round(avg_sum,3)) +"\n")
        
        new_output = pd.DataFrame ([[(i*dr), avg_rdf, avg_sum]], columns=cols_output)
        output_data = output_data.append(new_output)
       
        #print (new_output)
        
    # Processing completed here
    
    #print (output_data)

    out_file.close()                    # closing the out file
    input_file.close()                    # closing the opend trj file
    
except ValueError:
    print ("Cannot find the file or wrong file format. Please check and try again.")
    system.exit(os.EX_CONFIG)

#----- Graph plotting -----#

print ("Ploting Radial Distribution Function ....\n")

plt.plot(output_data['R'], output_data['AVG_RDF'], color='g')
plt.plot(output_data['R'], output_data['AVG_SUM'], color='b')

plt.ylabel('G(R)')
plt.xlabel('R')
plt.plot
plt.show()
