# -*- coding: utf-8 -*-

# Copyright (C) 2021 Gaspar Carrasco-Huertas
# gasparcarrascohuertas@gmail.com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import warnings

import matplotlib.pyplot as plt # import matplotlib library
import matplotlib.gridspec as gridspec # import matplotlib library
import matplotlib as mpl # import matplotlib library
import numpy as np # import numpy library
import glob # import glob library
import sys # import sys library
import os  # import os library
import shutil # import shutil library
import math # import math library
from astropy.table import Table, Column # import astropy library
from astropy.io import ascii # import astropy library
import matplotlib._layoutbox as layoutbox # import numpy library
from mpl_toolkits.axes_grid1.inset_locator import inset_axes # import matplotlib library
import matplotlib.ticker as ptick # import matplotlib library

np.warnings.filterwarnings('ignore')

column_0 =  0 # Start potential
column_1 =  1 # Upper vertex potential
column_2 =  2 # Lower vertex potential
column_3 =  3 # Stop potential
column_4 =  4 # Scan rate
column_5 =  5 # Step
column_6 =  6 # Scan
column_7 =  7 # Scan sourcer windower
column_8 =  8 # Columna 9 del .txt
column_9 =  9 # Potential applied (V)
column_10 = 10 # Time (s)
column_11 = 11 # WE(1).Current (A)
column_12 = 12 # WE(1).Potential (V)
column_13 = 13 # Index

SIZE_Font = 10



#-------------------------------------PREAMBLE--------------------------------------------------
#load files
files = glob.glob('*.txt') #load every file with .txt extension 
print(files) #print in screen "files loaded"
#Raw text remarks
input_delimiter= "\t" #de raw data delimiter is tabular
input_skip_space=2 #we want 2 skip space from raw data
#Path in which we are working
path=os.getcwd() #get path 

#---------------------------------------------------------------------------------------
def new_directories(folder1,folder2):

    #This  function create new folders named as folder1 and folder2 in where will be stored figures and data separately and has been created 
    #by Gaspar Carrasco. gasparcarrascohuertas@gmail.com for contact


    folder_1= folder1 #folder 1
    folder_2= folder2 #folder 2
    path= os.getcwd() + "\\" #obtain path
    dir1 = path+folder_1  #'path_to_my_folder_1'
    dir2 = path+folder_2  #'path_to_my_folder_2'

    print(dir1)
    print(dir2)

    if not os.path.exists(dir1 and dir2): # if the directory does not exist
        os.makedirs(dir1 ) # make the directory 1
        os.makedirs(dir2)  # make the directory 2
    else: # the directory exists
        #removes all files in a folder
        for the_file in os.listdir(dir1 and dir2):
            file_path = os.path.join(dir1, the_file)
            file_path = os.path.join(dir2, the_file)
    print("-------------------------------------END CREATE NEW DIRECTORIES FUNCTION---------------------------------------------------")
#---------------------------------------------------------------------------------------
def removing_files_in_folder(folder1,folder2):

    #This  function remove every file in new folders named as folder1 and folder2  and has been created 
    #by Gaspar Carrasco. gasparcarrascohuertas@gmail.com for contact

    path= os.getcwd() + "\\" #Obtain path 
    folder_1= folder1 #folder 1
    folder_2= folder2 #folder 2

    list_files = os. listdir(path+folder1)  #list all the files in folder 1

    files_figures = glob.glob(path+folder1+ "/*") 
    for f in files_figures:  #for every file  file in folder 2 remove
     os.remove(f)  

    files_data = glob.glob(path+folder2+ "/*")
    for f in files_data:   #for every file  file in folder 2 remove
     os.remove(f) 

    print("-------------------------------------END REMOVING FILES IN FOLDERS FUNCTION---------------------------------------------------")
#---------------------------------------------------------------------------------------
def file_change_name(imput):

    #This  function change the name of O2 files in order to obtain ordered list in our directory and has been created by Gaspar Carrasco. gasparcarrascohuertas@gmail.com for contact

    for name_old in imput: #for every file in our O2 list 

        a=name_old[5:-6] # characters located between 3 and -7 position
        b = float(a) # convert characters selected to float
        c = '%04i' % b # add 0 up to 4 positions before numbers converted to float 
        print("The file name with 0 placed is: " + c)
        name_new= name_old.replace(a, c)
        print("Renamed file is: " + name_new)
        print(name_old, ' ------- ',   name_new)
        os.rename(name_old,name_new)
    print("-------------------------------------END CHANGE NAME OF OXYGEN FILES FUNCTION---------------------------------------------------")
#---------------------------------------------------------------------------------------
def electrochemical_parameters(imput):
    """This function obtain electrochemical info of input analysis files performed in ORR experiments. Colum association must be set in equipment previously   
    has been created by Gaspar, gasparcarrascohuertas@gmail.com for contact"""
    
    cycles= float(input("Introduce numer of scan to perform: ")) 

    for file in imput:

        f = np.genfromtxt(file, delimiter=input_delimiter, skip_header=input_skip_space)

        #------------------
        #Electrochemical info 

        start_potential = f[:1, column_0] #Define our column and position associated to start potential selected in equipment
        upper_vertex_potential = f[:1, column_1]  #Define our column and position associated to upper vertex potential  selected in equipment
        stop_potential = f[:1, column_3]  #Define our column and position associated to stop potential selected in equipment
        lower_vertex_potential = f[:1, column_2]  #Define our column and position associated to lower vertex potential selected in equipment
        scan_rate = f[:1, column_4]  #Define our column and position associated to scan rate selected in equipment
        step = f[:1, column_5]  #Define our column and position associated to analysis step selected in equipment
        cycle_selected= f[:1, column_7] #Define our column and position associated to start cycle we have selected in equipment

        print("\n")
        print(str(file))
        print("Start potential is : " + str(start_potential) + "V")
        print("Upper_vertex_potential is: " + str(upper_vertex_potential) + "V")
        print("Stop potential is: " + str(stop_potential) + "V")
        print("Lower vertex potential is: " + str(lower_vertex_potential) + "V")
        print("Scan rate is: " + str(scan_rate)  + "V/seg")
        print("Step is: " + str(step)  + "V")
        print("Scan selected is : " + str(cycle_selected) )
        print("\n")
     
        #------------------
        #Operations with upper & lower vertex potential 
        rango=upper_vertex_potential-lower_vertex_potential
        print("Potential range scaned for 1 branch is : " + str(rango)+ "V")
        total_range=rango*2
        print("Potential range scaned for 2 branch is : " + str(total_range)+ "V")

        #------------------
        brach_values=rango
        step_values=brach_values/step
        print("Step values for one brach are : " + str(step_values) + " positions") # This value is the number of calculated positions (a.u.) in one branch based on step value selected
        print("Step values for both branches are : " + str(step_values*2) + " positions") # This value is the number of calculated positions (a.u.) in both branches based on step value selected

        #------------------
        #Analysis time  (seg)
        #------------------
        interval_scan_time=step/scan_rate
        print("Interval time : " + str(interval_scan_time)+ "seg")
        branch_time=rango/scan_rate
        print("Time spend in measuring one brach is : " + str(branch_time) + "seg")
        time_both_braches=branch_time*2
        print("Scan time  " + str(cycle_selected) + " is : " + str(time_both_braches) + "seg")
        time_analysis=(branch_time*2)*cycles
        print("Time used for  " + str(cycles) + " is : " + str(time_analysis) + "seg")
        print("\n")

    print("-------------------------------------END ELECTROCHEMICAL PARAMETERS FUNCTION---------------------------------------------------")
#---------------------------------------------------------------------------------------
def ECSA():

    #This  function plot all the scan rates performed in electrochemical surface area analysis  and has been created 
    #by Gaspar Carrasco. gasparcarrascohuertas@gmail.com for contact


    counter=0 #counter set to 0
    label = ['5 mV/seg',"10 mV/seg",'20 mV/seg',"40 mV/seg","60 mV/seg","80 mV/seg","100 mV/seg","125 mV/seg","150 mV/seg","200 mV/seg"] #label set for graph
    color = ['black','red','blue',"green","orange","purple","yellow","gray","pink","brown"] #colors set for graphs
    fig, ax = plt.subplots()
    for file in files:   ## For each file in files

        f = np.genfromtxt(file, delimiter="\t", skip_header=1) #charge the file with tabular delimiter and skip header of 1 line
        WE_potential = f[:, column_9] #Load from the "file" the column asociated to WE potential (V)
        WE_current = f[:, column_11]  #Load from the "file" the column asociated to WE current (A)
        WE_current_corrected=WE_current*1000000 #change WE current to microamperes

        

        ax.plot(WE_potential,WE_current_corrected,color=color[counter],label=label[counter]) #Plot the "WE potential (V)" vs  "WE current (A)"
        counter=counter+1 #we add 1 to counter for colors and labels

    #--------------------PLOTS OPTIONS--------------------------


    ax.tick_params(axis="y", right=True, direction='in')
    ax.tick_params(axis="x",top=True , direction='in')
    ax.legend(loc='lower right', prop={'size':8}) #graph legend

    ax.set_xlabel('Potencial vs. Ag/AgCl (V)')
    ax.set_ylabel('Intensidad de corriente (\u00B5A) ')

    #plt.title("Cyclic voltammetry \n Electrode Area")# graph title
    ax.set_xlim(-0.1,0.8)
    ax.set_ylim(-60,60)  # Y axe limits
    #plt.grid()# paint a grid over the graph
    ax.figure.savefig("figure_ECSA_Dropsens.png")
    ax.figure.savefig("figure_ECSA_Dropsens.eps")



    shutil.move("figure_ECSA_Dropsens.png",  "figures") #Move the graph as "name.png" to folder figures
    shutil.move("figure_ECSA_Dropsens.eps",  "figures") #Move the graph as "name.png" to folder figures
    #plt.show() # Show graph 
    #ax.clf() #clear the figure - you can still paint another plot onto it
    print("-------------------END ECSA REDUCION BRANCH FUNCTION---------------------------")
#---------------------------------------------------------------------------------------
def ECSA_reduction_branch():

    #This  function plot all the reduction branch of scan rates performed in electrochemical surface area analysis  and has been created 
    #by Gaspar Carrasco. gasparcarrascohuertas@gmail.com for contact

    counter=0 #counter set to 0
    label = ['5 mV/seg',"10 mV/seg",'20 mV/seg',"40 mV/seg","60 mV/seg","80 mV/seg","100 mV/seg","125 mV/seg","125 mV/seg","200 mV/seg"] #label set for graph
    color = ['black','red','blue',"green","orange","purple","yellow","gray","pink","brown"] #colors set for graphs

    data_list_1=[] #empty data list
    data_list_2=[] #empty data list
    data_list_3=[] #empty data list
    data_list_4=[] #empty data list

    for file in files:   ## For each file in files

        f = np.genfromtxt(file, delimiter="\t", skip_header=1) #charge the file with tabular delimiter and skip header of 1 line
        WE_potential = f[:623, column_9] #Load from the file the column asociated to WE potential (V)
        WE_current = f[:623, column_11] #Load from the file the column asociated to WE current (A)
        WE_current_corrected=WE_current*1000000 #change WE current to microamperes
        index=  f[:623, column_13] #Load from the file the column asociated to index (a.u.)
        
        #--------------------LINEAR REGRESSION--------------------------
        #we need to define two analysis ranges for the regression
        #---------FIRST-----------
        time1=1625 #Start interval
        time2=1951 #End interval
        i_interval_1 = np.where( (index< time2) & (index > time1) )[0] # defining interval
        x_interval_1 = WE_potential[i_interval_1] #Find in my interval values asociated to x (WE potential (V))
        y_interval_1 = WE_current_corrected[i_interval_1]#Find in my interval values asociated to y (WE current (microA))
        #---------SECOND-----------
        tiempo3=1700 #Start interval
        tiempo4=1710 #End interval
        i_interval_2 = np.where( (index < tiempo4) & (index > tiempo3) )[0] # defining interval
        x_interval_2 = WE_potential[i_interval_2]#Find in my interval values asociated to x (WE potential (V))
        y_interval_2 = WE_current_corrected[i_interval_2]#Find in my interval values asociated to y (WE current (microA))
        #---------LINEAR REGRESSION-----------
        adjust = np.polyfit(x_interval_2, y_interval_2, deg=1) # linealice the x-values and y-values of the interval to polynomial function order 1
        y_adjust = np.polyval(adjust, x_interval_1) #applying the polynomial function which you got using polyfit
        #print("x-value and y-value fits the linear regression as A and B= "+  str(adjust))
        print(file)
        data_list_1.append(adjust)

        #--------------------FIND LOCAL MINIMUM--------------------------
        #we need to define analysis range in which minimum is located 
        x1_1=(-0.1) #Start interval
        x1_2=(0.4) #End interval
        i_interval_min = np.where( (WE_potential < x1_2) & (WE_potential > x1_1) )[0] # defining interval
        x_interval_min = WE_potential[i_interval_min] #Find in my interval x- values asociated to WE potential (V)
        y_interval_min = WE_current_corrected[i_interval_min] #Find in my interval y-values asociated to WE current (microA)
        min_y = np.min( y_interval_min )
        index1 = np.where(y_interval_min == min_y)[0]  ## Array in which positions y = ymin
        min_x = x_interval_min[index1][0]   ## xmin es el valor de x al que corresponde el ymin
        text_min = (str(min_x)+ ";"+ str(min_y)+ "\n")
        print("Minimum x-value of function is:", min_x, "with y-value:", min_y)
        data_list_2.append(min_x) 
        data_list_3.append(min_y) 
        
        #--------------------FIND INTERSECTION MINIMUM--------------------------
        #we need to know the intersection between the linear regression and the x-minimum  value
        y_intersec = np.polyval(adjust, min_x)
        print("y-value in intersection is " + str(y_intersec))
        print("----------------------------------------")
        data_list_4.append(y_intersec) 

        #--------------------PLOTS--------------------------
        plt.vlines(min_x, -100, 100,color='k', linestyle='--') # plot vertical lines in x-value for minimum from -100 to 100 color black
        plt.plot(min_x, y_intersec, 'Xg') #plot green point in the x-value for minimum and y-value for intersection 
        plt.plot(x_interval_1, y_adjust,'k',linestyle='--') #plot 
        plt.plot(x_interval_1, y_interval_1,'-',color='k',label='1',linewidth=1)
        plt.plot(WE_potential,WE_current_corrected,color=color[counter],label=label[counter]) #plot x-values for potential(V) and y-values for current(microA) 

        counter=counter+1

    #--------------------TABLES--------------------------
    data = (data_list_2,data_list_3,data_list_4) 
    data_array=np.array(data)
    tabla= Table(data_array.transpose())
    new_column=Column([0.005,0.010,0.020,0.040,0.060,0.080,0.100,0.125,0.150,0.200], name="Scan rate (V/seg)")
    tabla.add_column(new_column)
    tabla.rename_column("col0","WE Potential min (V)") # Rename column 0 to  max_x (V)
    tabla.rename_column("col1","WE Current min (\u00B5A)") # Rename column 1 to  max_y (microA)
    tabla.rename_column("col2","y-axe intersection (\u00B5A)") # Rename column 2 to  y_intersec (microA)
    print(tabla)
    tabla.write('data_info_table_reduction.txt', format='ascii')
    
    #--------------------PLOTS OPTIONS--------------------------
    plt.xlabel('Potencial (V) vs. Ag/AgCl')
    plt.ylabel('Intensidad de corriente (\u00B5A) ')
    plt.title("Cyclic voltammetry \n Cathodic peak current") # graph title
    plt.xlim(-0.1,0.6)  # X axe limits
    plt.ylim(-60,60)  # Y axe limits
    #plt.grid() # paint a grid over the graph
    plt.savefig("figure_reduction_Dropsens.png")# Save the graph as "name.png"
    plt.savefig("figure_reduction_Dropsens.eps")# Save the graph as "name.png"
    shutil.move("figure_reduction_Dropsens.png",  "figures")#Move the graph as "name.png" to folder figures
    shutil.move("figure_reduction_Dropsens.eps",  "figures")#Move the graph as "name.png" to folder figures
    #plt.show() # Show graph 
    plt.clf() #clear the figure - you can still paint another plot onto it
    print("---------------------------------------------------END ECSA REDUCION BRANCH FUNCTION---------------------------------------------------")

def ECSA_oxidation_branch():

    #This  function plot all the oxidation branch of scan rates performed in electrochemical surface area analysis  and has been created 
    #by Gaspar Carrasco. gasparcarrascohuertas@gmail.com for contact
    
    counter=0 #counter set to 0
    label = ['5 mV/seg',"10 mV/seg",'20 mV/seg',"40 mV/seg","60 mV/seg","80 mV/seg","100 mV/seg","125 mV/seg","150 mV/seg","200 mV/seg"] #label set for graph
    color = ['black','red','blue',"green","orange","purple","yellow","gray","pink","brown"] #colors set for graphs

    data_list_1=[] #empty data list
    data_list_2=[] #empty data list
    data_list_3=[] #empty data list
    data_list_4=[] #empty data list

    t=Table() #creates new empty table


    for file in files:   ## For each file in files

        f = np.genfromtxt(file, delimiter="\t", skip_header=1) #charge the file with tabular delimiter and skip header of 1 line
        WE_potential = f[:623, column_9] #Load from the file the column asociated to WE potential (V)
        WE_current = f[:623, column_11] #Load from the file the column asociated to WE current (A)
        WE_current_corrected=WE_current*1000000 #change WE current to microamperes
        index=  f[:623, column_13] #Load from the file the column asociated to index (a.u.)
        
        #--------------------LINEAR REGRESSION--------------------------
        #we need to define two analysis ranges for the regression
        #---------FIRST-----------
        time1=1304 #Start interval
        time2=1600 #End interval
        i_interval_1 = np.where( (index < time2) & (index > time1) )[0] # defining interval
        x_interval_1 = WE_potential[i_interval_1] #Find in my interval values asociated to x (WE potential (V))
        y_interval_1 = WE_current_corrected[i_interval_1] #Find in my interval values asociated to y (WE current (microA))
        #---------SECOND-----------
        #intervalo para regresion
        tiempo3=1400 #Start interval
        tiempo4=1414  #End interval
        i_interval_2 = np.where( (index < tiempo4) & (index > tiempo3) )[0] # defining interval
        x_interval_2 = WE_potential[i_interval_2] #Find in my interval values asociated to x (WE potential (V))
        y_interval_2 = WE_current_corrected[i_interval_2] #Find in my interval values asociated to y (WE current (microA))
        #---------LINEAR REGRESSION-----------
        adjust = np.polyfit(x_interval_2, y_interval_2, deg=1) # linealice the x-values and y-values of the interval to polynomial function order 1
        y_adjust = np.polyval(adjust, x_interval_1) #applying the polynomial function which you got using polyfit
        print(file)
        #print("x-value and y-value fits the linear regression as A and B= "+  str(adjust))
        data_list_1.append(adjust)
        
        #--------------------FIND LOCAL MAXIMUM--------------------------
        #we need to define analysis range in which minimum is located 
        x2_1=(0.0) #Start interval
        x2_2=(0.5) #End interval
        i_interval = np.where( (WE_potential < x2_2) & (WE_potential > x2_1) )[0] # defining interval
        x_interval = WE_potential[i_interval] #Find in my interval values asociated to x (WE potential (V))
        y_interval = WE_current_corrected[i_interval] #Find in my interval values asociated to y (WE current (microA))
        max_y = np.max( y_interval )
        index2 = np.where(y_interval == max_y)[0]  ## Array in which positions y = ymax
        max_x = x_interval[index2][0]  ## xmax es el valor de x al que corresponde el ymax
        #Files with data extracted from x-y-max  .txt
        text_max = (str(max_x)+ ";"+ str(max_y)+ "\n")        
        print("Maximum x-value of function is :", max_x, "with y-value:", max_y)
        data_list_2.append(max_x)   
        data_list_3.append(max_y)      
        
        #--------------------FIND INTERSECTION MAXIMUM--------------------------
        #we need to know the intersection between the linear regression and the x-maximum  value
        y_intersec = np.polyval(adjust, max_x)
        print("y-value in intersection is: " + str(y_intersec))
        print("----------------------------------------")
        data_list_4.append(y_intersec) 

        #--------------------PLOTS--------------------------
        plt.vlines(max_x, -50, 50,color='k', linestyle='--') # plot vertical lines in x-value for minimum from -100 to 100 color black
        plt.plot(max_x, y_intersec, 'Xg')  #plot green point in the x-value for minimum and y-value for intersection 
        plt.plot(x_interval_1, y_adjust,'k',linestyle='--') #plot 
        plt.plot(x_interval_1, y_interval_1,'-',color='k',label='1',linewidth=1)
        plt.plot(WE_potential,WE_current_corrected,color=color[counter],label=label[counter]) #plot x-values for potential(V) and y-values for current(microA) 
        
        counter=counter+1
        
    #--------------------TABLES--------------------------
    data = (data_list_2,data_list_3,data_list_4) #max_x , max_y , y_intersec
    data_array=np.array(data)
    tabla= Table(data_array.transpose())
    new_column=Column([0.005,0.010,0.020,0.040,0.060,0.080,0.100,0.125,0.150,0.200], name="Scan rate (V/seg)")
    tabla.add_column(new_column)
    tabla.rename_column("col0","WE Potential max (V)") # Rename column 0 to  max_x (V)
    tabla.rename_column("col1","WE Current max (\u00B5A)") # Rename column 1 to  max_y (microA)
    tabla.rename_column("col2","y-axe intersection (\u00B5A)") # Rename column 2 to  y_intersec (microA)
    print(tabla)
    tabla.write('data_info_table_oxidation.txt', format='ascii')
    
    #--------------------PLOTS OPTIONS--------------------------
    plt.xlabel('Potencial (V) vs. Ag/AgCl')
    plt.ylabel('Intensidad de corriente (\u00B5A) ')
    plt.title("Cyclic voltammetry \n Anodic peak current") # graph title
    plt.xlim(-0.1,0.6)  # X axe limits
    plt.ylim(-60,60)  # Y axe limits
    #plt.grid() # paint a grid over the graph
    plt.savefig("figure_oxidation_Dropsens.png") # Save the graph as "name.png"
    plt.savefig("figure_oxidation_Dropsens.eps") # Save the graph as "name.png"
    shutil.move("figure_oxidation_Dropsens.png",  "figures") #Move the graph as "name.png" to folder figures
    shutil.move("figure_oxidation_Dropsens.eps",  "figures") #Move the graph as "name.png" to folder figures
    #plt.show() # Show graph 
    plt.clf() #clear the figure - you can still paint another plot onto it
    print("-------------------------------------END ECSA REDUCION BRANCH FUNCTION---------------------------------------------------")
#---------------------------------------------------------------------------------------
def operation():

    f1 = np.genfromtxt("data_info_table_oxidation.txt", delimiter=" ", skip_header=1)
    intersect_oxidation = f1[:, 2]


    intensity_oxidation_peak = f1[:, 1]

    list1=[]
    ipa=intensity_oxidation_peak - intersect_oxidation 

    print("Intensity anodic peak current (micoA):  " +str(ipa))



    f2 = np.genfromtxt("data_info_table_reduction.txt", delimiter=" ", skip_header=1)
    intersect_reduction = f2[:, 2]

    intensity_reduction_peak = f2[:, 1]

    list2=[]
    ipc=intensity_reduction_peak - intersect_reduction 
    #print(ipc)


    data_ipa = (ipa) #
    array_ipa=np.array(data_ipa)
    np.savetxt("data_ipa.txt", array_ipa.transpose())


    data_ipc = (ipc) #
    array_ipc=np.array(data_ipc)
    np.savetxt("data_ipc.txt", array_ipc.transpose())




    print("-------------------------------------END OPERATION FUNCTION---------------------------------------------------")
#---------------------------------------------------------------------------------------
def Randles_Sevick_equation():

    """
    This  function calculate Area from Randle-Sevick equation acording to slope obtained of "ipa vs root square of scan rate" and has been created 
    by Gaspar Carrasco. gasparcarrascohuertas@gmail.com for contact
    """

    n = 1   #number of electrons transferred in the redox event (usually 1)
    #A  == electrode area in cm2
    #F = 96485.33289  #Faraday Constant in C mol−1 == 96485.33289 mol−1
    D = 0.0000076  #diffusion coefficient in cm2/s == A 25°C los valores Do  (coeficiente de difusión) para el K_3 [Fe(CN)_6] en KCl 0,10 M son 7.6∗〖10〗^(−6) cm2 / s, respectivamente
    C = 0.000001  #concentration in mol/cm3 ----0.001 mol/L
    R = 8.3144598   #Gas constant in J K−1 mol−1  ==8.3144598 J⋅mol−1⋅K−1
    T = 298  #temperature in K = 25ºC ==298 K
    #Randles–Sevcik equation (anodic peak)
    #ip=268600*(n**3/2)*(A)*(D**0.5)*(C)*(value**0.5)

    fig, ax1 = plt.subplots()
    fig, ax2 = plt.subplots()
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    list_1 = []
    list_2 = []

    #--------------------------------------CATHODIC PROCESS-------------------------------------

    f_1 = np.genfromtxt("data_ipc.txt")
    ipc = np.array(f_1)
    print("ipc= "+  str(ipc))

    ip_1 = ipc  # Cathodic peak current in microamps (microA)
    ip_corrected_1=ip_1*0.000001 #Cathodic  peak current  in amps (A) 

    list_1=(ipc,ip_corrected_1) #list with anodic peak current in amps and in microamps
    array_1=np.array(list_1) #array of list 
    #print(array_1)
    tabla_1= Table(array_1.transpose()) #table of array
    new_column_1=Column([0.005,0.010,0.020,0.040,0.060,0.080,0.100,0.125,0.150,0.200], name="Scan rate (V/seg)")
    tabla_1.add_column(new_column_1) #add new column 
    tabla_1.write("data_info_table_oxidation_global.txt", format='ascii')

    f_2 = np.genfromtxt("data_info_table_oxidation_global.txt", delimiter=" ", skip_header=1)
    scan_rate_list_1 = f_2[:, 2]  #scan rate in V/s
    print(scan_rate_list_1)

    #-------------------- ip vs v^1/2 LINEAR REGRESSION--------------------------
    #---------FIRST-----------
    w_1=0 #start value of range
    w_2=0.5 #end value of range
    i_interval_1 = np.where( ((scan_rate_list_1**0.5) < w_2) & ((scan_rate_list_1**0.5) > w_1) )[0] #range of regresion
    x_interval_1 = (scan_rate_list_1**0.5)[i_interval_1] #find in X axis the range of regresion
    y_interval_1 = ip_corrected_1[i_interval_1] #find in Y axis the range of regresion

    #---------LINEAR REGRESSION-----------
    adjust_1 = np.polyfit(x_interval_1, y_interval_1, deg=1)
    y_adjust_1 = np.polyval(adjust_1, x_interval_1)
    a_1 = adjust_1[0] #slope
    b_1 = adjust_1[1] #independ variable
    print("x-value and y-value fits the linear regression as A and B= "+  str(adjust_1))


    def lsqfity(X, Y):

        """
        Calculate a "MODEL-1" least squares fit.

        The line is fit by MINIMIZING the residuals in Y only.

        The equation of the line is:     Y = my * X + by.

        Equations are from Bevington & Robinson (1992)
        Data Reduction and Error Analysis for the Physical Sciences, 2nd Ed."
        pp: 104, 108-109, 199.

        Data are input and output as follows:

        my, by, ry, smy, sby = lsqfity(X,Y)
        X     =    x data (vector)
        Y     =    y data (vector)
        my    =    slope
        by    =    y-intercept
        ry    =    correlation coefficient
        smy   =    standard deviation of the slope
        sby   =    standard deviation of the y-intercept

        """

        X, Y = map(np.asanyarray, (X, Y))

        # Determine the size of the vector.
        n = len(X)

        # Calculate the sums.

        Sx = np.sum(X)
        Sy = np.sum(Y)
        Sx2 = np.sum(X ** 2)
        Sxy = np.sum(X * Y)
        Sy2 = np.sum(Y ** 2)

        # Calculate re-used expressions.
        num = n * Sxy - Sx * Sy
        den = n * Sx2 - Sx ** 2

        # Calculate my, by, ry, s2, smy and sby.
        my = num / den
        by = (Sx2 * Sy - Sx * Sxy) / den
        ry = num / (np.sqrt(den) * np.sqrt(n * Sy2 - Sy ** 2))

        diff = Y - by - my * X

        s2 = np.sum(diff * diff) / (n - 2)
        smy = np.sqrt(n * s2 / den)
        sby = np.sqrt(Sx2 * s2 / den)



        return my, by, ry, smy, sby    



    print(lsqfity(x_interval_1,y_interval_1))




    #-------------------- AREA CALCULATING--------------------------
    Area_final_cathodic=(a_1)/((269000)*(1)*(C)*((D**0.5))) #calculate the area
    print("Area (cm2) Cathodic " + str(Area_final_cathodic))



    text_file_1 = open("data_area_cathodic_convolution.txt", "w")
    n_1 = text_file_1.write(str(Area_final_cathodic))
    text_file_1.close()


    #--------------------------------------ANODIC PROCESS-------------------------------------

    f_3 = np.genfromtxt("data_ipa.txt")
    ipa = np.array(f_3)
    print("ipa= "+  str(ipa))
    ip_2 = ipa  # Cathodic peak current in microamps (microA)
    ip_corrected_2=ip_2*0.000001 #Cathodic  peak current  in amps (A) 

    list_2=(ipa,ip_corrected_2) #list with anodic peak current in amps and in microamps
    array_2=np.array(list_2) #array of list 
    #print(array_2)
    tabla_2= Table(array_2.transpose()) #table of array
    new_column_2=Column([0.005,0.010,0.020,0.040,0.060,0.080,0.100,0.125,0.150,0.200], name="Scan rate (V/seg)")
    tabla_2.add_column(new_column_2) #add new column 
    tabla_2.write("data_info_table_reduction_global.txt", format='ascii')

    f_4 = np.genfromtxt("data_info_table_reduction_global.txt", delimiter=" ", skip_header=1)
    scan_rate_list_2 = f_4[:, 2]  #scan rate in V/s
    print(scan_rate_list_2)

    #-------------------- ip vs v^1/2 LINEAR REGRESSION--------------------------
    #---------FIRST-----------
    w_1_2=0 #start value of range
    w_2_2=0.5 #end value of range
    i_interval_2 = np.where( ((scan_rate_list_2**0.5) < w_2_2) & ((scan_rate_list_2**0.5) > w_1_2) )[0] #range of regresion
    x_interval_2 = (scan_rate_list_2**0.5)[i_interval_2] #find in X axis the range of regresion
    y_interval_2 = ip_corrected_2[i_interval_2] #find in Y axis the range of regresion

    #---------LINEAR REGRESSION-----------
    adjust_2 = np.polyfit(x_interval_2, y_interval_2, deg=1)
    y_adjust_2 = np.polyval(adjust_2, x_interval_2)
    a_2 = adjust_2[0] #slope
    b_2 = adjust_2[1] #independ variable
    print("x-value and y-value fits the linear regression as A and B= "+  str(adjust_2))

    def lsqfity(X, Y):

        """
        Calculate a "MODEL-1" least squares fit.

        The line is fit by MINIMIZING the residuals in Y only.

        The equation of the line is:     Y = my * X + by.

        Equations are from Bevington & Robinson (1992)
        Data Reduction and Error Analysis for the Physical Sciences, 2nd Ed."
        pp: 104, 108-109, 199.

        Data are input and output as follows:

        my, by, ry, smy, sby = lsqfity(X,Y)
        X     =    x data (vector)
        Y     =    y data (vector)
        my    =    slope
        by    =    y-intercept
        ry    =    correlation coefficient
        smy   =    standard deviation of the slope
        sby   =    standard deviation of the y-intercept

        """

        X, Y = map(np.asanyarray, (X, Y))

        # Determine the size of the vector.
        n = len(X)

        # Calculate the sums.

        Sx = np.sum(X)
        Sy = np.sum(Y)
        Sx2 = np.sum(X ** 2)
        Sxy = np.sum(X * Y)
        Sy2 = np.sum(Y ** 2)

        # Calculate re-used expressions.
        num = n * Sxy - Sx * Sy
        den = n * Sx2 - Sx ** 2

        # Calculate my, by, ry, s2, smy and sby.
        my = num / den
        by = (Sx2 * Sy - Sx * Sxy) / den
        ry = num / (np.sqrt(den) * np.sqrt(n * Sy2 - Sy ** 2))

        diff = Y - by - my * X

        s2 = np.sum(diff * diff) / (n - 2)
        smy = np.sqrt(n * s2 / den)
        sby = np.sqrt(Sx2 * s2 / den)



        return my, by, ry, smy, sby    



    print(lsqfity(x_interval_2,y_interval_2))


    #-------------------- AREA CALCULATING--------------------------
    Area_final_anodic=abs((a_2)/((269000)*(1)*(C)*((D**0.5)))) #calculate the area
    print("Area (cm2) Anodic " + str(Area_final_anodic))

    text_file_2 = open("data_area_anodic_convolution.txt", "w")
    n_2 = text_file_2.write(str(Area_final_anodic))
    text_file_2.close()


    #--------------------PLOTS--------------------------


    ax1.plot((scan_rate_list_1**0.5),ip_corrected_1,"o", color='black') #"plot root square of scan rate" vs "ipc"
    ax2.plot((scan_rate_list_2**0.5),ip_corrected_2,"o", color='red') #"plot root square of scan rate" vs "ipc"


    #--------------------PLOTS OPTIONS--------------------------
    #plt.title("Randles Sevick plot area ") # graph title
    #plt.ylim(-50,50)   # Y axe limits
    #plt.grid() # paint a grid over the graph
    #plt.show() # Show graph 

    ax1.tick_params(axis="y", right=True, direction='in', labelcolor="black")
    ax1.tick_params(axis="x",top=True , direction='in')
    #ax.legend(loc='lower right', prop={'size':8}) #graph legend
    ax1.set_xlabel(r"Velocidad de barrido $v^\frac{1}{2}$ (V/s) ",  fontsize= SIZE_Font)
    ax1.set_ylabel('Ip cathodic(A)',  fontsize= SIZE_Font,  color="black")
    #plt.title("Cyclic voltammetry \n Electrode Area")# graph title
    ax1.set_xlim(0,0.5)
    #plt.grid()# paint a grid over the graph
    # Change the y ticklabel format to scientific format
    ax1.set_ylim(-0.00005,0.00005)  # Y axe limits
    ax1.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))



    ax2.tick_params(axis="y", right=True, direction='in',  labelcolor="red")
    ax2.tick_params(axis="x",top=True , direction='in')
    #ax.legend(loc='lower right', prop={'size':8}) #graph legend
    ax2.set_xlabel(r"Velocidad de barrido $v^\frac{1}{2}$ (V/s) ",  fontsize= SIZE_Font)
    ax2.set_ylabel('Ip anodic (A)',  fontsize= SIZE_Font,  color="red")
    #plt.title("Cyclic voltammetry \n Electrode Area")# graph title
    ax2.set_xlim(0,0.5)
    ax2.set_ylim(-0.00005,0.00005)  # Y axe limits
    #plt.grid()# paint a grid over the graph
    # Change the y ticklabel format to scientific format
    #ax1.set_ylim(0,0.0005)  # Y axe limits
    ax2.ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))

    fig.tight_layout()  # otherwise the right y-label is slightly clipped



    ax1.figure.savefig("figure_randles_sevick_plot_cathodic.png")
    ax1.figure.savefig("figure_randles_sevick_plot_cathodic.eps")
    shutil.move("figure_randles_sevick_plot_cathodic.png",  "figures") #Move the graph as "name.png" to folder figures
    shutil.move("figure_randles_sevick_plot_cathodic.eps",  "figures") #Move the graph as "name.png" to folder figures




    plt.clf() #clear the figure - you can still paint another plot onto it
    print("-------------------------------------END RANDLES-SEVICK FUNCTION---------------------------------------------------")
#---------------------------------------------------------------------------------------
def move_data():

    """ This function move every data file .txt to folder data and
    has been created by Gaspar, gasparcarrascohuertas@gmail.com for contact
    """

    files = glob.glob(path+ "/data*.txt")
    print(files)
    print("Files in moved to directories are: "+str(files))
    for f in files:
     shutil.move(f, "data")
    print("-------------------------------------END MOVING DATA .TXT  FUNCTION---------------------------------------------------")
#---------------------------------------------------------------------------------------


new_directories("figures", "data")
removing_files_in_folder("figures","data")
file_change_name(files)
electrochemical_parameters(files)
ECSA()
ECSA_reduction_branch()
ECSA_oxidation_branch()
operation()
Randles_Sevick_equation()
move_data()