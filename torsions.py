#Stuff for calculation of SASA by RDKit
import os
import os.path as op
import math
import numpy as np                              #math functions
from operator import itemgetter
import pandas as pd                             #data set handling
import matplotlib.pyplot as plt                 #plotting
# Extracts a set of torsion environments from a torenv dat file from fraggen,
# and stores it in a list of lists where each entry is [TORS, FREQ]
def torenv_from_file(input_filename):
    full_list=[]
    with open(input_filename,"r") as f:
        for line in f:
            linesplit_list=line.rsplit("-",1)
            linesplit_list[1]=linesplit_list[1].strip()
            full_list.append([linesplit_list[0],int(linesplit_list[1])])
    return(full_list)
 


# Takes in the "parent" list and gives each entry a unique ID based on the
# highest frequency. This is to compare against later to keep the same
# ordering as this parent_list
def prep_parent_list(parent_list):
    #sorts the lists by frequency, largest frequency first
    parent_list = sorted(parent_list, key=itemgetter(1), reverse=True)

    #Gives each entry in dock_list an identifier, with largest frequency=1
    entry_number=1
    for entry in parent_list:
        entry.append(entry_number)
        entry_number+=1
    return(parent_list)

# Normalizes a torenv list against itself and stores entries in X and Y lists
# for graphing. MUST BE DONE AFTER IDing
def normalize_torenv_list(input_list):
    x=[]
    y=[]

    for entry in range(0, len(input_list)):
        x.append(input_list[entry][2])
        y.append(input_list[entry][1])
    y_summer=0
    for entry in y:
        y_summer+=entry
    for val in range(0, len(y)):
        y[val]=y[val]/y_summer

    return((x,y))

# Compares a "child" environment against the parent, filling in holes where
# there are no listed torsions
def compare_against_parent(parent_list, compare_list):
    for i in parent_list:
        found = False
        for y in compare_list:
            if i[0] == y[0]:
                y.append(i[2])
                found = True
                break
        if found == False:
            compare_list.append([i[0],0,i[2]])
    return(compare_list)


def plot_envs_w_parent(x_list, plot_data, plot_names, graph_title,
                            colors, file_out_name, axes):
    if len(plot_data) != len(plot_names):
        print("The number of data sets does not equal the number of names.")
        exit(1)

    plt.style.use('default')
    fig = plt.figure(num=None, \
                    figsize=(8,2), \
                    dpi=400, \
                    facecolor='w', \
                    edgecolor='k')
    plt.axis(axes)
    
    plt.xlabel('Torsion ID', size=5)
    plt.ylabel('Rel. Frequency', size=5)
    plt.title(graph_title, size=8)
    plt.tick_params(axis='both', which='major', labelsize=6)

    for i in range(1,len(plot_data)):
        plt.plot(x_list,plot_data[i],color=colors[i],label=plot_names[i],\
            linewidth=0.6)
        
    plt.plot(x_list, plot_data[0],color=colors[0],label=plot_names[0],\
            linewidth=0.6)
    plt.legend(loc="upper right")
    plt.savefig(file_out_name,dpi=300,transparent=False)
    plt.close()

def plot_top_and_bottom(x_list, plot_data, plot_names, graph_title,
                    colors, file_out_name, axes):
    fig, (ax1, ax2) = plt.subplots(2,gridspec_kw={'hspace': 0},\
                    dpi=150,\
                    figsize=(8,2))

    ax1.axis(axes)
    ax2.axis(axes)

    for ax in fig.get_axes():
        ax.label_outer()

    fig.suptitle(graph_title)
    fig.axes[1].set_xlabel('Torsion ID',size=5)
    ax1.tick_params(axis='both',which='major',labelsize=5)
    ax2.tick_params(axis='both',which='major',labelsize=5)


    line_gr0=ax1.plot(x_list, plot_data[1], color=colors[1],\
        label=plot_names[1], linewidth=0.6)
    line_other=ax2.plot(x_list, plot_data[2], color=colors[2],\
        label=plot_names[2], linewidth=0.6)
    line_orig=ax1.plot(x_list, plot_data[0], color=colors[0],\
        label=plot_names[0], linewidth=0.4)
    ax2.plot(x_list, plot_data[0],color='black',\
        label="Input Library",linewidth=0.4)

    ax1.legend(handles=[line_orig[0],line_gr0[0],line_other[0]],\
        loc="upper right", fontsize="x-small")
    plt.savefig(file_out_name,dpi=300,transparent=False)
    plt.close()








    
