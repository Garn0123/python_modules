import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt # plotting package
import os
import sys
from tqdm import tqdm

#extracts mol_props data and converts all data to floats
#returns list of lists for all experiments provided
def extract_molecular_graphing_data(experiment_list,workdir):
    data_collected=[]
    for experiment in tqdm(experiment_list):
        data_list=[]
        #data extraction from previously saved data
        with open(f"{workdir}/{experiment}/{experiment}_props.dat") as f:
            for line in f:
                data_list.append(line.strip().split(","))
        #converts all values to floats
        for i,entry in enumerate(data_list):
            data_list[i] = [float(x) for x in entry]
        data_collected.append(data_list)
    return(data_collected)

def graph_9panel_mol_props_by_sasa(experiment_list,workdir,data_collected,deliniation):
    graph_values=["MW","LOGP","HBA","HBD","ROT","ARO","PAINS","QED","SASA","SPIR","STER","PSA"]
    graph_titles=["Molecular Weight",\
                "ALogP",\
                "HBond Acceptors",\
                "HBond Donors",\
                "# Rotatable Bonds",\
                "# Aromatic Rings",\
                "# Structural Alerts",\
                "QED",\
                "Novartis Synthetic Accessibility",\
                "Spiro Atoms",\
                "Stereo Centers",\
                "PSA"]

    #this function graphs by SASA, which here is 0 to 10
    for i in tqdm(range(0,10)):
        fig,axs = plt.subplots(4,3,figsize=(20,15),dpi=500,sharey=True)
        pos_x=0
        pos_y=0
        graph_val_position=0
        pos = np.arange(len(experiment_list)+1,1,-1)
        for val_ind, value in enumerate(graph_values):
            data_list=[]
            #separates out all molecule data by the SynthA
            for exp_ind, experiment in enumerate(experiment_list):
                #data_list.append([x[val_ind] for x in data_collected[exp_ind]])
                data_list.append([x[val_ind] for x in data_collected[exp_ind] \
                                if x[8] > i and x[8] < i+1])
            for index, entry in enumerate(data_list):
                #Just in the case there are NO MOLECULES MADE for this value, 
                if len(entry) == 0:
                    data_list[index]=[float('nan'), float('nan')]
            if value == "MW":
                plt.setp(axs[pos_x,pos_y], xlim=(0,850))
                axs[pos_x,pos_y].violinplot(data_list,pos,vert=False,showmeans=True)
            if value == "LOGP":
                plt.setp(axs[pos_x,pos_y], xlim=(-8,10))
                axs[pos_x,pos_y].violinplot(data_list,pos,vert=False,showmeans=True)
            if value == "HBA" or value == "ROT":
                plt.setp(axs[pos_x,pos_y], xlim=(0,16))
                axs[pos_x,pos_y].boxplot(data_list,positions=pos,vert=False)
            if (value == "HBD") or (value == "ARO"):
                plt.setp(axs[pos_x,pos_y], xlim=(0,10))
                axs[pos_x,pos_y].boxplot(data_list,positions=pos,vert=False)
            if value == "SASA":
                plt.setp(axs[pos_x,pos_y], xlim=(0,10))
                axs[pos_x,pos_y].violinplot(data_list,pos,vert=False,showmeans=True)
            if (value == "STER") or (value =="SPIR"):
                plt.setp(axs[pos_x,pos_y], xlim=(0,5))
                axs[pos_x,pos_y].boxplot(data_list,positions=pos,vert=False)
            if (value == "PAINS"):
                plt.setp(axs[pos_x,pos_y], xlim=(0,8))
                axs[pos_x,pos_y].boxplot(data_list,positions=pos,vert=False)
            if value == "QED":
                plt.setp(axs[pos_x,pos_y], xlim=(0,1))
                axs[pos_x,pos_y].violinplot(data_list,pos,vert=False,showmeans=True)
            if value == "PSA":
                plt.setp(axs[pos_x,pos_y], xlim=(0,300))
                axs[pos_x,pos_y].violinplot(data_list,pos,vert=False,showmeans=True)
            axs[pos_x,pos_y].set_title(graph_titles[graph_val_position])
            axs[pos_x,pos_y].set_yticks(pos)
            axs[pos_x,pos_y].set_yticklabels(experiment_list)
            if (pos_y > 0) and (pos_y % 2 == 0):
                pos_x+=1
                pos_y=0
            else:
                pos_y+=1
            graph_val_position+=1
            fig.subplots_adjust(top=0.94,hspace = 0.2)
        st = fig.suptitle(f"{experiment_list[0]}_{deliniation}_{i}_to_{i+1}",fontsize=20)
        st.set_y(1.00)
        plt.savefig(f"{workdir}/{experiment_list[0]}/zzz.{experiment_list[0]}_{deliniation}_sasa_{i}.png")
        plt.close()




def graph_9panel_mol_props(experiment_list,workdir,data_collected,deliniation):
    graph_values=["MW","LOGP","HBA","HBD","ROT","ARO","PAINS","QED","SASA","SPIR","STER","PSA"]
    graph_titles=["Molecular Weight",\
                "ALogP",\
                "HBond Acceptors",\
                "HBond Donors",\
                "# Rotatable Bonds",\
                "# Aromatic Rings",\
                "# Structural Alerts",\
                "QED",\
                "Novartis Synthetic Accessibility",\
                "Spiro Atoms",\
                "Stereo Centers",\
                "PSA"]
    fig,axs = plt.subplots(4,3,figsize=(20,15),dpi=500,sharey=True)
    pos_x=0
    pos_y=0
    pos = np.arange(len(experiment_list)+1,1,-1)
    for val_ind, value in enumerate(graph_values):
        data_list=[]
        #separates out all molecule data by the SynthA
        for exp_ind, experiment in enumerate(experiment_list):
            data_list.append([x[val_ind] for x in data_collected[exp_ind]])
        for index, entry in enumerate(data_list):
            #Just in the case there are NO MOLECULES MADE for this value, 
            if len(entry) == 0:
                data_list[index]=[float('nan'), float('nan')]
        if value == "MW":
            plt.setp(axs[pos_x,pos_y], xlim=(0,850))
            axs[pos_x,pos_y].violinplot(data_list,pos,vert=False,showmeans=True)
        if value == "LOGP":
            plt.setp(axs[pos_x,pos_y], xlim=(-8,10))
            axs[pos_x,pos_y].violinplot(data_list,pos,vert=False,showmeans=True)
        if value == "HBA" or value == "ROT":
            plt.setp(axs[pos_x,pos_y], xlim=(0,16))
            axs[pos_x,pos_y].boxplot(data_list,positions=pos,vert=False)
        if (value == "HBD") or (value == "ARO"):
            plt.setp(axs[pos_x,pos_y], xlim=(0,10))
            axs[pos_x,pos_y].boxplot(data_list,positions=pos,vert=False)
        if value == "SASA":
            plt.setp(axs[pos_x,pos_y], xlim=(0,10))
            axs[pos_x,pos_y].violinplot(data_list,pos,vert=False,showmeans=True)
        if (value == "STER") or (value =="SPIR"):
            plt.setp(axs[pos_x,pos_y], xlim=(0,5))
            axs[pos_x,pos_y].boxplot(data_list,positions=pos,vert=False)
        if (value == "PAINS"):
                plt.setp(axs[pos_x,pos_y], xlim=(0,8))
                axs[pos_x,pos_y].boxplot(data_list,positions=pos,vert=False)
        if value == "QED":
            plt.setp(axs[pos_x,pos_y], xlim=(0,1))
            axs[pos_x,pos_y].violinplot(data_list,pos,vert=False,showmeans=True)
        if value == "PSA":
            plt.setp(axs[pos_x,pos_y], xlim=(0,300))
            axs[pos_x,pos_y].violinplot(data_list,pos,vert=False,showmeans=True)
        axs[pos_x,pos_y].set_title(graph_values[val_ind])
        axs[pos_x,pos_y].set_yticks(pos)
        axs[pos_x,pos_y].set_yticklabels(experiment_list)
        if (pos_y > 0) and (pos_y % 2 == 0):
            pos_x+=1
            pos_y=0
        else:
            pos_y+=1
        #print(pos_x)
        #print(pos_y)
        fig.subplots_adjust(top=0.94,hspace = 0.2)
    st = fig.suptitle(f"{experiment_list[0]}_{deliniation}",fontsize=20)
    st.set_y(1.00)
    plt.savefig(f"{workdir}/{experiment_list[0]}/zzz.{experiment_list[0]}_{deliniation}.png")
    plt.close()

def correlation_plot(dataset1,dataset2,name1,name2,axes,title,fileout):
    correlation_matrix = np.corrcoef(dataset1,dataset2)
    correlation_xy = correlation_matrix[0,1]
    r_squared = correlation_xy**2

    fix, ax = plt.subplots()
    ax.scatter(dataset1, dataset2,s=1)
    ax.set_xlabel(name1)
    ax.set_ylabel(name2)
    ax.axis(axes)
    plt.title(title)
    plt.text(0,0.95, 
        str("RSquared = %.4f" % r_squared), 
        fontsize = 12, 
        transform = ax.transAxes,
        horizontalalignment='left')
    plt.savefig(f"{fileout}_correlation.png", dpi=300)
    plt.close()
    return

def histogram_plot(data,x_label,title,fileout,bins):
    fix, ax = plt.subplots()
    ax.hist(data,density=True)
    ax.set_xlabel(x_label)
    plt.title(title)
    ax.set_ylim([0,1])
    plt.savefig(f"{fileout}_histogram.png",dpi=300)
    plt.close()
    return
