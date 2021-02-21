import os 
import sys
import pandas as pd 
import collections
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from  Bio import SeqIO
from matplotlib.transforms import Bbox

matplotlib.rcParams['font.sans-serif']   = ["Helvetica","Arial","Lucida Sans","DejaVu Sans","Lucida Grande","Verdana"]
matplotlib.rcParams['font.family']       = 'sans-serif'
matplotlib.rcParams['font.sans-serif']   = ["Helvetica","Arial","DejaVu Sans","Lucida Grande","Verdana"]
matplotlib.rcParams['font.size']         = 12.0
matplotlib.rcParams["axes.labelcolor"]   = "#000000"
matplotlib.rcParams["axes.linewidth"]    = 1.0
matplotlib.rcParams["xtick.major.width"] = 1.0
matplotlib.rcParams["ytick.major.width"] = 1.0
matplotlib.rcParams['xtick.major.pad']   = 6
matplotlib.rcParams['ytick.major.pad']   = 6
matplotlib.rcParams['xtick.major.size']  = 6
matplotlib.rcParams['ytick.major.size']  = 6

def map_feat(fig, ax, feats, length, head_length=np.pi * 0.02, unvisible_types=["source"], visible_types=[], enlarge=1.0):
    y_list    = [] 
    ty_list   = [] 
    visible   = 1
    unvisible = 1 
    gene_position_matrix = [[0] * length]
    text_position_matrix = [[0] * length] 
    label_position_list = [] 
    if len(unvisible_types) == 0:
        visible   = 1
        unvisible = 1
    else:
        unvisible = 0
    
    if len(visible_types) == 0:
        visible = 1
    else:
        visible   = 0 
        unvisible = 0
        unvisible_types = []

    for i, feat in enumerate(feats):
        if feat.type in unvisible_types:
            pass
        
        elif visible == 1 or feat.type in visible_types:
            if "label" in feat.qualifiers:
                if type(feat.qualifiers["label"]) == list:
                    label = feat.qualifiers["label"][0]
                else:
                    label = feat.qualifiers["label"]
            else:
                label = feat.type
            strand = feat.location.strand
            if strand == 1:
                gs_origin = feat.location.parts[0].start.position 
                ge_origin = feat.location.parts[-1].end.position
                gs = gs_origin * 2 * np.pi / length 
                ge = ge_origin * 2 * np.pi / length
            else:
                gs_origin = feat.location.parts[-1].start.position 
                ge_origin = feat.location.parts[0].end.position
                gs = gs_origin * 2 * np.pi / length
                ge = ge_origin * 2 * np.pi / length 
            
            if i > 0:
                flag = 0
                if gs_origin < ge_origin:
                    for y, row in enumerate(gene_position_matrix):
                        if 1 in row[gs_origin:ge_origin]:
                            flag = 1
                        else:
                            flag = 0
                            break    
                    if flag == 1:
                        y += 1 
                        gene_position_matrix.append([0] * length)
                    else:
                        pass 
                else:
                    for y, row in enumerate(gene_position_matrix):
                        if 1 in row[gs_origin:]:
                            flag = 1
                        else:
                            flag = 0
                            break    

                    for y, row in enumerate(gene_position_matrix):
                        if 1 in row[:ge_origin]:
                            flag = 1
                        else:
                            flag = 0
                            break    

                    if flag == 1:
                        y += 1 
                        gene_position_matrix.append([0] * length)
                    else:
                        pass 

            if abs(ge-gs) < head_length * 1.2:
                hl  = abs(ge-gs)
                margin1 = head_length * 0.15 /enlarge
                margin2 = head_length * 0.25 /enlarge
                hl2 = hl - (margin1+margin2)
                if hl2 < 0:
                    hl2 = 0
                    wd2 = 0 
                else:
                    wd2 = 0.8 * ((hl2/hl) ** 1.0)
            else:
                hl  = head_length 
                hl2 = hl*0.85  
                margin1 = hl * 0.15 /enlarge
                margin2 = hl * 0.25 /enlarge
                wd2     = 0.65

            if "facecolor_dbrick" in feat.qualifiers:
                if type(feat.qualifiers["facecolor_dbrick"]) == list:
                    facecolor = feat.qualifiers["facecolor_dbrick"][0] 
                else:
                    facecolor = feat.qualifiers["facecolor_dbrick"]
            else:
                facecolor = "#ffffec" 
            
            if "edgecolor_dbrick" in feat.qualifiers:
                if type(feat.qualifiers["edgecolor_dbrick"]) == list:
                    edgecolor = feat.qualifiers["edgecolor_dbrick"][0] 
                else:
                    edgecolor = feat.qualifiers["edgecolor_dbrick"]
            else:
                if strand == 1:
                    edgecolor = "#FACAC8" 
                elif strand == -1:
                    edgecolor = "#C8DFFA" 
                else:
                    edgecolor = "#CCCCCC"
                
            if "note_dbrick" in feat.qualifiers:
                note   = feat.qualifiers["note_dbrick"]
                pos_s  = int(note.split(":")[0].split("..")[0])
                pos_e  = int(note.split(":")[0].split("..")[1]) 
                feat_length = int(note.split(":")[1])
                if (pos_s != 1 or pos_e != feat_length):
                    label = label + ":" + str(pos_s) + ".." + str(pos_e) + ":" + str(feat_length)
            
            if gs < ge: 
                width = ge-gs
            else:
                width = ge + (2*np.pi-gs) 
            
            margin = np.pi*0.004
            if width > head_length:
                if strand == 1:
                    ax.bar([gs], [60], bottom=y*75+370, width=width-head_length*0.98, align="edge", fc=edgecolor, lw=0.0)
                    ax.bar([gs+margin], [45], bottom=y*75+377.5, width=width-margin-head_length*0.98, align="edge", fc=facecolor, lw=0.0)
                    ax.arrow(x=ge-head_length, y=y*75+400, dx=head_length, dy=0, width=60, head_width=60, head_length=head_length, length_includes_head=True, fc=edgecolor, lw=0.0)
                    ax.arrow(x=ge-head_length, y=y*75+400, dx=head_length-margin*1.5, dy=0, width=45, head_width=45, head_length=head_length-margin*1.5, length_includes_head=True, fc=facecolor, lw=0.0)
                    #ax.arrow(x=gs+margin1, y=y*50+500, dx=ge-gs-(margin1+margin2), dy=0, width=wd2, head_width=wd2, head_length=hl2, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                elif strand == -1:
                    ax.bar([gs+head_length*0.95], [60], bottom=y*75+370, width=width-head_length*0.95, align="edge", fc=facecolor, ec=edgecolor, lw=1.0)
                    ax.arrow(x=gs+head_length, y=y*75+400, dx=-1*head_length, dy=0, width=60, head_width=60, head_length=head_length, length_includes_head=True, fc=facecolor, ec=edgecolor)
                    #ax.arrow(x=ge-margin1, y=y*50+500, dx=gs-ge+(margin1+margin2), dy=0, width=wd2, head_width=wd2, head_length=hl2, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                else:
                    ax.bar([gs], [60], bottom=y*75+380, width=width, align="edge", fc=facewcolor, edgecolor=edgecolor, lw=1.0)
                    #ax.arrow(x=gs, y=y*50+500, dx=ge-gs, dy=0, width=50, head_width=50, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                    #ax.arrow(x=gs+margin1, y=y*50+500, dx=ge-gs-(margin1+margin2), dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
            else:
                if strand == 1:
                    ax.arrow(x=ge-width, y=y*75+400, dx=width, dy=0, width=60, head_width=60, head_length=width, length_includes_head=True, color='k', fc=facecolor, ec=edgecolor, lw=1.0)
                    #ax.arrow(x=gs+margin1, y=y*50+500, dx=ge-gs-(margin1+margin2), dy=0, width=wd2, head_width=wd2, head_length=hl2, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                
                elif strand == -1:
                    ax.arrow(x=gs+width, y=y*75+400, dx=-1*width, dy=0, width=60, head_width=60, head_length=width, length_includes_head=True, color='k', fc=facecolor, ec=edgecolor, lw=1.0)
                    #ax.arrow(x=ge-margin1, y=y*50+500, dx=gs-ge+(margin1+margin2), dy=0, width=wd2, head_width=wd2, head_length=hl2, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                
                else:
                    ax.bar([gs], [60], bottom=y*75+400, width=width, align="edge", fc=facecolor, ec=edgecolor, lw=0.0)
                    #ax.arrow(x=gs+margin1, y=y*50+500, dx=ge-gs-(margin1+margin2), dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)

            label_position_list.append((label,(gs+ge)/2, y, gs, ge, hl, edgecolor, facecolor))
            if gs_origin < ge_origin:
                for j in range(gs_origin,ge_origin):
                    gene_position_matrix[y][j] = 1  
            else:
                for j in range(gs_origin,length):
                    gene_position_matrix[y][j] = 1 
                for j in range(0, ge_origin):
                    gene_position_matrix[y][j] = 1

            y_list.append(y)

    return ax, y_list, ty_list

def visualize(brick):
    figure = plt.figure(figsize=(5,5))
    ax     = plt.subplot(111, polar=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.spines['polar'].set_visible(False)
    ax.xaxis.set_ticks([])
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticks([])
    ax.yaxis.set_ticklabels([])
    ax, y_list, ty_list = map_feat(figure, ax, brick.features, len(brick.seq), unvisible_types=["source"], visible_types=[], enlarge=1.0)
    if max(y_list)*75+500 < 600:
        ax.set_ylim(0,600)
    else:
        ax.set_ylim(0,max(y_list)*50+450)
    
    return figure 

if __name__ == "__main__":
    from dbrick import *
    brick = Dbrick(record="./notebook/tutorial/pUC19.gbk") 
    fig   = visualize(brick) 
    fig.patch.set_alpha(0.0) 
    fig.savefig("test_circular.pdf", bbox_inches="tight")
