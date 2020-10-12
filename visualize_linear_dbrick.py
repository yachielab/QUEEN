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

def map_feat(fig, ax, feats, length, head_length, unvisible_types=["source"], visible_types=[], enlarge=1.0):
    y = 0
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
                gs = feat.location.parts[0].start.position 
                ge = feat.location.parts[-1].end.position
            else:
                gs = feat.location.parts[-1].start.position 
                ge = feat.location.parts[0].end.position
            
            if i > 0:
                flag = 0
                for y, row in enumerate(gene_position_matrix):
                    if 1 in row[gs:ge]:
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
                if type(note) == list:
                    note   = feat.qualifiers["note_dbrick"][0]
                if (pos_s != 1 or pos_e != feat_length):
                    label = note

            if strand == 1:
                if "note_dbrick" in feat.qualifiers:
                    if (gs == 0 or ge >= length) and (pos_s != 1 or pos_e != feat_length):
                        if gs == 0 and ge >= length:
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=0.8, head_width=0.8, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                        elif gs == 0:
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=0.8, head_width=0.8, head_length=hl, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs-margin2, dy=0, width=wd2, head_width=wd2, head_length=hl2, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                        else:
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=0.8, head_width=0.8, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=gs+margin1, y=-1*y, dx=ge-gs-margin1, dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                    else:
                        ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=0.8, head_width=0.8, head_length=hl, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                        ax.arrow(x=gs+margin1, y=-1*y, dx=ge-gs-(margin1+margin2), dy=0, width=wd2, head_width=wd2, head_length=hl2, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                else:
                    ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=0.8, head_width=0.8, head_length=hl, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                    ax.arrow(x=gs+margin1, y=-1*y, dx=ge-gs-(margin1+margin2), dy=0, width=wd2, head_width=wd2, head_length=hl2, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
            
            elif strand == -1:
                if "note_dbrick" in feat.qualifiers:
                    if (gs == 0 or ge >= length) and (pos_s != 1 or pos_e != feat_length):
                        if gs == 0 and ge >= length:
                            ax.arrow(x=ge, y=-1*y, dx=gs-ge, dy=0, width=0.8, head_width=0.8, head_length=hl, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=ge, y=-1*y, dx=gs-ge, dy=0, width=wd2, head_width=wd2, head_length=hl2, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                        elif gs == 0:
                            ax.arrow(x=ge, y=-1*y, dx=gs-ge, dy=0, width=0.8, head_width=0.8, head_length=hl, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=ge-margin1, y=-1*y, dx=gs-ge+margin1, dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                        else:
                            ax.arrow(x=ge, y=-1*y, dx=gs-ge, dy=0, width=0.8, head_width=0.8, head_length=hl, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=ge, y=-1*y, dx=gs-ge+margin2, dy=0, width=wd2, head_width=wd2, head_length=hl2, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                    else:
                        ax.arrow(x=ge, y=-1*y, dx=gs-ge, dy=0, width=0.8, head_width=0.8, head_length=hl, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                        ax.arrow(x=ge-margin1, y=-1*y, dx=gs-ge+(margin1+margin2), dy=0, width=wd2, head_width=wd2, head_length=hl2, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                else:
                    ax.arrow(x=ge, y=-1*y, dx=gs-ge, dy=0, width=0.8, head_width=0.8, head_length=hl, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                    ax.arrow(x=ge-margin1, y=-1*y, dx=gs-ge+(margin1+margin2), dy=0, width=wd2, head_width=wd2, head_length=hl2, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
            
            else:
                if "note_dbrick" in feat.qualifiers:
                    if (gs == 0 or ge >= length) and (pos_s != 1 or pos_e != feat_length):
                        if gs == 0 and ge >= length:
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=0.8, head_width=0.8, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                        elif gs == 0:
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=0.8, head_width=0.8, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs-margin2, dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                        else:
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=0.8, head_width=0.8, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=gs+margin1, y=-1*y, dx=ge-gs-margin1, dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                    else:
                        ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=0.8, head_width=0.8, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                        ax.arrow(x=gs+margin1, y=-1*y, dx=ge-gs-(margin1+margin2), dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                else: 
                    ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=0.8, head_width=0.8, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                    ax.arrow(x=gs+margin1, y=-1*y, dx=ge-gs-(margin1+margin2), dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
            
            label_position_list.append((label,(gs+ge)/2, -1*y, gs, ge, hl, edgecolor, facecolor))
            for j in range(gs,ge):
                gene_position_matrix[y][j] = 1    
            y_list.append(y)
    
    renderer     = fig.canvas.get_renderer()
    coordinate  = ax.transData.inverted() 
    for label, tx, ty, gs, ge, hl, ec, fc in label_position_list:
        text      = ax.text(tx, ty, label, ha="center", va="center")
        bbox_text = text.get_window_extent(renderer=renderer)
        bbox_text = Bbox(coordinate.transform(bbox_text))
        width = bbox_text.width*length
        ts    = int(tx-int(width/2)*1.2) 
        ts    = 0 if ts < 0 else ts 
        te    = int(tx+int(width/2)*1.2)
        te    = length if te > length else te
        if width > abs(ge-gs) * 0.9:
            text.set_visible(False)
            for t, trow in enumerate(text_position_matrix):
                if 1 not in trow[ts:te]:
                    flag = 0 
                    break
                else:
                    flag = 1
            
            if flag  == 1:
                text_position_matrix.append([0] * length)
                t += 1       
            
            if t % 2 == 0:
                ax.text(tx, t//2+1, label, ha="center",va="center", bbox=dict(boxstyle="round", ec=ec, fc=fc), zorder=100)
                ax.plot([(gs+ge)/2, (gs+ge)/2], [t//2+1, ty], lw=0.5, color="k", zorder=0)
                ty_list.append(t//2+1)
            
            else:
                ax.text(tx, -1*max(y_list)-2.0-t//2, label, ha="center",va="center",  bbox=dict(boxstyle="round", ec=ec, fc=fc), zorder=100)
                ax.plot([(gs+ge)/2, (gs+ge)/2], [-1*max(y_list)-2.0-t//2, ty], lw=0.5, color="k", zorder=0)
                ty_list.append(-1*max(y_list)-2.0-t//2)
            
            for j in range(ts, te):
                text_position_matrix[t][j] = 1
        else:
            ty_list.append(ty) 
    return ax, y_list, ty_list
    
def visualize(brick, start=0, end=None, aspect=None, unvisible_types=["source"], visible_types=[], enlarge=1.0):
    if end == None:
        end = len(brick.seq) 
    sub_brick  = brick[start:end]

    if start == 0 and end == len(brick.seq):
        pass 
    else:
        start, end = sub_brick.start, sub_brick.end 
    
    fig = plt.figure(figsize=(3,0.35))
    matplotlib.rcParams['font.size'] = 12.0
    if len(sub_brick.seq) < 251:
        std = 50 
        head_length = 3
    elif len(sub_brick.seq) < 501:
        std = 100
        head_length = 6
    elif len(sub_brick.seq) < 1001:
        std = 200 
        head_length = 12
    elif len(sub_brick.seq) < 2001:
        std = 400
        head_length = 24
    elif len(sub_brick.seq) < 4001:
        std = 800
        head_length = 48
    else:
        std = 1600
        head_length = 48
    
    zero_position = start + 1 
    ax  = fig.add_axes([0, 0, enlarge*len(sub_brick.seq)/std, 1.0]) 
    ax, y_list, ty_list = map_feat(fig, ax, sub_brick.features, len(sub_brick.seq), head_length, unvisible_types=unvisible_types, visible_types=visible_types, enlarge=enlarge) 
    ty_list.append(0) 
    
    ax.spines["top"].set_visible(False) 
    ax.spines["right"].set_visible(False) 
    ax.spines["left"].set_visible(False) 
    ax.spines["bottom"].set_position(("data",-1.0 * max(y_list)-0.5))
    positions   = [p+zero_position-len(brick.seq) if p+zero_position > len(brick.seq) else p+zero_position for p in range(0, len(sub_brick.seq))] 
    ticks       = [x for p, x in zip(positions, list(range(0, len(sub_brick.seq)))) if p % int(std/(2*enlarge)) == 0 or p == 1] 
    tick_labels = [str(p) for p, x in zip(positions, list(range(0, len(sub_brick.seq)))) if p % int(std/(2*enlarge)) == 0 or  p == 1]
    ax.set_xticks(ticks) 
    ax.set_xticklabels(tick_labels)
    ax.set_xlim(0, len(sub_brick.seq))
    ax.set_ylim(min(ty_list)-0.5, max(ty_list)+0.5) 
    ax.set_yticks([])
    fig.set_size_inches(3, 0.35*(abs(max(ty_list)-min(ty_list))+1.0))
    return fig 
