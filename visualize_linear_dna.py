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
matplotlib.rcParams['xtick.major.pad']   = 4
matplotlib.rcParams['ytick.major.pad']   = 4
matplotlib.rcParams['xtick.major.size']  = 4
matplotlib.rcParams['ytick.major.size']  = 4
color_dict = {"G":"#f2f059", "C":"#74b2d7", "A":"#79E5B7", "T":"#ff776c", "N":"#FFFFFF", "-":"#FFFFFF"}

feature_color_dict = collections.defaultdict(list)
feature_color_dict["CDS"]          = [('#92c6ff', '#4c72b0'), ('#97f0aa', '#55a868'), ('#ff9f9a', '#c44e52'), ('#d0bbff', '#8172b2'), ('#fffea3', '#ccb974')]  
feature_color_dict["primer_bind"]  = [('#a6cee3', '#1f78b4'), ('#b2df8a', '#33a02c'), ('#fb9a99', '#e31a1c'), ('#fdbf6f', '#ff7f00'), ('#cab2d6', '#6a3d9a')] 
feature_color_dict["primer"]       = [('#a6cee3', '#1f78b4'), ('#b2df8a', '#33a02c'), ('#fb9a99', '#e31a1c'), ('#fdbf6f', '#ff7f00'), ('#cab2d6', '#6a3d9a')] 
feature_color_dict["promoter"]     = [('#b0e0e6', '#64b5cd'), ('#92c6ff', '#4c72b0'), ('#97f0aa', '#55a868')]  
feature_color_dict["rep_origin"]   = [('#fff2ae', '#ffd92f')]
feature_color_dict["misc_feature"] = [('#fbb4ae', '#e41a1c'), ('#b3cde3', '#377eb8'), ('#ccebc5', '#4daf4a'), ('#decbe4', '#984ea3'), ('#fed9a6', '#ff7f00')]  
misc_colors = [('#ffffcc', '#d9d927'), ('#e5d8bd', '#a65628'), ('#fddaec', '#f781bf'), ('#f2f2f2', '#999999'), ('#fbb4ae', '#e41a1c')]

def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i+lv//3], 16) for i in range(0, lv, lv//3))

def colorbar(ax,color_dict,ref_seq,char=False):
    bars = ax.bar(list(range(len(ref_seq))), [0.9] * (len(ref_seq)), width=1.0, edgecolor="#A0A0A0", linewidth=0.2, align="edge",bottom=0.05)
    ax.set_xlim(0,len(ref_seq))
    ax.set_ylim(0,1.00)
    p = 0
    for bar, c in zip(bars,ref_seq):
        color = color_dict[c]
        if char == True:
            color = hex_to_rgb(color)
            color = color + (0.75,)
            bar.set_facecolor(color)
            ax.text(p+0.5,0.38,c,va="center",ha="center",fontsize=9,zorder=100)     
        else:
            bar.set_facecolor(color)
        p += 1 
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.patch.set_alpha(0.0)
    return bars

def map_feat(fig, ax, feats, length, head_length, unvisible_types=["source"], visible_types=[], enlarge_w=1.0, enlarge_h=1.0, annotation_loc="both"):
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
                margin1 = head_length * 0.15 /enlarge_w
                margin2 = head_length * 0.25 /enlarge_w 
                hl2 = hl - (margin1+margin2)
                if hl2 < 0:
                    hl2 = 0
                    wd2 = 0 
                else:
                    wd2 = 0.8 * ((hl2/hl) ** 1.0)
            else:
                hl  = head_length 
                hl2 = hl*0.85  
                margin1 = hl * 0.15 /enlarge_w
                margin2 = hl * 0.25 /enlarge_w
                wd2     = 0.65

            if "facecolor_dna.py" in feat.qualifiers:
                if type(feat.qualifiers["facecolor_dna.py"]) == list:
                    facecolor = feat.qualifiers["facecolor_dna.py"][0] 
                else:
                    facecolor = feat.qualifiers["facecolor_dna.py"]
            else:
                facecolor = "#ffffec" 
            
            if "edgecolor_dna.py" in feat.qualifiers:
                if type(feat.qualifiers["edgecolor_dna.py"]) == list:
                    edgecolor = feat.qualifiers["edgecolor_dna.py"][0] 
                else:
                    edgecolor = feat.qualifiers["edgecolor_dna.py"]
            else:
                if strand == 1:
                    edgecolor = "#FACAC8" 
                elif strand == -1:
                    edgecolor = "#C8DFFA" 
                else:
                    edgecolor = "#CCCCCC"
                
            if "note_dna.py" in feat.qualifiers:
                note   = feat.qualifiers["note_dna.py"]
                if type(note) == list:
                    note   = feat.qualifiers["note_dna.py"][0]
                
                pos_s  = int(note.split(":")[1].split("..")[0]) 
                pos_e  = int(note.split(":")[1].split("..")[1])
                feat_length = int(note.split(":")[2])
                if (pos_s != 1 or pos_e != feat_length):
                    label = note

            if strand == 1:
                if "note_dna.py" in feat.qualifiers:
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
                if "note_dna.py" in feat.qualifiers:
                    if (gs == 0 or ge >= length) and (pos_s != 1 or pos_e != feat_length):
                        if gs == 0 and ge >= length:
                            ax.arrow(x=ge, y=-1*y, dx=gs-ge, dy=0, width=0.8, head_width=0.8, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=ge, y=-1*y, dx=gs-ge, dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                        elif gs == 0:
                            ax.arrow(x=ge, y=-1*y, dx=gs-ge, dy=0, width=0.8, head_width=0.8, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
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
                if "note_dna.py" in feat.qualifiers:
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
        fc = "#ffffec" #facecolor 
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
            
            if annotation_loc == "both": 
                if t % 2 == 0:
                    ax.text(tx, t//2+1, label, ha="center",va="center", bbox=dict(boxstyle="round", ec=ec, fc=fc), zorder=100)
                    ax.plot([(gs+ge)/2, (gs+ge)/2], [t//2+1, ty], lw=0.5, color="k", zorder=0)
                    ty_list.append(t//2+1)
                
                else:
                    ax.text(tx, -1*max(y_list)-2.0-t//2, label, ha="center",va="center",  bbox=dict(boxstyle="round", ec=ec, fc=fc), zorder=100)
                    ax.plot([(gs+ge)/2, (gs+ge)/2], [-1*max(y_list)-2.0-t//2, ty], lw=0.5, color="k", zorder=0)
                    ty_list.append(-1*max(y_list)-2.0-t//2)
            
            elif annotation_loc == "bottom":
                ax.text(tx, -1*max(y_list)-2.0-t, label, ha="center",va="center",  bbox=dict(boxstyle="round", ec=ec, fc=fc), zorder=100)
                ax.plot([(gs+ge)/2, (gs+ge)/2], [-1*max(y_list)-2.0-t, ty], lw=0.5, color="k", zorder=0)
                ty_list.append(-1*max(y_list)-2.0-t)
                
            else:
                ax.text(tx, t+1, label, ha="center",va="center", bbox=dict(boxstyle="round", ec=ec, fc=fc), zorder=100)
                ax.plot([(gs+ge)/2, (gs+ge)/2], [t+1, ty], lw=0.5, color="k", zorder=0)
                ty_list.append(t+1)

            for j in range(ts, te):
                text_position_matrix[t][j] = 1
        else:
            ty_list.append(ty) 
    return ax, y_list, ty_list

def colorbar(ax, color_dict, ref_seq, char=False):
    bars = ax.bar(list(range(len(ref_seq))), [0.9] * (len(ref_seq)), width=1.0, edgecolor="#BBBBBB", linewidth=0.2, align="edge",bottom=0.05)
    ax.set_xlim(0,len(ref_seq))
    ax.set_ylim(0,1.00)
    p = 0
    for bar, c in zip(bars,ref_seq):
        color = color_dict[c]
        bar.set_facecolor(color)
        if char == True:
            bar.set_alpha(0.7)
            ax.text(p+0.5,0.45,c,va="center",ha="center",fontsize=10,zorder=100)     
        p += 1 
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.patch.set_alpha(0.0)
    return bars

def visualize(brick, start=0, end=None, wrap_width=None, annotation_loc=None, unvisible_types=["source"], visible_types=[], enlarge_w=1.0, enlarge_h=1.0, fontsize=12, with_seq=False, nucl_char=None, nucl_color_dict=None):
    width = wrap_width 
    if nucl_color_dict == None:
        nucl_color_dict = color_dict  

    if annotation_loc is None:
        if with_seq == True:
            annotation_loc = "top"
        else:
            annotation_loc = "both" 

    if end == None:
        end = len(brick.seq) 

    if end < start:
        end = end + len(brick.seq) 
    
    if width == None:
        width = abs(end-start)  

    ceil = 0   
    fig       = plt.figure(figsize=(3,0.35))
    axes_list = [] 
    
    visible   = 1
    unvisible = 1 
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
    
    misc_color_count    = 0 
    label_color_dict    = {} 
    feature_color_count = collections.defaultdict(int) 
    for i, feat in enumerate(brick.dnafeature):
        if feat.type in unvisible_types:
            pass 
        
        elif visible == 1 or feat.type in visible_types:
            flag = 0 
            if "label" in feat.qualifiers:
                if type(feat.qualifiers["label"]) == list:
                    label = feat.qualifiers["label"][0]
                else:
                    label = feat.qualifiers["label"]
                flag = 1
            else:
                label = feat.type

            if "facecolor_dna.py" not in feat.qualifiers and "edgecolor_dna.py" not in feat.qualifiers:
                if label in label_color_dict:
                    feat.qualifiers["edgecolor_dna.py"] = [label_color_dict[label][1]]
                    feat.qualifiers["facecolor_dna.py"] = [label_color_dict[label][0]] 
                
                else:
                    cflag = 0 
                    for _type in feature_color_dict:
                        if feat.type == _type:
                            feat.qualifiers["edgecolor_dna.py"] = [feature_color_dict[_type][feature_color_count[_type]%len(feature_color_dict[_type])][1]]
                            feat.qualifiers["facecolor_dna.py"] = [feature_color_dict[_type][feature_color_count[_type]%len(feature_color_dict[_type])][0]]
                            feature_color_count[_type] += 1 
                            cflag = 1
                            break
                        else:
                            pass
                    if cflag == 0:
                        feat.qualifiers["edgecolor_dna.py"] = [misc_colors[misc_color_count%len(misc_colors)][1]]
                        feat.qualifiers["facecolor_dna.py"] = [misc_colors[misc_color_count%len(misc_colors)][0]] 
                        misc_color_count += 1 
                
                if flag == 1:
                    label_color_dict[label] = (feat.qualifiers["facecolor_dna.py"][0], feat.qualifiers["edgecolor_dna.py"][0]) 

    for num, sub_start in enumerate(list(range(start, end, width))):
        sub_end = sub_start + width
        if sub_end >= end:
            sub_end = end
        
        if sub_start >= len(brick.seq):
            sub_start = sub_start - len(brick.seq) 

        if sub_end > len(brick.seq):
            sub_end = sub_end - len(brick.seq) 
        
        sub_brick  = brick[sub_start:sub_end]        
        matplotlib.rcParams['font.size'] = fontsize
        if width < 101:
            std = 25
            head_length = 2
        elif width < 251:
            std = 50 
            head_length = 3
        elif width < 501:
            std = 100
            head_length = 6
        elif width < 1001:
            std = 250 
            head_length = 12
        elif width < 2001:
            std = 500
            head_length = 24
        elif width < 4001:
            std = 900
            head_length = 48
        elif width < 5001:
            std = 1200
            head_length = 56
        elif width < 6001:
            std = 1600
            head_length = 64
        else:
            std = 2000
            head_length = 72

        zero_position = sub_start + 1
        ax  = fig.add_axes([0, 0, enlarge_w*len(sub_brick.seq)/std, 1.0*enlarge_h], label=str(num)) 
        ax, y_list, ty_list = map_feat(fig, ax, sub_brick.dnafeature, len(sub_brick.seq), head_length, unvisible_types=unvisible_types, visible_types=visible_types, enlarge_w=enlarge_w, enlarge_h=enlarge_h, annotation_loc=annotation_loc)
        
        ty_list.append(0)  
        y_list.append(0) 
        
        ax.spines["top"].set_visible(False) 
        ax.spines["right"].set_visible(False) 
        ax.spines["left"].set_visible(False) 
        ax.spines["bottom"].set_position(("data",-1.0 * max(y_list)-0.5))
       
        y_list = list(-1 * np.array(y_list)) 
        ytop    = max(ty_list+y_list)+0.5
        ybottom = min(ty_list+y_list)-0.5

        ax.set_ylim(ybottom, ytop) 
        ax.set_yticks([])
       
        positions   = [p+zero_position-len(brick.seq) if p+zero_position > len(brick.seq) else p+zero_position for p in range(0, len(sub_brick.seq))] 
        ticks       = [x for p, x in zip(positions, list(range(0, len(sub_brick.seq)))) if p % int(std/(2*enlarge_w)) == 0 or p == 1] 
        tick_labels = [str(p) if p != 1 else "" for p, x in zip(positions, list(range(0, len(sub_brick.seq)))) if p % int(std/(2*enlarge_w)) == 0 or  p == 1]
        if with_seq == True: 
            ax.set_xticks([]) 
            ax.set_xticklabels([])
            ax.set_xlim(0, len(sub_brick.seq))
            ax.spines["bottom"].set_visible(False)  
            ax_seq = fig.add_axes([0, -0.65*enlarge_h/(abs(ytop-ybottom)), enlarge_w*len(sub_brick.seq)/std, 0.6*enlarge_h/(abs(ytop-ybottom))], label="{}_seq".format(num))
            
            if nucl_char != True and nucl_char != False:
                if width > 100:
                    colorbar(ax_seq, nucl_color_dict, sub_brick.seq, char=False)
                else:
                    colorbar(ax_seq, nucl_color_dict, sub_brick.seq, char=True)
            else:
                colorbar(ax_seq, nucl_color_dict, sub_brick.seq, char=nucl_char)

            ticks = np.array(ticks) + 0.5
            ax_seq.set_xticks(ticks) 
            ax_seq.set_xticklabels(tick_labels)
            ax_seq.set_xlim(0, len(sub_brick.seq)) 
            if num == 0:
                ax.set_position([0, 0, enlarge_w*len(sub_brick.seq)/std, 1.0*enlarge_h*(abs(ytop-ybottom))]) 
                ax_seq.set_position([0, -0.65*enlarge_h, enlarge_w*len(sub_brick.seq)/std, 0.6*enlarge_h]) 
                ceil = -0.65*enlarge_h
            else:
                ax.set_position([0, ceil-1.2*enlarge_h-1.0*enlarge_h*(abs(ytop-ybottom)), enlarge_w*len(sub_brick.seq)/std, 1.0*enlarge_h*(abs(ytop-ybottom))]) 
                ax_seq.set_position([0, ceil-1.2*enlarge_h-1.0*enlarge_h*(abs(ytop-ybottom))-0.65*enlarge_h, enlarge_w*len(sub_brick.seq)/std, 0.6*enlarge_h]) 
                ceil = ceil-1.2*enlarge_h-1.0*enlarge_h*(abs(ytop-ybottom))-0.65*enlarge_h
        else:
            ax.set_xticks(ticks) 
            ax.set_xticklabels(tick_labels)
            ax.set_xlim(0, len(sub_brick.seq))
            ax_seq = None
            if num == 0:
                ax.set_position([0, 0, enlarge_w*len(sub_brick.seq)/std, 1.0*enlarge_h*(abs(ytop-ybottom))]) 
                ceil = 0
            else:
                ax.set_position([0, ceil-1.2*enlarge_h-1.0*enlarge_h*(abs(ytop-ybottom)), enlarge_w*len(sub_brick.seq)/std, 1.0*enlarge_h*(abs(ytop-ybottom))]) 
                ceil = ceil-1.0*enlarge_h-1.2*enlarge_h*(abs(ytop-ybottom))
        
        axes_list.append((ax, ax_seq)) 
    #fig.set_size_inches(3, 0.35*(abs(ytop-ybottom)))
    return fig, axes_list  
