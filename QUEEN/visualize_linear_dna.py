import os 
import sys
import copy
import collections
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from  Bio import SeqIO
from matplotlib.transforms import Bbox

#sys.path.append("/".join(__file__.split("/")[:-1]))
#from queen import * 
matplotlib.rcParams["figure.max_open_warning"] = 0
matplotlib.rcParams['ps.fonttype']  = 42
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif']   = ["Arial","Lucida Sans","DejaVu Sans","Lucida Grande","Verdana"]
matplotlib.rcParams['font.family']       = 'sans-serif'
matplotlib.rcParams['font.size']         = 12.0
matplotlib.rcParams['font.weight']       = 500
matplotlib.rcParams["axes.labelcolor"]   = "#000000"
matplotlib.rcParams["axes.linewidth"]    = 1.0
matplotlib.rcParams["xtick.major.width"] = 1.0
matplotlib.rcParams["ytick.major.width"] = 1.0
matplotlib.rcParams['xtick.major.pad']   = 4
matplotlib.rcParams['ytick.major.pad']   = 4
matplotlib.rcParams['xtick.major.size']  = 4
matplotlib.rcParams['ytick.major.size']  = 4

color_dict = {"G":"#f2f059", "C":"#74b2d7", "A":"#79E5B7", "T":"#ff776c", "N":"#FFFFFF", "-":"#FFFFFF"}

#colorblind (facecolor_set)
colorblind=["#0173B2", "#DE8F05", "#029E73", "#D55E00", "#CC78BC", "#CA9161", "#FBAFE4", "#949494", "#ECE133", "#56B4E9"]

#pastel (edgecolor_set) 
pastel=["#A1C9F4", "#FFB482", "#8DE5A1", "#FF9F9B", "#D0BBFF", "#DEBB9B", "#FAB0E4", "#CFCFCF", "#FFFEA3", "#B9F2F0"]

for i in range(len(pastel)):    
    fc = pastel[i]
    fc = matplotlib.colors.to_rgba(fc) 
    fc = list(fc) 
    fc = list(map(lambda x: 1-0.7+0.7*x, fc[0:3]))
    pastel[i] = matplotlib.colors.to_hex(fc) 

feature_color_dict = collections.defaultdict(list)
feature_color_dict["misc_feature"]  = list(zip(pastel, colorblind))
feature_color_dict["CDS"]           = list(zip(pastel, colorblind))
misc_colors                         = list(zip(pastel, colorblind))
#feature_color_dict["CDS"]          = [('#ff9d9a', '#e15759'), ('#f1ce63', '#b6992d'), ('#86bcb6', '#499894'), ('#fabfd2', '#D37295'), ("#A0CBE8", "#4E79A7")]
#feature_color_dict["CDS"]          = [('#ff9f9a', '#c44e52'), ('#92c6ff', '#4c72b0'), ('#97f0aa', '#55a868'), ('#d0bbff', '#8172b2'), ('#fffea3', '#ccb974')]
#feature_color_dict["primer_bind"]  = [('#a6cee3', '#1f78b4'), ('#b2df8a', '#33a02c'), ('#fb9a99', '#e31a1c'), ('#fdbf6f', '#ff7f00'), ('#cab2d6', '#6a3d9a')] 
#feature_color_dict["primer"]       = [('#a6cee3', '#1f78b4'), ('#b2df8a', '#33a02c'), ('#fb9a99', '#e31a1c'), ('#fdbf6f', '#ff7f00'), ('#cab2d6', '#6a3d9a')] 
#feature_color_dict["promoter"]     = [('#b0e0e6', '#64b5cd'), ('#92c6ff', '#4c72b0'), ('#97f0aa', '#55a868')]  
#feature_color_dict["rep_origin"]   = [('#fff2ae', '#ffd92f')]
#feature_color_dict["misc_feature"] = [('#fbb4ae', '#e41a1c'), ('#b3cde3', '#377eb8'), ('#ccebc5', '#4daf4a'), ('#decbe4', '#984ea3'), ('#fed9a6', '#ff7f00')]  
#misc_colors = [('#ffffcc', '#d9d927'), ('#e5d8bd', '#a65628'), ('#fddaec', '#f781bf'), ('#f2f2f2', '#999999'), ('#fbb4ae', '#e41a1c')]

def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i+lv//3], 16) for i in range(0, lv, lv//3))

def map_feat(fig, ax, feats, length, head_length, unvisible_types=["source"], visible_types=[], enlarge_w=1.0, enlarge_h=1.0, annotation_loc="either", label_visible=True, fontsize=12, project="", title_visible=True, axis_visible=True, tick_space="auto", labelcolor="k", titlename=None):
    y = 0
    y_list      = [] 
    ty_list     = [] 
    visible     = 1
    unvisible   = 1 
    head_length = head_length * 1.0/enlarge_w 
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
            
            wd1 = 0.82
            if abs(ge-gs) < head_length * 1.2:
                hl  = abs(ge-gs)
                margin1 = head_length * 0.15 
                margin2 = head_length * 0.25 
                hl2 = hl - (margin1+margin2)
                if hl2 < 0:
                    hl2 = 0
                    wd2 = 0 
                else:
                    wd2 = 0.8 * ((hl2/hl) ** 1.0) * wd1/0.8
            else:
                hl  = head_length 
                hl2 = hl*0.85  
                margin1 = hl * 0.15 
                margin2 = hl * 0.25 
                wd2     = 0.65 * wd1/0.8

            if "facecolor_queen" in feat.qualifiers:
                if type(feat.qualifiers["facecolor_queen"]) == list:
                    facecolor = feat.qualifiers["facecolor_queen"][0] 
                else:
                    facecolor = feat.qualifiers["facecolor_queen"]
            else:
                facecolor = "#ffffec" 
            
            if "edgecolor_queen" in feat.qualifiers:
                if type(feat.qualifiers["edgecolor_queen"]) == list:
                    edgecolor = feat.qualifiers["edgecolor_queen"][0] 
                else:
                    edgecolor = feat.qualifiers["edgecolor_queen"]
            else:
                if strand == 1:
                    edgecolor = "#FACAC8" 
                elif strand == -1:
                    edgecolor = "#C8DFFA" 
                else:
                    edgecolor = "#CCCCCC"
            
            if "labelcolor_queen" in feat.qualifiers:
                if type(feat.qualifiers["labelcolor_queen"]) == list:
                    lc = feat.qualifiers["labelcolor_queen"][0] 
                else:
                    lc = feat.qualifiers["labelcolor_queen"]
            else:
                lc = None
            
            if "labelweight_queen" in feat.qualifiers:
                if type(feat.qualifiers["labelweight_queen"]) == list:
                    fw = feat.qualifiers["labelweight_queen"][0] 
                else:
                    fw = feat.qualifiers["labelweight_queen"]
            else:
                fw = None

            if "broken_feature" in feat.qualifiers:
                note   = feat.qualifiers["broken_feature"]
                if type(note) == list:
                    note   = feat.qualifiers["broken_feature"][0]
                feat_length = int(note.split(":")[-4]) 
                pos_s  = int(note.split(":")[-1].split("..")[0]) 
                pos_e  = int(note.split(":")[-1].split("..")[1])
                #if (pos_s != 1 or pos_e != feat_length):
                #    label = note.split(":")[1]

            if strand == 1:
                if "broken_feature" in feat.qualifiers:
                    if (gs == 0 or ge >= length) and (pos_s != 1 or pos_e != feat_length):
                        if gs == 0 and ge >= length:
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=wd1, head_width=wd1, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                        elif gs == 0:
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=wd1, head_width=wd1, head_length=hl, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs-margin2, dy=0, width=wd2, head_width=wd2, head_length=hl2, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                        else:
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=wd1, head_width=wd1, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=gs+margin1, y=-1*y, dx=ge-gs-margin1, dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                    else:
                        ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=wd1, head_width=wd1, head_length=hl, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                        ax.arrow(x=gs+margin1, y=-1*y, dx=ge-gs-(margin1+margin2), dy=0, width=wd2, head_width=wd2, head_length=hl2, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                else:
                    ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=wd1, head_width=wd1, head_length=hl, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                    ax.arrow(x=gs+margin1, y=-1*y, dx=ge-gs-(margin1+margin2), dy=0, width=wd2, head_width=wd2, head_length=hl2, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
            
            elif strand == -1:
                if "broken_feature" in feat.qualifiers:
                    if (gs == 0 or ge >= length) and (pos_s != 1 or pos_e != feat_length):
                        if gs == 0 and ge >= length:
                            ax.arrow(x=ge, y=-1*y, dx=gs-ge, dy=0, width=wd1, head_width=wd1, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=ge, y=-1*y, dx=gs-ge, dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                        elif gs == 0:
                            ax.arrow(x=ge, y=-1*y, dx=gs-ge, dy=0, width=wd1, head_width=wd1, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=ge-margin1, y=-1*y, dx=gs-ge+margin1, dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                        else:
                            ax.arrow(x=ge, y=-1*y, dx=gs-ge, dy=0, width=wd1, head_width=wd1, head_length=hl, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=ge, y=-1*y, dx=gs-ge+margin2, dy=0, width=wd2, head_width=wd2, head_length=hl2, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                    else:
                        ax.arrow(x=ge, y=-1*y, dx=gs-ge, dy=0, width=wd1, head_width=wd1, head_length=hl, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                        ax.arrow(x=ge-margin1, y=-1*y, dx=gs-ge+(margin1+margin2), dy=0, width=wd2, head_width=wd2, head_length=hl2, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                else:
                    ax.arrow(x=ge, y=-1*y, dx=gs-ge, dy=0, width=wd1, head_width=wd1, head_length=hl, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                    ax.arrow(x=ge-margin1, y=-1*y, dx=gs-ge+(margin1+margin2), dy=0, width=wd2, head_width=wd2, head_length=hl2, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
            
            else:
                hl = 0 
                if "broken_feature" in feat.qualifiers:
                    if (gs == 0 or ge >= length) and (pos_s != 1 or pos_e != feat_length):
                        if gs == 0 and ge >= length:
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=wd1, head_width=wd1, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                        elif gs == 0:
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=wd1, head_width=wd1, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs-margin2, dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                        else:
                            ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=wd1, head_width=wd1, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                            ax.arrow(x=gs+margin1, y=-1*y, dx=ge-gs-margin1, dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                    else:
                        ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=wd1, head_width=wd1, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                        ax.arrow(x=gs+margin1, y=-1*y, dx=ge-gs-(margin1+margin2), dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
                else: 
                    ax.arrow(x=gs, y=-1*y, dx=ge-gs, dy=0, width=wd1, head_width=wd1, head_length=0, length_includes_head=True, color='k', fc=edgecolor, lw=0.0)
                    ax.arrow(x=gs+margin1, y=-1*y, dx=ge-gs-(margin1+margin2), dy=0, width=wd2, head_width=wd2, head_length=0, length_includes_head=True, color='k', fc=facecolor, lw=0.0)
            if strand  == -1:   
                label_position_list.append((label,(gs+ge+hl*0.5)/2, -1*y, gs, ge, hl, edgecolor, facecolor, lc, fw))
            else:
                label_position_list.append((label,(gs+ge-hl*0.5)/2, -1*y, gs, ge, hl, edgecolor, facecolor, lc, fw))
            for j in range(gs,ge):
                gene_position_matrix[y][j] = 1    
            y_list.append(y)
    
    renderer     = fig.canvas.get_renderer()
    coordinate  = ax.transData.inverted() 
    for label, tx, ty, gs, ge, hl, ec, fc, lc, fw in label_position_list:
        #fc = "#ffffec" #facecolor
        if lc is None:
            lc = labelcolor
        else:
            pass

        if fw is None:
            fw="bold" if lc=="w" else None
        else:
            pass 

        text = ax.text(tx, ty-0.05, label, ha="center", va="center", fontsize=fontsize, color=lc, fontweight=fw)
        bbox_text = text.get_window_extent(renderer=renderer)
        bbox_text = Bbox(coordinate.transform(bbox_text))
        width = bbox_text.width*length
        ts    = int(tx-int(width/2)*1.2) 
        ts    = 0 if ts < 0 else ts 
        te    = int(tx+int(width/2)*1.2)
        te    = length if te > length else te
        if width > abs(ge-gs) * 0.90:
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
            
            if label_visible == 2:
                if annotation_loc == "either": 
                    if t % 2 == 0:
                        ax.text(tx, t//2+1, label, ha="center",va="center", bbox=dict(boxstyle="round", ec=ec, fc=fc, pad=0.3, lw=1.4), zorder=100, fontsize=fontsize, color=lc, fontweight=fw)
                        ax.plot([(gs+ge)/2, (gs+ge)/2], [t//2+1, ty], lw=0.5, color="k", zorder=0)
                        ty_list.append(t//2+1)
                    
                    else:
                        if axis_visible == True:
                            bufspace = 0.8 * fontsize/12
                        else:
                            bufspace = 0 
                        ax.text(tx, -1*max(y_list)-1-bufspace-t//2, label, ha="center",va="center",  bbox=dict(boxstyle="round", ec=ec, fc=fc, pad=0.3, lw=1.4), zorder=100, fontsize=fontsize, color=lc, fontweight=fw)
                        ax.plot([(gs+ge)/2, (gs+ge)/2], [-1*max(y_list)-1-bufspace-t//2, ty], lw=0.5, color="k", zorder=0)
                        ty_list.append(-1*max(y_list)-1-bufspace-t//2)
                
                elif annotation_loc == "bottom":
                    if axis_visible == True:
                        bufspace = 0.8 * fontsize/12
                    else:
                        bufspace = 0 
                    ax.text(tx, -1*max(y_list)-1-bufspace-t, label, ha="center",va="center",  bbox=dict(boxstyle="round", ec=ec, fc=fc, pad=0.3, lw=1.4), zorder=100, fontsize=fontsize, color=lc, fontweight=fw)
                    ax.plot([(gs+ge)/2, (gs+ge)/2], [-1*max(y_list)-1-bufspace-t, ty], lw=0.5, color="k", zorder=0)
                    ty_list.append(-1*max(y_list)-1-bufspace-t)
                    
                else:
                    ax.text(tx, t+1, label, ha="center",va="center", bbox=dict(boxstyle="round", ec=ec, fc=fc, pad=0.3, lw=1.4), zorder=100, fontsize=fontsize, color=lc, fontweight=fw)
                    ax.plot([(gs+ge)/2, (gs+ge)/2], [t+1, ty], lw=0.5, color="k", zorder=0)
                    ty_list.append(t+1)

                for j in range(ts, te):
                    text_position_matrix[t][j] = 1
        else:
            if label_visible == 0:
                text.set_visible(False)
            ty_list.append(ty) 
    return ax, y_list, ty_list

def colorbar(ax, color_dict, ref_seq, char=False, fontsize=10):
    if len(ref_seq) > 200:
        bars = ax.bar(list(range(len(ref_seq))), [0.9] * (len(ref_seq)), width=1.0, edgecolor="#BBBBBB", linewidth=0.0, align="edge",bottom=0.05)
    else:
        bars = ax.bar(list(range(len(ref_seq))), [0.9] * (len(ref_seq)), width=1.0, edgecolor="#BBBBBB", linewidth=0.2, align="edge",bottom=0.05)
    ax.set_xlim(0,len(ref_seq))
    ax.set_ylim(0,1.00)
    p = 0
    for bar, c in zip(bars,ref_seq):
        color = color_dict[c]
        bar.set_facecolor(color)
        if char == True:
            bar.set_alpha(0.7)
            ax.text(p+0.5,0.45,c,va="center",ha="center",fontsize=fontsize,zorder=100)     
        p += 1 
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.patch.set_alpha(0.0)
    return bars

def visualize(brick, start=0, end=None, wrap_width=None, annotation_loc=None, label_visible=True, feature_list=None, unvisible_types=["source"], visible_types=[], enlarge_w=1.0, enlarge_h=1.0, scale="fix", fontsize=12, with_seq=False, nucl_char=None, nucl_color_dict=None, title_visible=True, axis_visible=True, tick_space="auto", labelcolor="k", titlename=None, fig=None):
    if titlename is None:
        titlename = brick.project
    else:
        pass
    brick = copy.deepcopy(brick) 
    width = wrap_width 
    if nucl_color_dict == None:
        nucl_color_dict = color_dict  

    if annotation_loc is None:
        if with_seq == True:
            annotation_loc = "top"
        else:
            annotation_loc = "either" 

    if end == None:
        end = len(brick.seq) 

    if end < start:
        end = end + len(brick.seq) 
    
    if width == None:
        width = abs(end-start)  

    ceil = 0   
    if fig is None:
        fig = plt.figure(figsize=(3,0.30))
        basenum = 0
    else:
        fig = fig
        fig.set_size_inches(3,0.30) 
        basenum = len(fig.axes)+1 

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
    
    if feature_list is None:
        feature_list = brick.dnafeatures    
    else:
        brick._dnafeatures = copy.deepcopy(feature_list)
        brick._features_dict = dict(list(map(lambda x:(x._id, x), brick.dnafeatures)))

    for i, feat in enumerate(brick.dnafeatures):
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

            if "facecolor_queen" not in feat.qualifiers and "edgecolor_queen" not in feat.qualifiers:
                if label in label_color_dict:
                    feat.qualifiers["edgecolor_queen"] = [label_color_dict[label][1]]
                    feat.qualifiers["facecolor_queen"] = [label_color_dict[label][0]] 
                
                else:
                    cflag = 0 
                    for _type in feature_color_dict:
                        if feat.type == _type:
                            feat.qualifiers["edgecolor_queen"] = [feature_color_dict[_type][feature_color_count[_type]%len(feature_color_dict[_type])][1]]
                            feat.qualifiers["facecolor_queen"] = [feature_color_dict[_type][feature_color_count[_type]%len(feature_color_dict[_type])][0]]
                            feature_color_count[_type] += 1 
                            cflag = 1
                            break
                        else:
                            pass
                    if cflag == 0:
                        feat.qualifiers["edgecolor_queen"] = [misc_colors[misc_color_count%len(misc_colors)][1]]
                        feat.qualifiers["facecolor_queen"] = [misc_colors[misc_color_count%len(misc_colors)][0]] 
                        misc_color_count += 1 
                
                if flag == 1:
                    label_color_dict[label] = (feat.qualifiers["facecolor_queen"][0], feat.qualifiers["edgecolor_queen"][0]) 
            
            if "label" in feat.qualifiers:
                if type(feat.qualifiers["label"]) == list:
                    label = feat.qualifiers["label"][0]
                else:
                    label = feat.qualifiers["label"]
            else:
                label = feat.type

            if "broken_feature" in feat.qualifiers:
                posinfo     = feat.qualifiers["broken_feature"][0].split(":")[-1]
                feat_length = feat.qualifiers["broken_feature"][0].split(":")[-4] 
                label = label + ":" + posinfo + ":" +  feat_length
                feat.qualifiers["label"] = [label] 

    for num, sub_start in enumerate(list(range(start, end, width))):
        sub_end = sub_start + width
        if sub_end >= end:
            sub_end = end
        
        if sub_start >= len(brick.seq):
            sub_start = sub_start - len(brick.seq) 

        if sub_end > len(brick.seq):
            sub_end = sub_end - len(brick.seq) 
        
        sub_brick  = brick[sub_start:sub_end]
        matplotlib.rcParams['font.size'] = fontsize if fontsize >=12 else 12
        """
        if scale == "fix":
            std = 500
            head_length = 24
        else:
            if scale == "auto":
                scale = width
            else:
                pass 
            if scale < 101:
                std = 25
                head_length = 1.5
            elif scale < 251:
                std = 50 
                head_length = 3
            elif scale < 501:
                std = 100
                head_length = 6
            elif scale < 1001:
                std = 200 
                head_length = 12
            elif scale < 2001:
                std = 500
                head_length = 24
            elif scale < 4001:
                std = 1000
                head_length = 48
            elif scale < 5001:
                std = 1200
                head_length = 56
            elif scale < 6001:
                std = 1600
                head_length = 64
            else:
                std = 2000
                head_length = 72
        """
        if scale == "fix":
            std = 1000
            head_length = 48
        #0.003 inch per nucleotide 
        if enlarge_w == "auto":
            if with_seq == True:
                enlarge_w = 40
            else:  
                enlarge_w = 1.0

        if enlarge_h == "auto":
            enlarge_h = 1.0
        zero_position = sub_start + 1
        ax  = fig.add_axes([0, 0, enlarge_w*len(sub_brick.seq)/std, 1.0*enlarge_h], label=str(num+basenum)) 
        ax, y_list, ty_list = map_feat(fig, ax, sub_brick.dnafeatures, len(sub_brick.seq), head_length, unvisible_types=unvisible_types, visible_types=visible_types, enlarge_w=enlarge_w, enlarge_h=enlarge_h, annotation_loc=annotation_loc, label_visible=label_visible, fontsize=fontsize, project=sub_brick.project, title_visible=title_visible, axis_visible=axis_visible, labelcolor=labelcolor, titlename=titlename)
        
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
        if tick_space == "auto":
            ticks       = [x for p, x in zip(positions, list(range(0, len(sub_brick.seq)))) if p % int(std/(2*enlarge_w)) == 0 or p == 1] 
            tick_labels = [str(p) if p != 1 else "1" for p, x in zip(positions, list(range(0, len(sub_brick.seq)))) if p % int(std/(2*enlarge_w)) == 0 or  p == 1]
        elif tick_space > 0:
            ticks       = [x for p, x in zip(positions, list(range(0, len(sub_brick.seq)))) if p % tick_space == 0 or p == 1] 
            tick_labels = [str(p) if p != 1 else "1" for p, x in zip(positions, list(range(0, len(sub_brick.seq)))) if p % tick_space == 0 or  p == 1]
        else:
            ticks = [] 
            tick_labels = [] 

        if with_seq == True: 
            ax.set_xticks([]) 
            ax.set_xticklabels([])
            ax.set_xlim(0, len(sub_brick.seq))
            ax.spines["bottom"].set_visible(False)
            ax_seq = fig.add_axes([0, -0.60*enlarge_h/(abs(ytop-ybottom)), enlarge_w*len(sub_brick.seq)/std, 0.6*enlarge_h/(abs(ytop-ybottom))], label="{}_seq".format(num+basenum))
            if nucl_char != True and nucl_char != False:
                if enlarge_w < 30:
                    colorbar(ax_seq, nucl_color_dict, sub_brick.seq, char=False, fontsize=fontsize-2)
                else:
                    colorbar(ax_seq, nucl_color_dict, sub_brick.seq, char=True, fontsize=fontsize-2)
            else:
                colorbar(ax_seq, nucl_color_dict, sub_brick.seq, char=nucl_char)

            ticks = np.array(ticks) + 0.5
            ax_seq.set_xticks(ticks) 
            ax_seq.set_xticklabels(tick_labels)
            ax_seq.set_xlim(0, len(sub_brick.seq)) 
            if num == 0:
                if title_visible == True:
                    ax.set_title(titlename, fontsize=fontsize * 1.125, pad=15)
                ax.set_position([0, 0, enlarge_w*len(sub_brick.seq)/std, 1.0*enlarge_h*(abs(ytop-ybottom))]) 
                ax_seq.set_position([0, -0.65*enlarge_h, enlarge_w*len(sub_brick.seq)/std, 0.6*enlarge_h]) 
                ceil = -0.65*enlarge_h
            else:
                ax.set_position([0, ceil-1.2*enlarge_h-1.0*enlarge_h*(abs(ytop-ybottom)), enlarge_w*len(sub_brick.seq)/std, 1.0*enlarge_h*(abs(ytop-ybottom))]) 
                ax_seq.set_position([0, ceil-1.2*enlarge_h-1.0*enlarge_h*(abs(ytop-ybottom))-0.65*enlarge_h, enlarge_w*len(sub_brick.seq)/std, 0.6*enlarge_h]) 
                ceil = ceil-1.2*enlarge_h-1.0*enlarge_h*(abs(ytop-ybottom))-0.65*enlarge_h
            if axis_visible == False:
                ax_seq.spines["bottom"].set_visible(False) 
                ax_seq.set_xticks([]) 
        else:
            ax.set_xticks(ticks) 
            ax.set_xticklabels(tick_labels)
            ax.set_xlim(0, len(sub_brick.seq))
            ax_seq = None
            if num == 0:
                if title_visible == True:
                    ax.set_title(titlename, fontsize=fontsize * 1.125)
                ax.set_position([0, 0, enlarge_w*len(sub_brick.seq)/std, 1.0*enlarge_h*(abs(ytop-ybottom))]) 
                ceil = 0
            else:
                ax.set_position([0, ceil-1.2*enlarge_h-1.0*enlarge_h*(abs(ytop-ybottom)), enlarge_w*len(sub_brick.seq)/std, 1.0*enlarge_h*(abs(ytop-ybottom))]) 
                ceil = ceil-1.0*enlarge_h-1.2*enlarge_h*(abs(ytop-ybottom))  
            
            if axis_visible == False:
                ax.spines["bottom"].set_visible(False) 
                ax.set_xticks([])
        
        axes_list.append((ax, ax_seq)) 
    #fig.set_size_inches(3, 0.35*(abs(ytop-ybottom)))
    return fig, axes_list  
