import os 
import sys
import copy
import collections
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from  Bio import SeqIO
from matplotlib.transforms import Bbox

matplotlib.rcParams["figure.max_open_warning"] = 0
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.sans-serif']   = ["Arial","Lucida Sans","DejaVu Sans","Lucida Grande","Verdana"]
matplotlib.rcParams['font.family']       = 'sans-serif'
matplotlib.rcParams['font.size']         = 8.0
matplotlib.rcParams['font.weight']       = 500
matplotlib.rcParams["axes.labelcolor"]   = "#000000"
matplotlib.rcParams["axes.linewidth"]    = 1.0
matplotlib.rcParams["xtick.major.width"] = 1.0
matplotlib.rcParams["ytick.major.width"] = 1.0
matplotlib.rcParams['xtick.major.pad']   = 6
matplotlib.rcParams['ytick.major.pad']   = 6
matplotlib.rcParams['xtick.major.size']  = 6
matplotlib.rcParams['ytick.major.size']  = 6

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

def map_feat(fig, ax, ax2, feats, length, head_length=np.pi * 0.030, unvisible_types=["source"], visible_types=[], enlarge=1.0, format=1, bottom=300, fontsize=10, label_visible=True, axis_visible=True, tick_space="auto", labelcolor="k"):
    if format == 0 or format == 1:
        outer    = 62 * 1.2
        inner    = 45 * 1.2
        if bottom is None and format == 1:
            bottom_h = 600
        elif bottom is None and format == 0:
            bottom_h = 300
        else:
            bottom_h = bottom 
        lane_h   = 70 * 1.2
        normal_w = 1000
        matplotlib.rcParams['font.size'] = 8.0
    
    else:
        outer    = 60
        inner    = 44
        if bottom is None:
            bottom_h = 400
        else:
            bottom_h = bottom 
        bottom_h = bottom
        lane_h   = 250
        normal_w = 1000
        matplotlib.rcParams['font.size'] = 7.0

    fig_width  = normal_w
    renderer   = fig.canvas.get_renderer()
    coordinate = ax2.transData.inverted() 
    y_list     = [] 
    ty_list    = [] 
    visible    = 1
    unvisible  = 1 
    gene_position_matrix = [[0] * length]
    text_position_matrix = [np.array([0] * length)] + [np.array([0] * length)] + [np.array([0] * length)]
    label_position_list  = [] 
    
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
    for i, feat in enumerate(feats):
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
            
            if label == "":
                label = " "

            strand = feat.location.strand
            if strand == -1:
                gs_origin = feat.location.parts[-1].start.position 
                ge_origin = feat.location.parts[0].end.position
                gs = gs_origin * 2 * np.pi / length
                ge = ge_origin * 2 * np.pi / length 
            else:
                gs_origin = feat.location.parts[0].start.position 
                ge_origin = feat.location.parts[-1].end.position
                gs = gs_origin * 2 * np.pi / length 
                ge = ge_origin * 2 * np.pi / length

            
            y = 0
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
                        text_position_matrix.append(np.array([0] * length))
                        text_position_matrix.append(np.array([0] * length))
                        text_position_matrix.append(np.array([0] * length))
                    else:
                        pass 
                else:
                    for y, row in enumerate(gene_position_matrix):
                        if 1 in row[gs_origin:]:
                            flag1 = 1
                        else:
                            flag1 = 0
                            break 
                    
                    for y, row in enumerate(gene_position_matrix):
                        if 1 in row[:ge_origin]:
                            flag2 = 1
                        else:
                            flag2 = 0
                            break 

                    if flag1 == 1 or flag2 == 1:
                        y += 1 
                        gene_position_matrix.append([0] * length)
                        text_position_matrix.append(np.array([0] * length))
                        text_position_matrix.append(np.array([0] * length))
                        text_position_matrix.append(np.array([0] * length))
                    else:
                        pass 
            
            
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
                

            if gs < ge: 
                width  = ge-gs
                middle = (ge+gs)/2
            else:
                width  = ge + (2*np.pi-gs) 
                middle = (ge+gs+2*np.pi)/2  
            
            if format == 1 or format == 0:
                margin = np.pi*(0.006)

            else:
                margin = np.pi*(0.005)
            
            if abs(ge-gs) < head_length * 1.2:
                hl  = abs(ge-gs)
            else:
                hl  = head_length * bottom_h/(y*lane_h+bottom_h)
            mg  = margin * bottom_h/(y*lane_h+bottom_h)

            if format == 0: 
                w, x, fc, ec = width, middle, facecolor, edgecolor
                if gs < ge: 
                    gmiddle = (ge+gs)/2
                else:
                    gmiddle = (ge+gs+2*np.pi)/2  
                slide = 0 
                pos_list   = []  
                width_list = []
                for char in label:
                    text        = ax2.text(slide, 0, char, ha="right", va="center", fontsize=fontsize)
                    bbox_text   = text.get_window_extent(renderer=renderer)
                    bbox_text   = Bbox(coordinate.transform(bbox_text))
                    text.set_visible(False)
                    width_list.append(bbox_text.width) 
                    pos_list.append(slide+bbox_text.width/2) 
                    slide += bbox_text.width
                
                pos_list = [-1*(p-0.5*slide) * 2*fig_width for p in pos_list]
                new_pos_list = [] 
                for pos, btw in zip(pos_list, width_list):
                    new_pos_list.append((np.arccos(pos/(y*lane_h+bottom_h))-0.5*np.pi+x, y*lane_h+bottom_h, y, pos, x, btw))
               
                tflag = 0 
                t_width = (new_pos_list[-1][0] - new_pos_list[0][0]) 
                if t_width < w-2*head_length:
                    pass 
                else:
                    shifted_pos_list = [] 
                    if feat.location.strand == -1:
                        tflag = -1
                        for pos in new_pos_list:
                            pos = list(pos)
                            new_pos    = pos
                            new_pos[0] = pos[0] - new_pos_list[0][0] + new_pos_list[0][-1] + ge
                            shifted_pos_list.append(new_pos) 
                        new_pos_list = shifted_pos_list
                    else:
                        tflag = 1
                        for pos in new_pos_list:
                            pos = list(pos)
                            new_pos    = pos
                            new_pos[0] = pos[0] - new_pos_list[-1][0] - new_pos_list[-1][-1] + gs
                            shifted_pos_list.append(new_pos) 
                        new_pos_list = shifted_pos_list

                if new_pos_list[len(new_pos_list) // 2][0] < 0.5 * np.pi or new_pos_list[len(new_pos_list) // 2][0] > 1.5 * np.pi:
                    rotation = lambda x:(-1.0*x)*180/np.pi 
                else:
                    label    = label[::-1]
                    slide = 0 
                    pos_list   = [] 
                    width_list = []
                    for char in label:
                        text        = ax2.text(slide, 0, char, ha="right", va="center", fontsize=fontsize)
                        bbox_text   = text.get_window_extent(renderer=renderer)
                        bbox_text   = Bbox(coordinate.transform(bbox_text))
                        text.set_visible(False)
                        pos_list.append(slide+bbox_text.width/2) 
                        width_list.append(bbox_text.width) 
                        slide += bbox_text.width
                    
                    pos_list = [-1*(p-0.5*slide) * 2 * fig_width for p in pos_list]
                    new_pos_list = [] 
                    for pos, btw in zip(pos_list, width_list):
                        new_pos_list.append((np.arccos(pos/(y*lane_h+bottom_h))-0.5*np.pi+x, y*lane_h+bottom_h, y, pos, x, btw))
                    
                    rotation = lambda x:(-1.0*x)*180/np.pi+180
                    t_width = (new_pos_list[-1][0] - new_pos_list[0][0]) 
                    if t_width < w-2*head_length:
                        pass 
                    else:
                        shifted_pos_list = [] 
                        if feat.location.strand == -1:
                            tflag = -1
                            for pos in new_pos_list:
                                pos = list(pos)
                                new_pos    = pos
                                new_pos[0] = pos[0] - new_pos_list[0][0] + new_pos_list[0][-1] + ge
                                shifted_pos_list.append(new_pos) 
                            new_pos_list = shifted_pos_list
                        else:
                            tflag =1
                            for pos in new_pos_list:
                                pos = list(pos)
                                new_pos    = pos
                                new_pos[0] = pos[0] - new_pos_list[-1][0] - new_pos_list[-1][-1] + gs
                                shifted_pos_list.append(new_pos) 
                            new_pos_list = shifted_pos_list

                new_origins = [pos * length / (2*np.pi) for pos in [p[0] for p in new_pos_list]] 
                new_gs = gs * length / (2*np.pi)  
                new_ge = ge * length / (2*np.pi) 
                if tflag == 0:
                    new_gs_origin = new_gs 
                    new_ge_origin = new_ge
                
                elif tflag == -1:
                    new_gs_origin = new_gs
                    new_ge_origin = new_origins[-1]
                
                elif tflag == 1:
                    new_gs_origin = new_origins[0] 
                    new_ge_origin = new_ge

                while new_gs_origin < 0:
                    new_gs_origin = new_gs_origin + length 

                while new_ge_origin > length:
                    new_ge_origin = new_ge_origin - length 
                
                gs_origin, ge_origin = int(new_gs_origin), int(new_ge_origin+1)
                flag = 0 
                if i > 0:
                    if gs_origin < ge_origin:
                        for yy, row in enumerate(gene_position_matrix[y:]):
                            if 1 in row[gs_origin:ge_origin]:
                                flag = 1
                            else:
                                flag = 0
                                break     
                        if flag == 1:
                            y = y + yy + 1
                            gene_position_matrix.append([0] * length)
                        else:
                            y = y + yy
                    else:
                        for yy, row in enumerate(gene_position_matrix[y:]):
                            if 1 in row[gs_origin:] or 1 in row[:ge_origin]:
                                flag = 1
                            else:
                                flag = 0

                        if flag == 1:
                            y = y + yy + 1
                            gene_position_matrix.append([0] * length)
                        else:
                            y = y + yy
                
                if abs(ge-gs) < head_length * 1.2:
                    hl  = abs(ge-gs)
                else:
                    hl  = head_length * bottom_h/(y*lane_h+bottom_h)
                mg = margin * bottom_h/(y*lane_h+bottom_h)
                       
                new_new_pos_list = []
                for pos_info in new_pos_list:
                    _, _, _, pos, x, btw = pos_info 
                    new_new_pos_list.append((np.arccos(pos/(y*lane_h+bottom_h))-0.5*np.pi+x, y*lane_h+bottom_h, y, pos, x, btw))
                new_pos_list = new_new_pos_list
                
                tflag = 0 
                t_width = (new_pos_list[-1][0] - new_pos_list[0][0]) 
                if t_width < w-1.5*head_length or label_visible <= 1:
                    if t_width < w-1.5*head_length:
                        tflag = 0
                    else:
                        pass 
                else:
                    shifted_pos_list = [] 
                    if feat.location.strand == -1:
                        tflag = -1
                        for pos in new_pos_list:
                            pos = list(pos)
                            new_pos    = pos
                            new_pos[0] = pos[0] - new_pos_list[0][0] + new_pos_list[0][-1] + 2*mg + ge
                            shifted_pos_list.append(new_pos) 
                        new_pos_list = shifted_pos_list
                    else:
                        tflag = 1
                        for pos in new_pos_list:
                            pos = list(pos)
                            new_pos    = pos
                            new_pos[0] = pos[0] - new_pos_list[-1][0] - new_pos_list[-1][-1] - 2*mg + gs
                            shifted_pos_list.append(new_pos) 
                        new_pos_list = shifted_pos_list

                new_origins   = [pos * length / (2*np.pi) for pos in ([p[0] for p in new_pos_list])] 
                if tflag == 0:
                    new_gs_origin = new_gs 
                    new_ge_origin = new_ge
                
                elif tflag == -1:
                    new_gs_origin = new_gs
                    new_ge_origin = new_origins[-1]
                
                elif tflag == 1:
                    new_gs_origin = new_origins[0] 
                    new_ge_origin = new_ge

                new_gs_origin = new_gs_origin - 2 * mg * length 
                new_ge_origin = new_ge_origin + 2 * mg * length
                
                while new_gs_origin < 0:
                    new_gs_origin = new_gs_origin + length 

                while new_ge_origin >= length:
                    new_ge_origin = new_ge_origin - length 
                
                gs_origin, ge_origin = int(new_gs_origin), int(new_ge_origin) 
                for char, (theta, height, y, pos, x, btw) in zip(label, new_pos_list):
                    if tflag == 1 or label_visible>=1:
                        ax.text(theta, height, char, ha="center", va="center", rotation=rotation(theta), zorder=10, fontsize=fontsize, color=labelcolor)
              
                if width > 1.2 * head_length:
                    if strand == 1:
                        ax.bar([gs], [outer], bottom=y*lane_h+bottom_h-outer/2, width=width-hl*0.98, align="edge", fc=edgecolor, lw=0.0, zorder=1)
                        ax.bar([gs+mg], [inner], bottom=y*lane_h+bottom_h-inner/2, width=width-mg-hl*0.98, align="edge", fc=facecolor, lw=0.0, zorder=2)
                        ax.arrow(x=ge-hl, y=y*lane_h+bottom_h, dx=hl*1.05, dy=0, width=outer, head_width=outer, head_length=hl, length_includes_head=True, fc=edgecolor, lw=0.0, zorder=3)
                        ax.arrow(x=ge-1.1*hl, y=y*lane_h+bottom_h, dx=1.1*hl-mg*1.4, dy=0, width=inner, head_width=inner, head_length=hl-mg*1.4, length_includes_head=True, fc=facecolor, lw=0.0, zorder=4)
                    elif strand == -1:
                        ax.bar([gs+hl*0.98], [outer], bottom=y*lane_h+bottom_h-outer/2, width=width-hl*0.98, align="edge", fc=edgecolor, lw=0.0, zorder=1)
                        ax.bar([gs+hl*0.98], [inner], bottom=y*lane_h+bottom_h-inner/2, width=width-mg-hl*0.98, align="edge", fc=facecolor, lw=0.0, zorder=2)
                        ax.arrow(x=gs+hl, y=y*lane_h+bottom_h, dx=-1*hl*1.05, dy=0, width=outer, head_width=outer, head_length=hl, length_includes_head=True, fc=edgecolor, lw=0.0, zorder=3) 
                        ax.arrow(x=gs+1.1*hl, y=y*lane_h+bottom_h, dx=-1*(1.1*hl-mg*1.4), dy=0, width=inner, head_width=inner, head_length=hl-mg*1.4, length_includes_head=True, fc=facecolor, lw=0.0, zorder=4)
                    else:
                        ax.bar([gs], [outer], bottom=y*lane_h+bottom_h-outer/2, width=width, align="edge", fc=edgecolor, lw=0.0, zorder=1)
                        ax.bar([gs+mg], [inner], bottom=y*lane_h+bottom_h-inner/2, width=width-2*mg, align="edge", fc=facecolor, lw=0.0, zorder=2)
                else:
                    if strand == 1:
                        ax.arrow(x=ge-width, y=y*lane_h+bottom_h, dx=width, dy=0, width=outer, head_width=outer, head_length=width, length_includes_head=True, fc=edgecolor, lw=0.0)
                        if width > 2.4 * mg:
                            hw = outer * 1.0 * (width-2.4*mg)/width
                            if hw > outer:
                                hw = outer * (width-2.4*mg)/width
                            ax.arrow(x=ge-width+mg, y=y*lane_h+bottom_h, dx=width-2.4*mg, dy=0, width=hw, head_width=hw, head_length=width-2.4*mg, length_includes_head=True, fc=facecolor, lw=0.0)
                    
                    elif strand == -1:
                        ax.arrow(x=gs+width, y=y*lane_h+bottom_h, dx=-1*width, dy=0, width=outer, head_width=outer, head_length=width, length_includes_head=True, fc=edgecolor, lw=0.0)
                        if width > 2.4 * mg:
                            hw = outer * 1.0 * (width-2.4*mg)/width
                            if hw > outer:
                                hw = outer * (width-2.4*mg)/width
                            ax.arrow(x=gs+width-mg, y=y*lane_h+bottom_h, dx=-1*(width-2.4*mg), dy=0, width=hw, head_width=hw, head_length=width-2.4*mg, length_includes_head=True, fc=facecolor, lw=0.0)
                    
                    else:
                        ax.bar([gs], [outer], bottom=y*lane_h+375, width=width, align="edge", fc=facecolor, ec=edgecolor, lw=0.0)

            else:
                if width > 1.2 * head_length:
                    if strand == 1:
                        ax.bar([gs], [outer], bottom=y*lane_h+bottom_h-outer/2, width=width-hl*0.98, align="edge", fc=edgecolor, lw=0.0, zorder=1)
                        ax.bar([gs+mg], [inner], bottom=y*lane_h+bottom_h-inner/2, width=width-mg-hl*0.98, align="edge", fc=facecolor, lw=0.0, zorder=2)
                        ax.arrow(x=ge-hl, y=y*lane_h+bottom_h, dx=hl, dy=0, width=outer, head_width=outer, head_length=hl, length_includes_head=True, fc=edgecolor, lw=0.0, zorder=3)
                        ax.arrow(x=ge-hl, y=y*lane_h+bottom_h, dx=hl-mg*1.4, dy=0, width=inner, head_width=inner, head_length=hl-mg*1.4, length_includes_head=True, fc=facecolor, lw=0.0, zorder=4)
                    
                    elif strand == -1:
                        ax.bar([gs+hl*0.98], [outer], bottom=y*lane_h+bottom_h-outer/2, width=width-hl*0.98, align="edge", fc=edgecolor, lw=0.0, zorder=1)
                        ax.bar([gs+hl*0.98], [inner], bottom=y*lane_h+bottom_h-inner/2, width=width-mg-hl*0.98, align="edge", fc=facecolor, lw=0.0, zorder=2)
                        ax.arrow(x=gs+hl, y=y*lane_h+bottom_h, dx=-1*hl, dy=0, width=outer, head_width=outer, head_length=hl, length_includes_head=True, fc=edgecolor, lw=0.0, zorder=3) 
                        ax.arrow(x=gs+hl, y=y*lane_h+bottom_h, dx=-1*(hl-mg*1.4), dy=0, width=inner, head_width=inner, head_length=hl-mg*1.4, length_includes_head=True, fc=facecolor, lw=0.0, zorder=4)
                    else:
                        ax.bar([gs], [outer], bottom=y*lane_h+bottom_h-outer/2, width=width, align="edge", fc=edgecolor, lw=0.0)
                        ax.bar([gs+mg], [inner], bottom=y*lane_h+bottom_h-inner/2, width=width-2*mg, align="edge", fc=facecolor, lw=0.0)
                else:
                    if strand == 1:
                        ax.arrow(x=ge-width, y=y*lane_h+bottom_h, dx=width, dy=0, width=outer, head_width=outer, head_length=width, length_includes_head=True, fc=edgecolor, lw=0.0)
                        if width > 2.4 * mg:
                            hw = outer * 1.0 * (width-2.4*mg)/width
                            if hw > outer:
                                hw = outer * (width-2.4*mg)/width
                            ax.arrow(x=ge-width+mg, y=y*lane_h+bottom_h, dx=width-2.4*mg, dy=0, width=hw, head_width=hw, head_length=width-2.4*mg, length_includes_head=True, fc=facecolor, lw=0.0)
                    
                    elif strand == -1:
                        ax.arrow(x=gs+width, y=y*lane_h+bottom_h, dx=-1*width, dy=0, width=outer, head_width=outer, head_length=width, length_includes_head=True, fc=edgecolor, lw=0.0)
                        if width > 2.4 * mg:
                            hw = outer * 1.0 * (width-2.4*mg)/width
                            if hw > outer:
                                hw = outer * (width-2.4*mg)/width
                            ax.arrow(x=gs+width-mg, y=y*lane_h+bottom_h, dx=-1*(width-2.4*mg), dy=0, width=hw, head_width=hw, head_length=width-2.4*mg, length_includes_head=True, fc=facecolor, lw=0.0)
                    
                    else:
                        ax.bar([gs], [outer], bottom=y*lane_h+375, width=width, align="edge", fc=facecolor, ec=edgecolor, lw=0.0)

                label_position_list.append((label, width,  gs, ge,  middle,  y, facecolor, edgecolor))
            
            if gs_origin < ge_origin:
                for j in range(gs_origin,ge_origin):
                    gene_position_matrix[y][j] = 1  
            
            else:
                for j in range(gs_origin,length):
                    gene_position_matrix[y][j] = 1 
                for j in range(0, ge_origin):
                    gene_position_matrix[y][j] = 1

            y_list.append(y)
             
    if format != 0:
        if max(y_list)*lane_h+bottom_h < normal_w:
            fig_width = normal_w
        else:
            fig_width = max(y_list)*lane_h+bottom_h

        for tnum, (label, w, gs, ge, x, y, fc, ec) in enumerate(label_position_list): 
            if gs < ge: 
                gmiddle = (ge+gs)/2
            else:
                gmiddle = (ge+gs+2*np.pi)/2  
            
            slide = 0 
            pos_list   = []  
            width_list = []
            for char in label:
                text        = ax2.text(slide, 0, char, ha="right", va="center")
                bbox_text   = text.get_window_extent(renderer=renderer)
                bbox_text   = Bbox(coordinate.transform(bbox_text))
                text.set_visible(False)
                width_list.append(bbox_text.width) 
                pos_list.append(slide+bbox_text.width/2) 
                slide += bbox_text.width
            
            pos_list = [-1*(p-0.5*slide) * 2*fig_width for p in pos_list]
            new_pos_list = [] 
            for pos, width in zip(pos_list, width_list):
                new_pos_list.append((np.arccos(pos/(y*lane_h+bottom_h))-0.5*np.pi+x, y*lane_h+bottom_h, y, pos, x, width))
            
            t_width = (new_pos_list[-1][0] - new_pos_list[0][0])
            if t_width < w-2*head_length:
                if new_pos_list[len(new_pos_list) // 2][0] < 0.5 * np.pi or new_pos_list[len(new_pos_list) // 2][0] > 1.5 * np.pi:
                    rotation = lambda x:(-1.0*x)*180/np.pi 
                else:
                    label    = label[::-1]
                    slide = 0 
                    pos_list   = [] 
                    width_list = []
                    for char in label:
                        text        = ax2.text(slide, 0, char, ha="right", va="center")
                        if label_visible == 0:
                            text.set_visible(False) 
                        bbox_text   = text.get_window_extent(renderer=renderer)
                        bbox_text   = Bbox(coordinate.transform(bbox_text))
                        text.set_visible(False)
                        pos_list.append(slide+bbox_text.width/2) 
                        width_list.append(bbox_text.width) 
                        slide += bbox_text.width
                    
                    pos_list = [-1*(p-0.5*slide) * 2 * fig_width for p in pos_list]
                    new_pos_list = [] 
                    for pos, width in zip(pos_list, width_list):
                        new_pos_list.append((np.arccos(pos/(y*lane_h+bottom_h))-0.5*np.pi+x, y*lane_h+bottom_h, y, pos, x, width))     
                    rotation = lambda x:(-1.0*x)*180/np.pi+180
                
                for char, (theta, height, y, pos, x, width) in zip(label,new_pos_list):
                    ax.text(theta, height, char, ha="center", va="center", rotation=rotation(theta), zorder=10, color=labelcolor)

            else:
                flag = 0 
                sign = 1
                if format == 1:
                    i = 0 
                else:
                    i = new_pos_list[0][2] * 3
                while i < len(text_position_matrix):
                    if format == 1:
                        if i < 3:
                            y = bottom_h - 58 - i * 55
                        else:
                            y = max(y_list)*lane_h+bottom_h + 58 + (i-3) * 55
                    else:
                        if i % 3 == 0:
                            y = (i//3)*lane_h+bottom_h+60
                        elif i % 3 == 1:
                            y = (i//3)*lane_h+bottom_h+113
                        else:
                            y = (i//3)*lane_h+bottom_h+166
                    
                    ts = np.arccos(new_pos_list[0][3]/y)-0.5*np.pi + x
                    te = np.arccos(new_pos_list[-1][3]/y)-0.5*np.pi + x 
                    if ts < 0:
                        ts = 2*np.pi + ts
                    if te > 2*np.pi:
                        te = te-2*np.pi 
                    tts     = int(length * ts / (2 * np.pi))
                    tte     = int(length * te / (2 * np.pi))
                    ts, te = tts, tte

                    sbuf = int(2.0 * length * (np.arcsin(new_pos_list[0][-1]*2*fig_width/y) / (2 * np.pi)))
                    ebuf = int(2.0 * length * (np.arcsin(new_pos_list[-1][-1]*2*fig_width/y) / (2 * np.pi)))
                    tmp_ts = ts - sbuf 
                    if tmp_ts < 0: 
                        tmp_ts = length + tmp_ts 
                    
                    tmp_te = te + ebuf
                    if tmp_te >= length:
                        tmp_te = tmp_te - length

                    if tmp_ts <= tmp_te:
                        middle = (tmp_ts + tmp_te) / 2
                        for j in range(0,int((middle-tmp_ts)*0.6)):
                            if tmp_te+j >= length:
                                break
                            if 1 in text_position_matrix[i][tmp_ts+j:tmp_te+j]:
                                pass 
                            else:
                                sign = 1
                                flag = 1
                                break 
                            
                            
                        if flag == 0:
                            for j in range(int((tmp_te-middle)*0.6),-1,-1):
                                if tmp_ts-j <= 0:
                                    break
                                if 1 in text_position_matrix[i][tmp_ts-j:tmp_te-j]:
                                    pass 
                                else:
                                    sign = -1
                                    flag = 1
                                    break
                                
                    else:
                        j = 0
                        if 1 in text_position_matrix[i][tmp_ts:] or 1 in text_position_matrix[i][0:tmp_te]:
                            pass
                        else:
                            flag = 1
                            break

                    if flag == 0:
                        pass 
                    else: 
                        break
                    i += 1 

                if flag == 0:
                    if format == 1:
                        if i < 3:
                            y = bottom_h - 58 - i * 55
                        else:
                            y = max(y_list)*lane_h+bottom_h + 58 + (i-3) * 55
                    else:
                        y = (i-len(gene_position_matrix)*3) * 60 + len(gene_position_matrix) * lane_h + bottom_h 
                    text_position_matrix.append(np.array([0] * length)) 
                    if tmp_ts <= tmp_te:
                        text_position_matrix[i][tmp_ts:tmp_te] = 1
                    else:
                        text_position_matrix[i][tmp_ts:] = 1
                        text_position_matrix[i][:tmp_te] = 1
                else:
                    if tmp_ts <= tmp_te:
                        text_position_matrix[i][tmp_ts+(j*sign):tmp_te+(j*sign)] = 1
                    else:
                        text_position_matrix[i][tmp_ts:] = 1
                        text_position_matrix[i][:tmp_te] = 1
                
                modified_pos_list = []
                for old_theta, old_height, old_y, pos, x, width in new_pos_list:
                    if flag == 1 and tmp_ts <= tmp_te:
                        new_theta = np.arccos(pos/y)-0.5*np.pi + x + (j*sign)/length * 2 * np.pi
                        target    = np.arccos(pos/y)-0.5*np.pi + x 
                    else:
                        new_theta = np.arccos(pos/y)-0.5*np.pi + x
                        target    = np.arccos(pos/y)-0.5*np.pi + x 
                    modified_pos_list.append([new_theta, target, old_theta, y, old_height, width]) 

                    
                if modified_pos_list[len(new_pos_list) // 2][0] < 0.5 * np.pi or modified_pos_list[len(new_pos_list) // 2][0] > 1.5 * np.pi:
                    direction = 1
                    rotation = lambda x:(-1.0*x)*180/np.pi 
                else:
                    direction = -1
                    label    = label[::-1]
                    slide = 0 
                    pos_list = [] 
                    for char in label:
                        text        = ax2.text(slide, 0, char, ha="right", va="center")
                        bbox_text   = text.get_window_extent(renderer=renderer)
                        bbox_text   = Bbox(coordinate.transform(bbox_text))
                        text.set_visible(False)
                        pos_list.append(slide+bbox_text.width/2) 
                        slide += bbox_text.width
                    pos_list = [-1*(p-0.5*slide) * 2*fig_width for p in pos_list]
                    
                    for p, pos in enumerate(pos_list):
                        if flag == 1:
                            new_theta = np.arccos(pos/y)-0.5*np.pi + x + (j*sign)/length * 2 * np.pi
                            target    = np.arccos(pos/y)-0.5*np.pi + x 
                        else:
                            new_theta = np.arccos(pos/y)-0.5*np.pi + x
                            target    = np.arccos(pos/y)-0.5*np.pi + x 
                        
                        modified_pos_list[p][0] = new_theta
                        modified_pos_list[p][1] = target
                    
                    rotation = lambda x:(-1.0*x)*180/np.pi+180
                
                for char, (theta, target, old_theta, y, old_height, width) in zip(label, modified_pos_list):
                    ax.text(theta, y, char, ha="center", va="center", rotation=rotation(theta), color=labelcolor)
                
                mid = modified_pos_list[len(new_pos_list) // 2]
                s = modified_pos_list[0]
                e = modified_pos_list[-1] 
                sbuf = sbuf*0.4*2*np.pi/length
                ebuf = ebuf*0.4*2*np.pi/length
                sbuf, ebuf = max([sbuf, ebuf]), max([sbuf, ebuf]) 
                ax.plot([gmiddle, gmiddle], [mid[3], mid[4]], lw=0.5, color="k", zorder=0)
                
                if direction == 1:
                    ax.bar([s[0]-sbuf], [40], width=e[0]-s[0]+sbuf+ebuf, bottom=s[3]-20, lw=0.5, ec=ec, fc=fc, align="edge")
                elif direction == -1:
                    ax.bar([s[0]-ebuf], [40], width=e[0]-s[0]+sbuf+ebuf, bottom=s[3]-20, lw=0.5, ec=ec, fc=fc, align="edge")
               
                ty_list.append(y) 

    y_set = list(set(y_list)) 
    y_set.sort() 
    
    if format == 0:
        ylim = max(y_list)*lane_h + bottom_h + 2.5 * lane_h
    else:
        ylim = max([max(y_list)*lane_h + bottom_h, max(ty_list)]) + 3 * lane_h
        
    if axis_visible == True:
        if tick_space == "auto":
            if   length < 1000:
                space = 250
            elif length < 2000:
                space = 500 
            elif length < 4000:
                space = 800
            elif length < 6000:
                space = 1000
            elif length < 10000:
                space = 2000
            else:
                space = 3000
        else:
            space = tick_space
        
        ax.bar([gs], [8], bottom=ylim-lane_h*0.80-4, width=2*np.pi, align="edge", fc="#606060", ec=edgecolor, lw=0.0, zorder=0)
        if space > 0:
            for pos in range(0, length, space):
                ax.plot([2*np.pi*(pos/length), 2*np.pi*(pos/length)], [(ylim-lane_h*0.80-4), (ylim-lane_h*0.80-4)-0.025*(normal_w-lane_h)], color="k", lw=1)
                if pos == 0:
                    ax.text(2*np.pi*(pos/length), ylim-lane_h*0.80-0.075*(normal_w-lane_h), str(1), ha="center", va="center",  rotation=-1*180/np.pi*(2*np.pi*(pos/length)), fontsize=fontsize)
                else:
                    if 2*np.pi*(pos/length) < 0.5 * np.pi or  2*np.pi*(pos/length) > 1.5 * np.pi: 
                        ax.text(2*np.pi*(pos/length), ylim-lane_h*0.80-0.075*(normal_w-lane_h), str(pos), ha="center", va="center",  rotation=-1*180/np.pi*(2*np.pi*(pos/length)), fontsize=fontsize)
                    else:
                        ax.text(2*np.pi*(pos/length), ylim-lane_h*0.80-0.06*(normal_w-lane_h), str(pos), ha="center", va="center",  rotation=(-1*180/np.pi*(2*np.pi*(pos/length)))+180, fontsize=fontsize)

    if ylim < normal_w:
        ylim = normal_w
        ax.set_ylim(0, normal_w)
    else:
        ax.set_ylim(0, ylim)
   
    ax.patch.set_alpha(0.0) 
    ax2.patch.set_alpha(0.0)
    
    if ylim > normal_w:
        fig.set_size_inches(6 * ylim/fig_width, 6 * ylim/fig_width)
    return ax, y_list, ty_list, fig_width, ylim, bottom_h

def visualize(brick, format=0, feature_list=None, unvisible_types=["source"], visible_types=[], bottom=None, fontsize=8, label_visible=True, axis_visible=True, title_visible=True, tick_space="auto", labelcolor="k", titlename=None, fig=None):
    if titlename is None:
        titlename = brick.project
    brick = copy.deepcopy(brick) 
    if fig is None:
        figure  = plt.figure(figsize=(6,6))
        basenum = 0
    else:
        figure  = fig 
        basenum = len(fig.axes) + 1
    ax      = figure.add_axes([0,0,1,1], polar=True, label="hoge"+str(basenum))
    ax2     = figure.add_axes([0,0,1,1], label="fuga"+str(basenum))
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.spines['polar'].set_visible(False)
    ax.xaxis.set_ticks([])
    ax.xaxis.set_ticklabels([])
    ax.yaxis.set_ticks([])
    ax.yaxis.set_ticklabels([])
    ax2.set_axis_off()
    if feature_list is None:
        feature_list = birck.dnafeatures
        feature_list.sort(key=lambda x:len(brick.getdnaseq(x.start, x.end)))
    else:
        pass
    ax, y_list, ty_list, fig_width, ylim, bottom = map_feat(figure, ax, ax2, feature_list, len(brick.seq), unvisible_types=unvisible_types, visible_types=visible_types, format=format, bottom=bottom, enlarge=1.0, label_visible=label_visible, fontsize=fontsize, axis_visible=axis_visible, tick_space=tick_space, labelcolor=labelcolor) 
   
    if title_visible == True:
        renderer    = figure.canvas.get_renderer()
        coordinate  = ax2.transData.inverted() 
        text        = ax2.text(0.5, 0.5, titlename, ha="center", va="center", fontsize=fontsize*1.125 if fontsize >= 8 else 10)
        bbox_text   = text.get_window_extent(renderer=renderer)
        bbox_text   = Bbox(coordinate.transform(bbox_text))
        bp_text     = ax2.text(0.5, 0.5-bbox_text.height-0.01, str(len(brick.seq)) + " bp", ha="center", va="center", fontsize=fontsize if fontsize >= 8 else 8)
        if bbox_text.width/2 > bottom:
            text.set_visible(False)
            bp_text.set_visible(False) 

    figure.patch.set_alpha(0.0)  
    return figure, ax

if __name__ == "__main__":
    from dna import *
    brick = Dbrick(record=sys.argv[1])
    brick.name = sys.argv[1].split("/")[-1].replace(".gbk","")  
    fig   = visualize(brick, format=2, unvisible_types=["primer_bind"], bottom=400) 
    fig.patch.set_alpha(0.0) 
    fig.savefig("{}.pdf".format(brick.name), bbox_inches="tight")
