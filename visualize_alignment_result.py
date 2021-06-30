import os 
import sys 
import visualize_linear_dna as vl 

def read_fasta(fasta_name):
    """Read fasa file
    """
    seq_dict = {}
    with open(fasta_name) as f:
        for line in f:
            if line[0] == ">":
                key = line.rstrip()[1:]
                seq_dict[key] = []
            else:
                seq_dict[key].append(line.rstrip()) 
    
    for key in seq_dict:
        seq_dict[key] = "".join(seq_dict[key])

    return seq_dict

def colorbar_align(ax, color_dict, template, sample, char=True, fontsize=10):
    bars = ax.bar(list(range(len(template))), [0.9] * (len(template)), width=1.0, edgecolor="#BBBBBB", linewidth=0.0, align="edge",bottom=0.05)
    ax.set_xlim(0,len(template))
    ax.set_ylim(0,1.00)
    
    p = 0
    for bar, t, s in zip(bars, template, sample):
        if t != s:
            color = color_dict[s]
        else:
            color = color_dict[t]
        bar.set_facecolor("#FFFFFF")
        if  t != s:
            if char == True:
                bar.set_alpha(0.7)
                ax.text(p+0.5, 0.45, s, va="center", ha="center", fontsize=fontsize, fontstyle="bold", zorder=100)
            bar.set_linewidth(0.4)
            bar.set_edgecolor("#FF0000")
        else:
            if char == True:
                ax.text(p+0.5, 0.45, s, va="center", ha="center", fontsize=fontsize, zorder=100)
            bar.set_linewidth(0.2)
            bar.set_edgecolor("#BBBBBB") 
        p += 1 
    
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.patch.set_alpha(0.0)
    return bars

def visualize_alignment(brick, msafasta, start=0, end=None, wrap_width=500, annotation_loc=None, label_box=True, feature_list=None, unvisible_types=["source"], visible_types=[], enlarge_w=1.0, enlarge_h=1.0, scale="fix", fontsize=12, nucl_char=None, nucl_color_dict=None, view_title=True, view_axis=True, tick_space="auto"):
    def generate_newbrick(brick, new_seq):
        pos            = 0 
        prechar        = "-"
        gaplength      = 0 
        positions      = [0]
        gaplength_dict = []
        for char in new_seq:
            if char != "-":
                if prechar == "-":
                    gaplength_dict[pos] = gaplength
                pos += 1
                gaplength = 0
            else:
                gaplength += 1
                if prechar != "-":
                    positions.append(pos)
            prechar = char
        
        if prechar == "-":
            gaplength_dict[pos] = gaplength
        else:
            positions.append(pos) 
    
        poscombi_list = [] 
        for i in range(len(positions)-1):
            poscombi_list.append((positions[i], positions[i+1])) 

        if gaplength_dict[poscombi[0][0]] == 0:
            new_brick = brick[poscombi[0][0]:poscombi[0][1]]
        else:
            new_brick = DNA(seq="-" * gaplength_dict[poscombi[0][0]]) 
            new_brick = new_brick + brick[poscombi[0][0]:poscombi[0][1]]
        
        for poscombi in poscombi_list[1]:
            if poscombi[0] in gaplength_dict:
                new_brick = new_brick + DNA(seq="-" * gaplength_dict[poscombi[0]]) 
            new_brick = new_brick + brick[poscombi[0]:poscombi[1]] 
        
        if positions[-1] in gaplength_dict:
            new_brick = new_brick + DNA(seq="-" * gaplength_dict[poscombi[0]]) 
    
    seq_dict = read_fasta(msafasta) 
    seq_keys = list(seq_dict.keys())
    template = seq_dict[seq_keys[0]]
    
    brick = generate_newbrick(brick, template)
    brick = copy.deepcopy(brick) 
    width = wrap_width 
    if nucl_color_dict == None:
        nucl_color_dict = color_dict  

    if annotation_loc is None:
        annotation_loc = "top"

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
    
    if feature_list is None:
        feature_list = brick.dnafeatures    
    else:
        brick.dnafeatures = copy.deepcopy(feature_list)
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
            
            if "label" in feat.qualifiers:
                if type(feat.qualifiers["label"]) == list:
                    label = feat.qualifiers["label"][0]
                else:
                    label = feat.qualifiers["label"]
            else:
                label = feat.type

            if "broken_feature" in feat.qualifiers:
                posinfo     = feat.qualifiers["broken_feature"][0].split("]")[-1]
                feat_length = feat.qualifiers["broken_feature"][0].split("]")[0].split(":")[-3] 
                label = label + posinfo + ":" +  feat_length
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
        if scale == "fix":
            std = 25
            head_length = 2
        else:
            if scale == "auto":
                scale = width
            else:
                pass 
            if scale < 101:
                std = 25
                head_length = 2
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

        zero_position = sub_start + 1
        ax  = fig.add_axes([0, 0, enlarge_w*len(sub_brick.seq)/std, 1.0*enlarge_h], label=str(num)) 
        ax, y_list, ty_list = vl.map_feat(fig, ax, sub_brick.dnafeatures, len(sub_brick.seq), head_length, unvisible_types=unvisible_types, visible_types=visible_types, enlarge_w=enlarge_w, enlarge_h=enlarge_h, annotation_loc=annotation_loc, label_box=label_box, fontsize=fontsize, project=sub_brick.project, view_title=view_title, view_axis=view_axis)
        
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
        ax.set_xticks([]) 
        ax.set_xticklabels([])
        ax.set_xlim(0, len(sub_brick.seq))
        ax.spines["bottom"].set_visible(False)  
        
        ax_seq = fig.add_axes([0, -0.65*enlarge_h/(abs(ytop-ybottom)), enlarge_w*len(sub_brick.seq)/std, 0.6*enlarge_h/(abs(ytop-ybottom))], label="{}_seq".format(num))
        if nucl_char != True and nucl_char != False:
            if width > 500:
                vl.colorbar(ax_seq, nucl_color_dict, sub_brick.seq, char=False, fontsize=fontsize-2)
            else:
                vl.colorbar(ax_seq, nucl_color_dict, sub_brick.seq, char=True, fontsize=fontsize-2)
        else:
            vl.colorbar(ax_seq, nucl_color_dict, sub_brick.seq, char=nucl_char)

        ticks = np.array(ticks) + 0.5
        #ax_seq.set_xticks(ticks) 
        #ax_seq.set_xticklabels(tick_labels)
        #ax_seq.set_xlim(0, len(sub_brick.seq)) 
        if num == 0:
            if view_title == True:
                ax.set_title(brick.project, fontsize=fontsize * 1.125)
            ax.set_position([0, 0, enlarge_w*len(sub_brick.seq)/std, 1.0*enlarge_h*(abs(ytop-ybottom))]) 
            ax_seq.set_position([0, -0.65*enlarge_h, enlarge_w*len(sub_brick.seq)/std, 0.6*enlarge_h]) 
            ceil = -0.65*enlarge_h
        else:
            ax.set_position([0, ceil-1.2*enlarge_h-1.0*enlarge_h*(abs(ytop-ybottom)), enlarge_w*len(sub_brick.seq)/std, 1.0*enlarge_h*(abs(ytop-ybottom))]) 
            ax_seq.set_position([0, ceil-1.2*enlarge_h-1.0*enlarge_h*(abs(ytop-ybottom))-0.65*enlarge_h, enlarge_w*len(sub_brick.seq)/std, 0.6*enlarge_h]) 
            ceil = ceil-1.2*enlarge_h-1.0*enlarge_h*(abs(ytop-ybottom))-0.65*enlarge_h
                        
        for samplenum, key in enumerate(seq_keys[1:]):
            sample_seq = seq_dict[key][sub_start:sub_end]
            ax_seq = fig.add_axes([0, -0.65*enlarge_h/(abs(ytop-ybottom)), enlarge_w*len(sub_brick.seq)/std, 0.6*enlarge_h/(abs(ytop-ybottom))], label="{}_{}_seq".format(num, key))
            colorbar_align(ax_seq, nucl_color_dict, sub_brick.seq, sample_seq, char=nucl_char)
            ax_seq.set_position([0, ceil-1.2*enlarge_h-1.0*enlarge_h*(abs(ytop-ybottom))-0.65*enlarge_h-0.62*samplenum, enlarge_w*len(sub_brick.seq)/std, 0.6*enlarge_h]) 
            ax_seq.set_xticks([]) 
            ax_seq.set_xlim(0, len(sub_brick.seq))

        ceil = ceil-1.2*enlarge_h-1.0*enlarge_h*(abs(ytop-ybottom))-0.65*enlarge_h-0.62*samplenum
        ax_seq.set_xticks(ticks) 
        ax_seq.set_xticklabels(tick_labels)
        ax_seq.set_xlim(0, len(sub_brick.seq)) 
        if view_axis == False:
            ax_seq.spines["bottom"].set_visible(False) 
            ax_seq.set_xticks([]) 
        axes_list.append((ax, ax_seq)) 
    #fig.set_size_inches(3, 0.35*(abs(ytop-ybottom)))
    return fig, axes_list  
