"""
Functions for visualizing identified TSMs.
Visualizes a phylogenetic tree and the related sequence alignment,
and indicates the child and parent sequences and structures
between which the TSM has occurred. In addition, visualizes
secondary structures and sequence for the given list of internal
and external tree nodes. Requires a newick-formatted tree file,
in which the intenal nodes are named as #name. In addition,
a sequence alignment file is required, which headers match
to node names in a tree file.

Usage:

Hit class
==========

    The following data classes has to be defined:


    hit(index_alignment_visualization_started,
        index_alignment_visualization_end,
        reference_node,
        query_node,
        tsm_target_start_index,
        tsm_target_end_index,
        tsm_source_start_index,
        tsm_source_end_index,
        aligned_struct_ref,
        aligned_struct_qry,
        tsm_parent_struct,
        tsm_child_struct,
        structural_change_as_string,
        sequence_to_show,
        structures_to_show)

    
    files(align_dict,
          outputfile,
          struct_files)

    make_fig(files_class,
             features_dict,
             hit_class,
             par,
             tree)


    Authors: Ari Löytynoja
             Heli Mönttinen
        
"""

from Bio import SeqIO
import re
from PIL import Image, ImageDraw, ImageFont, ImageColor
from dataclasses import dataclass
from typing import List, Dict

# parameters defining the size and spacing of everything
@dataclass
class params:
    vnuc: int = 15        # vertical size
    hnuc: int = 12        # horizontal size
    voff: int = 1         # offset for character
    hoff: int = int(hnuc/3)
    tmar: int = 80        # top margin for headers etc
    vmar: int = 2*vnuc    # vertical margin
    hmar: int = 2*hnuc    # horizontal margin
    lmar: int = 10*hnuc   # left margin
    vspc: int = 1         # vertical spacing
    hspc: int = 0         # horizontal spacing

    htre: int = 100       # tree width (pixels)
    hstp: int = 0         # horizontal step (in tree) = htre/count_levels(tree)
    vstp: int = vnuc+vspc # vertical step (in tree)
    hvstp: int = vstp/2   # half vertical step (in tree)

    alpha: int = 80       # nuc color shading
    noder: int = 3        # dot radius
    nodeshft: int = -2    # dot horizontal shift
    ancline: bool = True  # draw lines for shown anc nodes

# definition of a fpa hit
@dataclass
class hit:
    align_start: int
    align_end: int
    ancestor_node: str
    query_node: str
    qry_ts_start: int
    qry_ts_end: int
    ref_ts_start: int
    ref_ts_end: int
    structure_ref: str
    structure_qry: str
    ts_region_ref: str
    ts_region_qry: str
    structural_change: str
    showseq: List[str]
    showstruct: List[str]
    compensating: str
    new_loop: str
    inverted_loop: str

@dataclass
class files:
    fasta_dict: Dict[str,str]
    outputfile: str
    struct_files: Dict[str,str]
    
# count (nested) levels  
def count_levels(node):
    if node.is_leaf():
        return 0
    else:
        return max(count_levels(node.children[0]),count_levels(node.children[1]))+1

# count nodes in subtree with names on given list
def count_nodes_on_list(node,keeplist,structlist):
    if node.is_leaf():
        s = 0
        if any(node.name in n for n in keeplist):
            s = 1
            if any(node.name in n for n in structlist):
                s = 2
        return s
    else:
        s = 0
        if any(node.name in n for n in keeplist):
            s = 1
            if any(node.name in n for n in structlist):
                s = 2
        return count_nodes_on_list(node.children[0],keeplist,structlist) + count_nodes_on_list(node.children[1],keeplist,structlist) + s

def dashedhline(x1,y1,x2,len1,len2,f,w, dr):
    i=x1
    while(i<x2):
        dr.line((i,y1,i+len1,y1),fill=f,width=w)
        i+=len1+len2

# main drawing function: from a point, draw a horizontal branch and vertical lines
def draw_next_node(node,taxashown,x1,y1,structnodes, nodeposxy, par, dr):
    if par.ancline and not node.is_leaf() and any(node.name in n for n in taxashown):
        dashedhline(x1,y1,par.htre+par.hmar,2,2,"gray",1, dr)
        
    x2 = par.htre+par.hmar-count_levels(node)*par.hstp
    dr.line((x1,y1,x2,y1),fill="black",width=1)
    nodeposxy[node.name] = (x2,y1)

    vshf = 0
    vshf0 = 0
    vshf1 = 0
    # this anc node shown -> half vstep to both directions
    if any(node.name in n for n in taxashown):
        vshf+=par.vstp/2
        # its structure is also shown -> vstep more down
        if any(node.name in n for n in structnodes):
            vshf1+=par.vstp 
    
    # child0 is anc node and shown -> full vstep up
    if not node.children[0].is_leaf() and any(node.children[0].name in n for n in taxashown):
        vshf0+=par.vstp/2
        # its structure is also shown
        if any(node.children[0].name in n for n in structnodes):
            vshf0+=par.vstp
        
    # child1 is anc node and shown -> full vstep down
    if not node.children[1].is_leaf() and any(node.children[1].name in n for n in taxashown):
        vshf1+=par.vstp/2

    vshf0s = 0
    # child0 is leaf and has structure
    if node.children[0].is_leaf() and any(node.children[0].name in n for n in structnodes):
        vshf0s = par.vstp/2 

    vshf1s = 0
    # child1 is leaf and has structure
    if node.children[1].is_leaf() and any(node.children[1].name in n for n in structnodes):
        vshf1s = par.vstp/2 
    
    nc0 = 1
    nc00 = 1
    nc01 = 1
    nc0 = count_nodes_on_list(node.children[0],taxashown,structnodes)
    if(not node.children[0].is_leaf()):
        nc00 = count_nodes_on_list(node.children[0].children[0],taxashown,structnodes)
        nc01 = count_nodes_on_list(node.children[0].children[1],taxashown,structnodes)
    hbl0 = nc0*par.vstp
    
    nc1 = 1
    nc10 = 1
    nc11 = 1
    nc1 = count_nodes_on_list(node.children[1],taxashown,structnodes)
    if(not node.children[1].is_leaf()):
        nc10 = count_nodes_on_list(node.children[1].children[0],taxashown,structnodes)
        nc11 = count_nodes_on_list(node.children[1].children[1],taxashown,structnodes)
    hbl1 = nc1*par.vstp
        
    x1 = x2
    y2 = y1-nc01*par.vstp-vshf0-vshf-vshf0s
    if(nc0==1):
        y2+=par.hvstp
        
    dr.line((x1,y1,x1,y2),fill="black",width=1)

    if(not node.children[0].is_leaf()):
        draw_next_node(node.children[0],taxashown,x2,y2,structnodes,nodeposxy,  par, dr)
    else:
        dr.line((x2,y2,par.htre+par.hmar,y2),fill="black",width=1)
        nodeposxy[node.children[0].name] = (x2,y2)
        
    y2 = y1+nc10*par.vstp+vshf1+vshf-vshf1s
    if(nc1==1):
        y2-=par.hvstp

    dr.line((x1,y1,x1,y2),fill="black",width=1)

    if(not node.children[1].is_leaf()):
        draw_next_node(node.children[1],taxashown,x2,y2,structnodes, nodeposxy, par, dr)
    else:
        dr.line((x2,y2,par.htre+par.hmar,y2),fill="black",width=1)
        nodeposxy[node.children[1].name] = (x2,y2)

# draw a tree
def draw_tree(tree,taxashown,structnodes, nodeposxy, par, dr):
    nc = count_nodes_on_list(tree,taxashown,structnodes)
    nc0 = count_nodes_on_list(tree.children[0],taxashown,structnodes)
    nc1 = count_nodes_on_list(tree.children[1],taxashown,structnodes)

    hbl0 = nc0*par.vstp

    x1 = par.hmar
    y1 = par.tmar+par.vmar+hbl0

    if any(tree.name in n for n in taxashown):
        y1+=par.hvstp 

    dr.line((x1-3,y1,x1,y1),fill="black",width=1)

    draw_next_node(tree,taxashown,x1,y1,structnodes, nodeposxy, par, dr)

# split a string to chars
def split(word): 
    return [char for char in word]  

# give the color for a nucleotide
def nuccol(char): # https://stackoverflow.com/questions/54165439
    if(char == 'A'):
        return 'khaki'
    if(char == 'C'):
        return 'lightblue'
    if(char == 'G'):
        return 'lightcoral'
    if(char == 'T'):
        return 'teal'
    else: return 'white'

# draw sequences marked to be shown
def draw_sequences(spos,epos,names,isshown,structnodes, par, dr, fn, seqs, namesshown, species_names):
    k=0
    for i in range(len(names)):
        if(isshown[i]):
            dr.text((par.htre+2.5*par.hnuc,par.tmar+par.vmar+k*(par.vnuc+par.vspc)+par.voff), species_names[i][0:100], font=fn, fill=(20,20,20))
            nucs = split(seqs[i])
            l=0
            for j in range(spos,epos):
                cl = nuccol(nucs[j])
                cl = list(ImageColor.getrgb(cl))
                cl.append(par.alpha)
                cl = tuple(cl)
                x1 = par.htre+par.lmar+l*(par.hnuc+par.hspc)
                y1 = par.tmar+par.vmar+k*(par.vnuc+par.vspc)
                x2 = par.htre+par.lmar+l*(par.hnuc+par.hspc)-par.hspc+par.hnuc
                y2 = par.tmar+par.vmar+k*(par.vnuc+par.vspc)-par.vspc+par.vnuc
                dr.rectangle([x1,y1,x2,y2], fill=cl)
                dr.text((x1+par.hoff,y1+par.voff), nucs[j], font=fn, fill=(20,20,20))
                l+=1
            namesshown.append(names[i])
            k+=1
            if any(names[i] in n for n in structnodes):
                k+=1
                namesshown.append(names[i]+"struct")

    return namesshown

# draw boxes around specific regions and dots for those nodes
def draw_boxes_dots(hit,namesshown, seqs, names, nodeposxy, spos, par, fn, dr):
    draw_box(hit.query_node,hit.qry_ts_start,hit.qry_ts_end,namesshown, seqs, names, spos, par, fn, dr)
    draw_box(hit.ancestor_node,hit.ref_ts_start,hit.ref_ts_end,namesshown, seqs, names, spos, par, fn, dr)
    draw_dot(hit.query_node,nodeposxy,'orange', par, dr)
    draw_dot(hit.ancestor_node,nodeposxy,'red', par, dr)

# draw (multiple) boxes around specific regions
def draw_box(name,start,end,namesshown, seqs, names, spos, par, fn, dr):
    sp = start
    ep = end
    start-=spos
    end-=(spos+1)  # ETE cannot read hashes
    k = namesshown.index(name)
    nucs = split(seqs[names.index(name)])
    for j in range(sp,ep):
        cl = nuccol(nucs[j])
        x1 = par.htre+par.lmar+(j-spos)*(par.hnuc+par.hspc)
        y1 = par.tmar+par.vmar+k*(par.vnuc+par.vspc)
        x2 = par.htre+par.lmar+(j-spos)*(par.hnuc+par.hspc)-par.hspc+par.hnuc
        y2 = par.tmar+par.vmar+k*(par.vnuc+par.vspc)-par.vspc+par.vnuc
        dr.rectangle([x1,y1,x2,y2], fill=cl)
        dr.text((x1+par.hoff,y1+par.voff), nucs[j], font=fn, fill=(20,20,20))
    x1 = par.htre+par.lmar+start*(par.hnuc+par.hspc)-par.hspc
    y1 = par.tmar+par.vmar+k*(par.vnuc+par.vspc)-par.vspc
    x2 = par.htre+par.lmar+end*(par.hnuc+par.hspc)+par.hnuc
    y2 = par.tmar+par.vmar+k*(par.vnuc+par.vspc)+par.vnuc
    dr.rectangle([x1-1,y1-1,x2+1,y2+1], outline="white")
    dr.rectangle([x1+1,y1+1,x2-1,y2-1], outline="white")
    dr.rectangle([x1,y1,x2,y2], outline="navy")

# draw (multiple) boxes around specific regions
def draw_box_only(name,start,end,namesshown, color, spos, par, fn, dr):
    sp = start
    ep = end
    start-=spos
    end-=(spos+1)  # ETE cannot read hashes
    k = namesshown.index(name)

    x1 = par.htre+par.lmar+start*(par.hnuc+par.hspc)-par.hspc
    y1 = par.tmar+par.vmar+k*(par.vnuc+par.vspc)-par.vspc
    x2 = par.htre+par.lmar+end*(par.hnuc+par.hspc)+par.hnuc
    y2 = par.tmar+par.vmar+k*(par.vnuc+par.vspc)+par.vnuc
    dr.rectangle([x1-1,y1-1,x2+1,y2+1], outline="white")
    dr.rectangle([x1+1,y1+1,x2-1,y2-1], outline="white")
    dr.rectangle([x1,y1,x2,y2], outline=color)

# draw structure annotation for specific nde
def draw_struct(node,struct,spos,epos,namesshown, par, fn, dr):
    k = namesshown.index(node+"struct")
    nucs = split(struct)
    for j in range(spos,epos):
        x1 = par.htre+par.lmar+(j-spos)*(par.hnuc+par.hspc)
        y1 = par.tmar+par.vmar+k*(par.vnuc+par.vspc)
        dr.text((x1+par.hoff,y1+par.voff), nucs[j], font=fn, fill=(20,20,20))
        
def draw_dot(node_to_draw,nodeposxy, color, par, dr):
    # ETE cannot read hashes
    x,y = nodeposxy[node_to_draw]
    dr.ellipse((x-par.noder+par.nodeshft,y-par.noder,x+par.noder+par.nodeshft,y+par.noder),fill=color,outline="black")


# draw bottom ruler
def draw_ruler(spos, epos, par, dr, snum):
    fn = ImageFont.truetype('VeraMono',9) #/usr/share/fonts/truetype/ubuntu/UbuntuMono-R.ttf
    y1 = par.tmar+par.vmar+snum*(par.vnuc+par.vspc)
    for i in range(spos,epos):
        p = i-spos
        if(i%10==0):
            x1 = par.htre+par.lmar+p*(par.hnuc+par.hspc)
            dr.text((x1+par.hoff,y1+par.voff), str(i), font=fn, fill=(20,20,20))

# add headers (or other text)
def add_header(text,x,y,fontsize,color, dr):
    fn = ImageFont.truetype('VeraMono',fontsize) #/usr/share/fonts/truetype/ubuntu/UbuntuMono-R.ttf
    dr.text((x,y), text, font=fn, fill=color)


def draw_structures(data,names,seqs,spos,epos,namesshown, par, fn, dr):
    files = data.struct_files
    for node in files:
        with open(files[node], "r") as file:
            for struct in file:
                pass
        struct = struct.split(" ")[0]
        seq = split(seqs[names.index(node)])
        i=0
        struct_gaps = ""

        for c in seq:
            if c!="-":
                struct_gaps+=struct[i:(i+1)]
                i+=1

            else:
                struct_gaps+=" "
                
        draw_struct(node,struct_gaps,spos,epos,namesshown, par, fn, dr)
        
    
def make_fig(data,
             hit,
             par,
             tree):

    # ETE cannot read hashes

    # horizontal step in tree can only be known after seeing the tree
    par.hstp = par.htre/count_levels(tree)

    names = []         # all node/leaf names
    seqs = []          # all sequences
    isshown = []       # boolean: node/leaf shown or not
    taxashown = []     # names of shown nodes/leafs
    structnodes = []   # names of nodes/leafs with structure
    species_names = [] # species names
    nodeposxy = {}     # positions of nodes
    
    for fasta in data.fasta_dict:
        name, sequence = fasta, data.fasta_dict[fasta]
        if re.match("^#",name):
            if name==hit.ancestor_node or name==hit.query_node or name in hit.showstruct:
                #name = re.sub("#","A",name)
                isshown.append(True)
                taxashown.append(name)
                structnodes.append(name)
            elif any (name in n for n in hit.showseq):
                isshown.append(True)
                taxashown.append(name)
                if any (name in n for n in hit.showstruct):
                    structnodes.append(name)
            else:
                isshown.append(False)  # ETE cannot read hashes: names have to match
            species_names.append(name)
        else:
            isshown.append(True)
            taxashown.append(name)
            species_names.append(name)
            if name in hit.showstruct:
                structnodes.append(name)
        names.append(name)
        seqs.append(sequence)

    snum = sum(isshown)+len(structnodes)

    ###
 
    # start and end positions in alignment
    #
    spos = hit.align_start  # start position
    epos = hit.align_end    # end position
    slen = epos-spos+1


    aliw = (par.hnuc+par.hspc)*slen-par.hspc+par.lmar+par.hmar   # alignment width
    trew = par.htre+2*par.hmar                                   # tree width
    alih = par.tmar+(par.vnuc+par.vspc)*snum-par.vspc+2*par.vmar # alignment height

    ###

    canvas = (trew+aliw, alih)

    im = Image.new('RGBA', canvas, (255, 255, 255))
    dr = ImageDraw.Draw(im)

    draw_tree(tree, taxashown, structnodes, nodeposxy, par,dr)

    fn = ImageFont.truetype('VeraMono',12) #/usr/share/fonts/truetype/ubuntu/UbuntuMono-R.ttf

    namesshown = []
    draw_sequences(spos, epos, names, isshown, structnodes, par, dr, fn, seqs, namesshown, species_names)
 
    # note: positions as alignment columns!
    #
    draw_boxes_dots(hit, namesshown, seqs, names, nodeposxy, spos, par, fn, dr)
    draw_box_only(hit.ancestor_node,hit.qry_ts_start,hit.qry_ts_end,namesshown, "darkred", spos, par, fn, dr)
    draw_structures(data,names,seqs,spos,epos,namesshown, par, fn, dr)

    draw_ruler(spos,epos, par, dr, snum)

    ts_structures = 'Structure in template switch region:  ref - ' + hit.ts_region_ref + ' qry - ' + hit.ts_region_qry

    #
    add_header('Compensating mutations: ' +  hit.compensating,20,70,15,"black", dr)

    if 'pdf' not in data.outputfile:
        im.save(data.outputfile)

    else:

        rgb = Image.new('RGB', im.size, (255, 255, 255))  # white background
        rgb.paste(im, mask=im.split()[3])
        rgb.save(data.outputfile, 'PDF', resolution=300.0)

