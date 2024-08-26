# -*- coding: utf-8 -*-

#Directory
dir=''

#Project
Project='Projects/granar_example/' # 

inputs='in/'
Gen='General.xml'
Geom='Geometry.xml'
Hydr='Hydraulics.xml'
BC='BC.xml'
Horm='Hormones_Carriers.xml'
Cell_connec_max=50
Ncellperimeters=100
V_modifier=1.0
test_mass_balance=0
maxCell2ThickWalls=101
matrix_analysis=0


#Libraries
print('Importing libraries')
import time
t0 = time.perf_counter()
    
import numpy as np #NumPy is the fundamental package for scientific computing with Python.
                   #It contains among other things:
                   #- a powerful N-dimensional array object
                   #- sophisticated (broadcasting) functions
                   #- tools for integrating C/C++ and Fortran code
                   #- useful linear algebra, Fourier transform, and random number capabilities
#import sympy as smp
#import mpmath as mp
#mp.mp.dps = 20 #Floating point precision
from numpy import genfromtxt #Load data from a text file, with missing values handled as specified.
from numpy.random import *  # for random sampling

from scipy import sparse
from scipy.sparse import lil_matrix,csr_matrix,spdiags,identity
from scipy.sparse.linalg import spsolve
import scipy.linalg as slin #Linear algebra functions

import math

import pylab #Found in the package pyqt, also needs to be installed for proper use of matlab-python connectivity
from pylab import *  # for plotting

from decimal import Decimal

import networkx as nx #NetworkX is a Python language software package for the creation,
                      #manipulation, and study of the structure, dynamics, and functions of complex networks.
from lxml import etree #Tree element analysis module

import sys, os # On importe le module os qui dispose de variables 
               # et de fonctions utiles pour dialoguer avec votre 
               # système d'exploitation
t1 = time.perf_counter()
print(t1-t0, "seconds process time")

#Precision
dp = np.dtype((np.float64))

#Import General data
print('Importing geometrical data')
t0 = time.perf_counter()
OS=etree.parse(dir + Project + inputs + Gen).getroot().xpath('OS')[0].get("value")
Output_path=etree.parse(dir + Project + inputs + Gen).getroot().xpath('Output')[0].get("path")
Paraview=int(etree.parse(dir + Project + inputs + Gen).getroot().xpath('Paraview')[0].get("value"))
ParaviewWF=int(etree.parse(dir + Project + inputs + Gen).getroot().xpath('Paraview')[0].get("WallFlux"))
ParaviewMF=int(etree.parse(dir + Project + inputs + Gen).getroot().xpath('Paraview')[0].get("MembraneFlux"))
ParaviewPF=int(etree.parse(dir + Project + inputs + Gen).getroot().xpath('Paraview')[0].get("PlasmodesmataFlux"))
ParaviewWP=int(etree.parse(dir + Project + inputs + Gen).getroot().xpath('Paraview')[0].get("WallPot"))
ParaviewCP=int(etree.parse(dir + Project + inputs + Gen).getroot().xpath('Paraview')[0].get("CellPot"))
ParTrack=int(etree.parse(dir + Project + inputs + Gen).getroot().xpath('ParTrack')[0].get("value"))
Sym_Contagion=int(etree.parse(dir + Project + inputs + Gen).getroot().xpath('Sym_Contagion')[0].get("value"))
Apo_Contagion=int(etree.parse(dir + Project + inputs + Gen).getroot().xpath('Apo_Contagion')[0].get("value"))
color_threshold=float(etree.parse(dir + Project + inputs + Gen).getroot().xpath('color_threshold')[0].get("value"))
thickness_disp=float(etree.parse(dir + Project + inputs + Gen).getroot().xpath('thickness_disp')[0].get("value"))
thicknessJunction_disp=float(etree.parse(dir + Project + inputs + Gen).getroot().xpath('thicknessJunction_disp')[0].get("value"))
radiusPlasmodesm_disp=float(etree.parse(dir + Project + inputs + Gen).getroot().xpath('radiusPlasmodesm_disp')[0].get("value"))
UniXwalls=int(etree.parse(dir + Project + inputs + Gen).getroot().xpath('UniXwalls')[0].get("value"))
sparseM=int(etree.parse(dir + Project + inputs + Gen).getroot().xpath('sparse')[0].get("value"))

#Import Geometrical data
Plant=etree.parse(dir + Project + inputs + Geom).getroot().xpath('Plant')[0].get("value")
path_geom=etree.parse(dir + Project + inputs + Geom).getroot().xpath('path')[0].get("value")
im_scale=float(etree.parse(dir + Project + inputs + Geom).getroot().xpath('im_scale')[0].get("value"))
Maturityrange=etree.parse(dir + Project + inputs + Geom).getroot().xpath('Maturityrange/Maturity')
Printrange=etree.parse(dir + Project + inputs + Geom).getroot().xpath('Printrange/Print_layer')
Xwalls=float(etree.parse(dir + Project + inputs + Geom).getroot().xpath('Xwalls')[0].get("value")) #Transverse walls or not
PileUp=int(etree.parse(dir + Project + inputs + Geom).getroot().xpath('PileUp')[0].get("value"))
passage_cell_range=etree.parse(dir + Project + inputs + Geom).getroot().xpath('passage_cell_range/passage_cell')
aerenchyma_range=etree.parse(dir + Project + inputs + Geom).getroot().xpath('aerenchyma_range/aerenchyma')
passage_cell_ID=[]
for passage_cell in passage_cell_range:
    passage_cell_ID.append(int(passage_cell.get("id")))
PPP=list()
InterCid=list() #Aerenchyma is classified as intercellular space
for aerenchyma in aerenchyma_range:
    if not int(aerenchyma.get("id"))>9E5:
        InterCid.append(int(aerenchyma.get("id"))) #Cell id starting at 0
    else:
        print('InterCid #'+str(int(aerenchyma.get("id")))+' excluded')
InterC_perim_search=int(etree.parse(dir + Project + inputs + Geom).getroot().xpath('InterC_perim_search')[0].get("value"))
if InterC_perim_search==1:
    InterC_perim1=float(etree.parse(dir + Project + inputs + Geom).getroot().xpath('InterC_perim1')[0].get("value"))
    InterC_perim2=float(etree.parse(dir + Project + inputs + Geom).getroot().xpath('InterC_perim2')[0].get("value"))
    InterC_perim3=float(etree.parse(dir + Project + inputs + Geom).getroot().xpath('InterC_perim3')[0].get("value"))
    InterC_perim4=float(etree.parse(dir + Project + inputs + Geom).getroot().xpath('InterC_perim4')[0].get("value"))
    InterC_perim5=float(etree.parse(dir + Project + inputs + Geom).getroot().xpath('InterC_perim5')[0].get("value"))
kInterC=float(etree.parse(dir + Project + inputs + Geom).getroot().xpath('kInterC')[0].get("value"))
cell_per_layer=zeros((2,1))
cell_per_layer[0][0]=float(etree.parse(dir + Project + inputs + Geom).getroot().xpath('cell_per_layer')[0].get("cortex"))
cell_per_layer[1][0]=float(etree.parse(dir + Project + inputs + Geom).getroot().xpath('cell_per_layer')[0].get("stele"))
thickness=float(etree.parse(dir + Project + inputs + Geom).getroot().xpath('thickness')[0].get("value")) #micron
PD_section=float(etree.parse(dir + Project + inputs + Geom).getroot().xpath('PD_section')[0].get("value")) #micron^2
Xylem_pieces=False
if float(etree.parse(dir + Project + inputs + Geom).getroot().xpath('Xylem_pieces')[0].get("flag"))==1:
    Xylem_pieces=True
t1 = time.perf_counter()
print(t1-t0, "seconds process time")

#Import hormone properties
Degrad1=float(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Hormone_movement/Degradation_constant_H1')[0].get("value")) #Hormone 1 degradation constant (mol degraded / mol-day)
K_MB1=float(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Hormone_movement/MB_H1')[0].get("Partition_coef")) #(-)
dx_MB1=float(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Hormone_movement/MB_H1')[0].get("Thickness")) #(micron)
Diff_MB1=float(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Hormone_movement/MB_H1')[0].get("Diffusivity")) #(cm^2/day)
Diff_PD1=float(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Hormone_movement/Diffusivity_PD_H1')[0].get("value")) #Hormone 1 diffusivity constant (cm^2/day)
Diff_PW1=float(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Hormone_movement/Diffusivity_PW_H1')[0].get("value")) #Hormone 1 diffusivity constant (cm^2/day)
Diff_X1=float(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Hormone_movement/Diffusivity_X_H1')[0].get("value")) #Hormone 1 diffusivity constant (cm^2/day)
D2O1=int(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Hormone_movement/H1_D2O')[0].get("flag")) #Hormone 1 diffusivity constant (cm^2/day)
Active_transport_range=etree.parse(dir + Project + inputs + Horm).getroot().xpath('Hormone_active_transport/carrier_range/carrier')
Sym_source_ini_range=etree.parse(dir + Project + inputs + Horm).getroot().xpath('Sym_Contagion/source_range/Steady-state/source')
Sym_source_transi_range=etree.parse(dir + Project + inputs + Horm).getroot().xpath('Sym_Contagion/source_range/Transient/source')
contact_range=etree.parse(dir + Project + inputs + Horm).getroot().xpath('Contactrange/Contact')
Contact=[]
for contact in contact_range:
    Contact.append(int(contact.get("id")))

#Import cellset data
import xml.etree.ElementTree as ET
#import xml.dom.minidom as DOM
#xml = xml.dom.minidom.parse(xml_fname) # or xml.dom.minidom.parseString(xml_string)
#pretty_xml_as_string = xml.toprettyxml()
tree = etree.parse(dir + 'cellsetdata/' + path_geom) #Maize_Charles\\Maize_pip_cross4.xml') # #Parse literally decrypts the tree element data         SteleOK_high.xml
rootelt = tree.getroot()
Cell2Wall_loop = rootelt.xpath('cells/cell/walls') #Cell2Wall_loop contains cell wall groups info (one group by cell), searched by xpath ("Smart" element identifier)

#Set path
points = rootelt.xpath('walls/wall/points') #points contains the wall elements attributes
Walls_loop = rootelt.xpath('cells/cell/walls/wall') #Walls_loop contains the individual cell to wall associations
Walls_PD = rootelt.xpath('walls/wall')
Cells_loop=rootelt.xpath('cells/cell') #Cells_loop contains the cell attributes
newpath=dir+Project+Output_path+Plant+'/'
print('Outputs in '+newpath)
if not os.path.exists(newpath):
    os.makedirs(newpath)

#Initializes network structure
G = nx.Graph() #Full network

#Creates wall & junction nodes
print('Creating network nodes')
t0 = time.perf_counter()
Nwalls=len(points)
Ncells=len(Cells_loop)
NwallsJun=Nwalls #Will increment at each new junction node
Junction_pos={}
Junction2Wall={}
nJunction2Wall={}
position_junctions=empty((Nwalls,4)) #Coordinates of junctions associated to each wall
position_junctions[:]=NAN
min_x_wall=inf
max_x_wall=0
lengths_ini=zeros((Nwalls,1))
jid=0
for p in points: #Loop on wall elements (groups of points)
    wid= int((p.getparent().get)("id")) #wid records the current wall id number
    xprev=inf
    yprev=inf
    length=0.0 #Calculating total wall length
    for r in p: #Loop on points within the wall element to calculate their average X and Y coordinates 
        x= im_scale*float(r.get("x")) #X coordinate of the point
        y= im_scale*float(r.get("y")) #Y coordinate of the point
        if xprev==inf: #First point
            pos="x"+str(x)+"y"+str(y) #Position of the first point
            position_junctions[wid][0]=x
            position_junctions[wid][1]=y
            if pos in Junction_pos:
                ind=Junction_pos[pos]
                Junction2Wall[ind].append(wid) #Several cell wall ID numbers can correspond to the same X Y coordinate where they meet
                nJunction2Wall[ind]+=1
            else: #New junction node
                Junction_pos[pos]=int(jid)
                Junction2Wall[jid]=[wid] #Saves the cell wall ID number associated to the junction X Y coordinates
                nJunction2Wall[jid]=1
                G.add_node(Nwalls+jid, indice=Nwalls+jid, type="apo", position=(float(x),float(y)), length=0) #Nodes are added at walls junctions (previous nodes corresponded to walls middle points). By default, borderlink is 0, but will be adjusted in next loop
                jid+=1
        else:
            length+=hypot(x-xprev,y-yprev)
        xprev=x
        yprev=y
    #Last point in the wall
    pos="x"+str(x)+"y"+str(y) #Position of the last point
    position_junctions[wid][2]=x
    position_junctions[wid][3]=y
    if pos in Junction_pos: #Get the junction ID
        ind=Junction_pos[pos]
        Junction2Wall[ind].append(wid) #Several cell wall ID numbers can correspond to the same X Y coordinate where they meet
        nJunction2Wall[ind]+=1
    else: #New junction node
        Junction_pos[pos]=int(jid)
        Junction2Wall[jid]=[wid] #Saves the cell wall ID number associated to the junction X Y coordinates
        nJunction2Wall[jid]=1
        G.add_node(Nwalls+jid, indice=Nwalls+jid, type="apo", position=(float(x),float(y)), length=0) #Nodes are added at walls junctions (previous nodes corresponded to walls middle points). By default, borderlink is 0, but will be adjusted in next loop
        jid+=1
    #Second round, identifying the mid-point of the wall
    xprev=inf
    yprev=inf
    length2=0.0 #Calculating the cumulative wall length in order to obtain the exact position of the mid-length of the wall from known total length
    for r in p: #Second loop to catch the true middle position of the wall
        x= im_scale*float(r.get("x")) #X coordinate of the point
        y= im_scale*float(r.get("y")) #Y coordinate of the point
        if not xprev==inf:
            temp1=hypot(x-xprev,y-yprev) #length of the current piece of wall
            if temp1==0:
                print('Warning null wall segment length! wid:',wid,' x, xprev, y, yprev:',x,xprev,y,yprev)
                error('error')
            temp2=length2+temp1-length/2 #Cumulative length along the wall
            if temp2>=0: #If beyond the half length of the wall
                mx=x-(x-xprev)*temp2/temp1 #Middle X coordinate of the wall
                my=y-(y-yprev)*temp2/temp1 #Middle Y coordinate of the wall
                break #End the r in p loop
            length2+=temp1
        xprev=x
        yprev=y
    min_x_wall=min(min_x_wall,mx)
    max_x_wall=max(max_x_wall,mx)
    #Creation of the wall node
    G.add_node(wid, indice=wid, type="apo", position=(mx,my), length=length) #Saving wall attributes for graphical display (id, border, type, X and Y coordinates)
    lengths_ini[wid]=length

NwallsJun=Nwalls+jid
Ntot=NwallsJun+Ncells
position=nx.get_node_attributes(G,'position') #Nodes XY positions (micrometers)

#Junction nodes are pointwise by definition so their length is null, except for junctions at root surface, which are attributed a quarter of the length of each surface neighbouring wall for radial transport 
#lengths=nx.get_node_attributes(G,'length') #Walls lengths (micrometers)
lengths=zeros((NwallsJun,1))
lengths[:Nwalls]=lengths_ini

##Calculation of the cosine of the trigonometric orientation between horizontal and the junction-wall vector (radian)
#cos_angle_wall=empty((Nwalls,2))
#cos_angle_wall[:]=NAN
#for wid in range(Nwalls):
#    cos_angle_wall[wid][0]=(position_junctions[wid][0]-position[wid][0])/(hypot(position_junctions[wid][0]-position[wid][0],position_junctions[wid][1]-position[wid][1])) #Vectors junction1-wall
#    cos_angle_wall[wid][1]=(position_junctions[wid][2]-position[wid][0])/(hypot(position_junctions[wid][2]-position[wid][0],position_junctions[wid][3]-position[wid][1])) #Vectors junction2-wall

#Identifies soil-root interface walls
Borderlink=2*ones((NwallsJun,1))
Borderwall=[] #Soil-root interface wall
Borderaerenchyma=[] #Wall at the surface of aerenchyma
for w in Walls_loop: #Loop on walls, by cell - wall association, hence a wall can be repeated if associated to two cells
    wid= int(w.get("id")) #Wall id number
    Borderlink[wid]-=1
for w in Cell2Wall_loop: #Loop on cells. Cell2Wall_loop contains cell wall groups info (one group by cell)
    cgroup=int(w.getparent().get("group")) #Cell type (1=Exodermis;2=epidermis;3=endodermis;4=cortex;5=stele;16=pericycle)
    for r in w: #w points to the cell walls around the current cell
        wid= int(r.get("id")) #Wall id number
        if Borderlink[wid]==1 and cgroup==2: #Wall node at the interface with soil
            if wid not in Borderwall:
                Borderwall.append(wid)
        elif Borderlink[wid]==1:
            if wid not in Borderaerenchyma:
                Borderaerenchyma.append(wid)
#for wid in range(Nwalls):
    
Borderjunction=[]
jid=0
for Junction, Walls in Junction2Wall.items():
    count=0
    length=0
    for wid in Walls:
        if wid in Borderwall:
            count+=1
            length+=lengths[wid]/4.0
    #if count>2: #Should not happen
    #    print('What the count?')
    if count==2:
        Borderjunction.append(jid+Nwalls)
        Borderlink[jid+Nwalls]=1 #Junction node at the interface with soil
        lengths[jid+Nwalls]=length
    else:
        Borderlink[jid+Nwalls]=0
    jid+=1

Sym_target_range=etree.parse(dir + Project + inputs + Horm).getroot().xpath('Sym_Contagion/target_range/target')
Sym_Target=[]
for target in Sym_target_range:
    Sym_Target.append(int(target.get("id")))
Sym_immune_range=etree.parse(dir + Project + inputs + Horm).getroot().xpath('Sym_Contagion/immune_range/immune')
Sym_Immune=[]
for immune in Sym_immune_range:
    Sym_Immune.append(int(immune.get("id")))
Apo_source_ini_range=etree.parse(dir + Project + inputs + Horm).getroot().xpath('Apo_Contagion/source_range/Steady-state/source')
Apo_source_transi_range=etree.parse(dir + Project + inputs + Horm).getroot().xpath('Apo_Contagion/source_range/Transient/source')
Apo_target_range=etree.parse(dir + Project + inputs + Horm).getroot().xpath('Apo_Contagion/target_range/target')
Apo_Target=[]
for target in Apo_target_range:
    Apo_Target.append(int(target.get("id")))
Apo_immune_range=etree.parse(dir + Project + inputs + Horm).getroot().xpath('Apo_Contagion/immune_range/immune')
Apo_Immune=[]
for immune in Apo_immune_range:
    Apo_Immune.append(int(immune.get("id")))

#Get X and Y for Cell nodes and cell nodes
listsieve=[]
listxyl=[]
listxylwalls=[]
Apo_w_Target=[]
Apo_w_Immune=[]
for w in Cell2Wall_loop: #Loop on cells. Cell2Wall_loop contains cell wall groups info (one group by cell)
	totx=0.0 #Summing up cell walls X positions
	toty=0.0 #Summing up cell walls Y positions
	cellnumber1 = int(w.getparent().get("id")) #Cell ID number
	cgroup=int(w.getparent().get("group")) #Cell type (1=Exodermis;2=epidermis;3=endodermis;4=cortex;5=stele;16=pericycle)
	if cgroup==26:
	    error('Please use the label Columella2 for Companion cells in CellSet, MECHA update needed before using the official Companion cell label')
	    cgroup=24
	div=float(len(w)) #Total number of walls around the current cell
	for r in w: #w points to the cell walls around the current cell
	    wid= int(r.get("id")) #Wall ID number
	    totx += position[wid][0] #Contains the walls average X positions
	    toty += position[wid][1] #Contains the walls average Y positions
	finalx=totx/div #Average cell X position (from the average position of its walls)
	finaly=toty/div #Average cell Y position (from the average position of its walls)
	G.add_node(NwallsJun + cellnumber1, indice=(NwallsJun) + cellnumber1, type="cell", position = (finalx,finaly), cgroup=cgroup) #Adding cell nodes  borderlink=0,
	if cgroup==23: #Phloem sieve tube
	    error('Please use the label "Columella1" for Phloem in CellSet, MECHA update needed in order to use the official CellSet phloem label')
	if cgroup==11 or cgroup==23: #Phloem sieve tube
	    listsieve.append(NwallsJun+cellnumber1)
	elif cgroup==13 or cgroup==19 or cgroup==20: #Xylem vessel
	    listxyl.append(NwallsJun+cellnumber1)
	    for r in w: #w points to the cell walls around the current cell
	        wid= int(r.get("id")) #Wall ID number
	        listxylwalls.append(wid) #ghost walls crossing xylem vessels will appear twice
	if Apo_Contagion==1:
	    if cellnumber1 in Apo_Target:
	        for r in w: #w points to the cell walls around the current cell
	            wid= int(r.get("id")) #Wall ID number
	            if wid not in Apo_w_Target:
	                Apo_w_Target.append(wid)
	    if cellnumber1 in Apo_Immune:
	        for r in w: #w points to the cell walls around the current cell
	            wid= int(r.get("id")) #Wall ID number
	            if wid not in Apo_w_Immune:
	                Apo_w_Immune.append(wid)

Nxyl=len(listxyl)
position=nx.get_node_attributes(G,'position') #Updates nodes XY positions (micrometers)

t1 = time.perf_counter()
print(t1-t0, "seconds process time")

#add Edges
print('Creating network connections')
t0 = time.perf_counter()
lat_dists=zeros((Nwalls,1))
Nmb=0 #Total number of membranes
cellperimeter=np.linspace(0,0,Ncells)
cellarea=np.linspace(0,0,Ncells) #(micron^2)
CellWallsList=[] #Includes both walls & junctions ordered in a consecutive order
for w in Cell2Wall_loop: #Loop on cells. Cell2Wall_loop contains cell wall groups info (one group by cell)
    cellnumber1 = int(w.getparent().get("id")) #Cell ID number
    i=0
    for r in w: #Loop for wall elements around the cell
        wid= int(r.get("id")) #Cell wall ID
        d_vec=array([position[wid][0]-position[NwallsJun+cellnumber1][0],position[wid][1]-position[NwallsJun+cellnumber1][1]])
        dist_cell=hypot(d_vec[0],d_vec[1]) #distance between wall node and cell node (micrometers)
        d_vec/=dist_cell
        lat_dists[wid]+=dist_cell
        if i==0:
            wid0=wid
            wid1=wid
        else: #This algorithm only works if walls are ordered anti-clockwise around the cell center
            wid2=wid #new wall id
            #Find junction closest to wid1
            dist1=hypot(position[wid1][0]-position_junctions[wid2][0],position[wid1][1]-position_junctions[wid2][1])
            dist2=hypot(position[wid1][0]-position_junctions[wid2][2],position[wid1][1]-position_junctions[wid2][3])
            if dist1<dist2:
                j=0
            else:
                j=2
            cellarea[cellnumber1] += (position[wid1][0]+position_junctions[wid2][0+j])*(position[wid1][1]-position_junctions[wid2][1+j]) #Cell area loop (micron^2)
            cellarea[cellnumber1] += (position_junctions[wid2][0+j]+position[wid2][0])*(position_junctions[wid2][1+j]-position[wid2][1]) #Cell area loop (micron^2)
            wid1=wid2
        Nmb+=1
        G.add_edge(NwallsJun + cellnumber1, wid, path='membrane', length=lengths[wid], dist=dist_cell, d_vec=d_vec) #, height=height #Adding all cell to wall connections (edges) #kaqp=kaqp, kw=kw, kmb=kmb, 
        cellperimeter[cellnumber1]+=lengths[wid]
        i+=1
    dist1=hypot(position[wid1][0]-position_junctions[wid0][0],position[wid1][1]-position_junctions[wid0][1])
    dist2=hypot(position[wid1][0]-position_junctions[wid0][2],position[wid1][1]-position_junctions[wid0][3])
    if dist1<dist2:
        j=0
    else:
        j=2
    cellarea[cellnumber1] += (position[wid1][0]+position_junctions[wid0][0+j])*(position[wid1][1]-position_junctions[wid0][1+j]) #Back to the first node
    cellarea[cellnumber1] += (position_junctions[wid0][0+j]+position[wid0][0])*(position_junctions[wid0][1+j]-position[wid0][1]) #Back to the first node
    cellarea[cellnumber1] /= -2.0

Cell_connec=-ones((Ncells,Cell_connec_max),dtype=int) #Connected cells for further ranking
nCell_connec=zeros((Ncells,1),dtype=int) #Quantity of cell to cell connections
for i in range(0, len(Walls_loop)): #Loop on walls, by cell - wall association, hence a wall can be repeated if associated to two cells. Parent structure: Cell/Walls/Wall
           	r1 = Walls_loop[i] #Points to the current wall
           	cellid1 = r1.getparent().getparent().get("id") #Cell1 ID number
           	id1 = r1.get("id") #Wall1 ID number
           	for j in range(i + 1, len(Walls_loop) ): #Loop on cell-wall associations that are further down in the list
           	    r2 = Walls_loop[j] #Points to the further down wall in the list of cell-wall associations
           	    cellid2 = r2.getparent().getparent().get("id") #Cell2 ID number
           	    id2 = r2.get("id") #Wall2 ID number
           	    if id1 == id2: #If walls 1 and 2 are the same, then cells 1 and 2 are connected by plasmodesmata
           	        d_vec=array([position[NwallsJun+int(cellid2)][0]-position[NwallsJun+int(cellid1)][0],position[NwallsJun+int(cellid2)][1]-position[NwallsJun+int(cellid1)][1]])
           	        dist_cell=hypot(d_vec[0],d_vec[1]) #distance between wall node and cell node (micrometers)
           	        d_vec/=dist_cell
           	        G.add_edge(NwallsJun + int(cellid1), NwallsJun + int(cellid2), path='plasmodesmata', length=lengths[int(id1)], d_vec=d_vec) #, height=height #Adding all cell to cell connections (edges) #kpl=kpl, 
           	        Cell_connec[int(cellid1)][nCell_connec[int(cellid1)]]=int(cellid2)
           	        nCell_connec[int(cellid1)]+=1
           	        Cell_connec[int(cellid2)][nCell_connec[int(cellid2)]]=int(cellid1)
           	        nCell_connec[int(cellid2)]+=1

jid=0
for Junction, Walls in Junction2Wall.items(): #Loop on junctions between walls
    for wid in Walls: #Walls is the list of cell walls ID meeting at the junction pos
        d_vec=array([position[wid][0]-position[jid+Nwalls][0],position[wid][1]-position[jid+Nwalls][1]])
        dist_wall=hypot(d_vec[0],d_vec[1]) #distance between wall node and cell node (micrometers)
        d_vec/=dist_wall #As compared to lat_dist, dist_wall underestimates the actual path length between the wall and the junction. dist_wall is rather a straight distance.
        G.add_edge(jid+Nwalls, int(wid), path='wall', length=lengths[int(wid)]/2, lat_dist=lat_dists[int(wid)][0], d_vec=d_vec, dist_wall=dist_wall) #Adding junction to wall connections (edges)
    jid+=1

#And calculation of the centre of gravity of the endodermis
x_grav=0.0 # (micrometers)
y_grav=0.0 # (micrometers)
n_cell_endo=0 #Counting the total number of cells in the endodermis
for w in Cell2Wall_loop: #Loop on cells. Cell2Wall_loop contains cell wall groups info (one group by cell)
	cellnumber1 = int(w.getparent().get("id")) #Cell ID number
	if G.node[NwallsJun + cellnumber1]['cgroup']==3: #Endodermis
	        x_grav+=position[NwallsJun + cellnumber1][0]
	        y_grav+=position[NwallsJun + cellnumber1][1]
	        n_cell_endo+=1
if n_cell_endo>0:
    x_grav/=n_cell_endo
    y_grav/=n_cell_endo
else:
    x_grav=nan
    y_grav=nan

t1 = time.perf_counter()
print(t1-t0, "seconds process time")

#Calculation of a and b parameter for cortex AQP activity radial distribution
print('Identifying cell layers')
t0 = time.perf_counter()
#<group id="2" name="epidermis" />      is epidermis
#<group id="1" name="general" />        is exodermis
#<group id="4" name="cortex" />         is cortex
#<group id="3" name="endodermis" />     is endodermis
#<group id="16" name="pericycle" />     is pericycle
#<group id="5" name="stele" />          is stelar parenchyma
#<group id="11" name="columella1" />    is phloem sieve tube
#<group id="13" name="columella3" />    is xylem
#<group id="12" name="columella2" />    is companion cell
#<group id="19" name="Protoxylem" />    is xylem
#<group id="20" name="Metaxylem" />     is xylem
#<group id="21" name="XylemPolePericyle" /> is pericycle
#<group id="23" name="Phloem" />        is phloem
#<group id="26" name="CompanionCell" /> is companion cell changed to 24...

Cell_rank=zeros((Ncells,1)) #Ranking of cells (1=Exodermis, 2=Epidermis, 3=Endodermis, 4*=Cortex, 5*=Stele, 11=Phloem sieve tube, 12=Companion cell, 13=Xylem, 16=Pericycle), stars are replaced by the ranking within cortical cells and stele cells
Layer_dist=zeros((62,1)) #Average cell layers distances from center of gravity, by cells ranking 
nLayer=zeros((62,1)) #Total number of cells in each rank (indices follow ranking numbers)
xyl_dist=[] #List of distances between xylem and cross-section centre
#angle_dist_endo_grav=array([-4,0]) #array of distances and angles between endo cells and grav. Initializing the array with values that will eventualy be deleted
#angle_dist_exo_grav=array([-4,0]) #array of distances and angles between exo cells and grav. Initializing the array with values that will eventualy be deleted
for w in Cell2Wall_loop: #Loop on cells. Cell2Wall_loop contains cell wall groups info (one group by cell)
	cellnumber1 = int(w.getparent().get("id")) #Cell ID number
	celltype=G.node[NwallsJun + cellnumber1]['cgroup'] #Cell type
	if celltype==19 or celltype==20: #Proto- and Meta-xylem in new Cellset version
	    celltype=13
	elif celltype==21: #Xylem pole pericycle in new Cellset version
	    celltype=16
	elif celltype==23: #Phloem in new Cellset version
	    celltype=11
	elif celltype==26 or celltype==24: #Companion cell in new Cellset version
	    celltype=12
	Cell_rank[cellnumber1]=celltype #Later on, cell types 4 and 5 will be updated to account for their ranking within the cortex / stele
	x_cell=position[NwallsJun + cellnumber1][0] #Cell position (micrometers)
	y_cell=position[NwallsJun + cellnumber1][1]
	dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	Layer_dist[celltype]+=dist
	nLayer[celltype]+=1
	if celltype==13:
	    xyl_dist.append([dist])
if len(xyl_dist)>0:
    xyl80_dist=percentile(xyl_dist, 80)
else:
    xyl80_dist=nan

if nLayer[16]==0: #If there is no labelled pericycle
    stele_connec_rank=3 #Endodermis connected to stele cells
else:
    stele_connec_rank=16 #Pericycle connected to stele cells
if nLayer[1]==0: #If there is no labelled exodermis
    outercortex_connec_rank=2 #Cortex connected to epidermis cells
else:
    outercortex_connec_rank=1 #Cortex connected to exodermis cells
if InterC_perim_search==1:
    rank_cellperimeters_in=linspace(nan,nan,Ncellperimeters)
    rank_cellperimeters_out=linspace(nan,nan,Ncellperimeters)
listprotosieve=[]
mincid=99999
Layer_dist[16]=0
nLayer[16]=0
for w in Cell2Wall_loop: #Loop on cells. Cell2Wall_loop contains cell wall groups info (one group by cell)
    cellnumber1 = int(w.getparent().get("id")) #Cell ID number #celltype=G.node[NwallsJun + cellnumber1]['cgroup']
    celltype=Cell_rank[cellnumber1] #Cell types 4 and 5 updated to account for their ranking within the cortex / stele
    if celltype==16: #Pericycle
        temp=Cell_rank[Cell_connec[cellnumber1][0:nCell_connec[cellnumber1][0]]]
        if any(temp==5) or any(temp==11) or any(temp==12) or any(temp==13) or any(temp==50): #Cell to cell connection with endodermis
            Cell_rank[cellnumber1]=16 #pericycle ranks now span 16 to 18 instead of 16
            x_cell=position[NwallsJun + cellnumber1][0] #Cell position (micrometers)
            y_cell=position[NwallsJun + cellnumber1][1]
            dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
            Layer_dist[16]+=dist
            nLayer[16]+=1
        elif any(temp==3): #Cell to cell connection with endodermis
            Cell_rank[cellnumber1]=18 #pericycle ranks now span 16 to 18 instead of 16
            x_cell=position[NwallsJun + cellnumber1][0] #Cell position (micrometers)
            y_cell=position[NwallsJun + cellnumber1][1]
            dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
            Layer_dist[18]+=dist
            nLayer[18]+=1
        elif all(np.logical_or(np.logical_or(temp==16,temp==17),temp==18)):
            Cell_rank[cellnumber1]=17 #pericycle ranks now span 16 to 18 instead of 16
            x_cell=position[NwallsJun + cellnumber1][0] #Cell position (micrometers)
            y_cell=position[NwallsJun + cellnumber1][1]
            dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
            Layer_dist[17]+=dist
            nLayer[17]+=1
        else:
            error('Pericycle cell falls in no category')
    elif celltype==4: #Cortex
	    if any(Cell_rank[Cell_connec[cellnumber1][0:nCell_connec[cellnumber1][0]]]==3): #Cell to cell connection with endodermis
	        Cell_rank[cellnumber1]=25 #Cortex ranks now span 25 to 49 instead of 40 to 49
	        x_cell=position[NwallsJun + cellnumber1][0] #Cell position (micrometers)
	        y_cell=position[NwallsJun + cellnumber1][1]
	        dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	        Layer_dist[25]+=dist
	        nLayer[25]+=1
	        if InterC_perim_search==1:
	            rank_cellperimeters_in[int(nLayer[25]-1)]=cellperimeter[cellnumber1]
	            if cellperimeter[cellnumber1]<InterC_perim1:
	                InterCid.append(cellnumber1) #Cell id starting at 0
	    elif any(Cell_rank[Cell_connec[cellnumber1][0:nCell_connec[cellnumber1][0]]]==outercortex_connec_rank): #Cell to cell connection with exodermis
	        Cell_rank[cellnumber1]=49
	        x_cell=position[NwallsJun + cellnumber1][0] #Cell position (micrometers)
	        y_cell=position[NwallsJun + cellnumber1][1]
	        dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	        Layer_dist[49]+=dist
	        nLayer[49]+=1
	        if InterC_perim_search==1:
	            rank_cellperimeters_out[int(nLayer[49]-1)]=cellperimeter[cellnumber1]
	            if cellperimeter[cellnumber1]<InterC_perim5:
	                InterCid.append(cellnumber1) #Cell id starting at 0
    elif celltype==5 or celltype==11 or celltype==12 or celltype==13: #Stele
	    if any(Cell_rank[Cell_connec[cellnumber1][0:nCell_connec[cellnumber1][0]]]==stele_connec_rank): #Cell to cell connection with pericycle
	        Cell_rank[cellnumber1]=50
	        x_cell=position[NwallsJun + cellnumber1][0] #Cell position (micrometers)
	        y_cell=position[NwallsJun + cellnumber1][1]
	        dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	        Layer_dist[50]+=dist
	        nLayer[50]+=1
	        if G.node[NwallsJun + cellnumber1]['cgroup']==11 or G.node[NwallsJun + cellnumber1]['cgroup']==23:
	           listprotosieve.append(NwallsJun + cellnumber1)
Nsieve=len(listsieve)
Nprotosieve=len(listprotosieve)

if InterC_perim_search==1:
    cortex_cellperimeters_in=rank_cellperimeters_in #Inner part of the cortex (close to endodermis)
    cortex_cellperimeters_out=rank_cellperimeters_out #Outer part of cortex
for i in range(12):
    if InterC_perim_search==1:
        rank_cellperimeters_in=linspace(nan,nan,Ncellperimeters)
        rank_cellperimeters_out=linspace(nan,nan,Ncellperimeters)
    for w in Cell2Wall_loop: #Loop on cells. Cell2Wall_loop contains cell wall groups info (one group by cell)
        cellnumber1 = int(w.getparent().get("id")) #Cell ID number
        celltype=Cell_rank[cellnumber1] #Cell types 4 and 5 updated to account for their ranking within the cortex / stele
        if celltype==4 and i<12: #Cortex # if i<12: #Within 12 layers of cortical sides
            temp=Cell_connec[cellnumber1][0:nCell_connec[cellnumber1][0]]
            AloneC=True
            for i_connec in range(len(temp)):
                if not temp[i_connec] in InterCid:
                    AloneC=False #Isolated cell only in contact with intercellular spaces
            done=False
            for i_connec in range(len(temp)):
                if Cell_rank[temp[i_connec]]==(25+i) and not done: #Cell to cell connection with endodermis
                    if (not temp[i_connec] in InterCid) or AloneC:    
                        Cell_rank[cellnumber1]=26+i
                        x_cell=position[NwallsJun + cellnumber1][0] #Cell position (micrometers)
                        y_cell=position[NwallsJun + cellnumber1][1]
                        dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
                        Layer_dist[26+i]+=dist
                        nLayer[26+i]+=1
                        done=True
                        if InterC_perim_search==1:
                            rank_cellperimeters_in[int(nLayer[26+i]-1)]=cellperimeter[cellnumber1]
                            if i==0 and cellperimeter[cellnumber1]<InterC_perim2:
                                InterCid.append(cellnumber1)
                            elif i==1 and cellperimeter[cellnumber1]<InterC_perim3:
                                InterCid.append(cellnumber1)
                            elif i==2 and cellperimeter[cellnumber1]<InterC_perim4:
                                InterCid.append(cellnumber1)
                            elif i>2 and cellperimeter[cellnumber1]<InterC_perim5:
                                InterCid.append(cellnumber1)
                elif Cell_rank[temp[i_connec]]==(49-i) and not done: #Cell to cell connection with exodermis
                    if (not temp[i_connec] in InterCid) or AloneC:
                        Cell_rank[cellnumber1]=48-i
                        x_cell=position[NwallsJun + cellnumber1][0] #Cell position (micrometers)
                        y_cell=position[NwallsJun + cellnumber1][1]
                        dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
                        Layer_dist[48-i]+=dist
                        nLayer[48-i]+=1
                        done=True
                        if InterC_perim_search==1:
                            rank_cellperimeters_out[int(nLayer[48-i]-1)]=cellperimeter[cellnumber1]
                            if cellperimeter[cellnumber1]<InterC_perim5:
                                InterCid.append(cellnumber1)
        elif celltype==5 or celltype==11 or celltype==12 or celltype==13: #Stele
            if i<10:
                temp=Cell_connec[cellnumber1][0:nCell_connec[cellnumber1][0]]
                done=False
                for i_connec in range(len(temp)):
                    if Cell_rank[temp[i_connec]]==(50+i) and not done:
                        if (not temp[i_connec]+NwallsJun in listxyl):
                            #if any(Cell_rank[temp]==(50+i)): #Cell to cell connection with pericycle
                            Cell_rank[cellnumber1]=51+i
                            x_cell=position[NwallsJun + cellnumber1][0] #Cell position (micrometers)
                            y_cell=position[NwallsJun + cellnumber1][1]
                            dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
                            Layer_dist[51+i]+=dist
                            nLayer[51+i]+=1
                            done=True
            else: #No more than 11 stele cell layers
                Cell_rank[cellnumber1]=61
                x_cell=position[NwallsJun + cellnumber1][0] #Cell position (micrometers)
                y_cell=position[NwallsJun + cellnumber1][1]
                dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
                Layer_dist[61]+=dist
                nLayer[61]+=1
    if i<12:
        if InterC_perim_search==1:
            cortex_cellperimeters_in=vstack((cortex_cellperimeters_in,rank_cellperimeters_in))
            cortex_cellperimeters_out=vstack((rank_cellperimeters_out,cortex_cellperimeters_out))
if InterC_perim_search==1:
    cortex_cellperimeters=vstack((cortex_cellperimeters_in,cortex_cellperimeters_out))
#InterCid=InterCid[1:]

#Calculating cell surfaces at tissue interfaces (total and interfacing with a cell that is not an intercellular space)
indice=nx.get_node_attributes(G,'indice') #Node indices (walls, junctions and cells)
Length_outer_cortex_tot=0.0 #Total cross-section membrane length at the interface between exodermis and cortex
Length_cortex_cortex_tot=0.0 #Total cross-section membrane length at the interface between cortex and cortex
Length_cortex_endo_tot=0.0 #Total cross-section membrane length at the interface between cortex and endodermis
Length_outer_cortex_nospace=0.0 #Cross-section membrane length at the interface between exodermis and cortex not including interfaces with intercellular spaces
Length_cortex_cortex_nospace=0.0 #Cross-section membrane length at the interface between exodermis and cortex not including interfaces with intercellular spaces
Length_cortex_endo_nospace=0.0 #Cross-section membrane length at the interface between exodermis and cortex not including interfaces with intercellular spaces
for node, edges in G.adjacency_iter() :
    i=indice[node] #Node ID number
    if i>=NwallsJun: #Cell
        if G.node[i]['cgroup']==16 or G.node[i]['cgroup']==21:
            for neighbour, eattr in edges.items(): #Loop on connections (edges)
                if eattr['path'] == "plasmodesmata" and (G.node[indice[neighbour]]['cgroup']==11 or G.node[indice[neighbour]]['cgroup']==23): #Plasmodesmata connection  #eattr is the edge attribute (i.e. connection type)
                    PPP.append(i-NwallsJun)
        elif G.node[i]['cgroup']==outercortex_connec_rank or G.node[i]['cgroup']==4 or G.node[i]['cgroup']==3: #exodermis or cortex or endodermis (or epidermis if there is no exodermis)
            if i-NwallsJun not in InterCid: #The loop focuses on exo, cortex and endodermis cells that are not intercellular spaces
                for neighbour, eattr in edges.items(): #Loop on connections (edges)
                    if eattr['path'] == "plasmodesmata": #Plasmodesmata connection  #eattr is the edge attribute (i.e. connection type)
                        j = (indice[neighbour]) #neighbouring node number
                        l_membrane=eattr['length']
                        if (G.node[i]['cgroup']==outercortex_connec_rank and G.node[j]['cgroup']==4) or (G.node[j]['cgroup']==outercortex_connec_rank and G.node[i]['cgroup']==4):#Exodermis to cortex cell or vice versa (epidermis if no exodermis exists)
                            Length_outer_cortex_tot+=l_membrane
                            if j-NwallsJun not in InterCid:
                                Length_outer_cortex_nospace+=l_membrane
                        elif (G.node[i]['cgroup']==4 and G.node[j]['cgroup']==4):#Cortex to cortex cell
                            Length_cortex_cortex_tot+=l_membrane
                            if j-NwallsJun not in InterCid:
                                Length_cortex_cortex_nospace+=l_membrane
                        elif (G.node[i]['cgroup']==3 and G.node[j]['cgroup']==4) or (G.node[j]['cgroup']==3 and G.node[i]['cgroup']==4):#Cortex to endodermis cell or vice versa
                            Length_cortex_endo_tot+=l_membrane
                            if j-NwallsJun not in InterCid:
                                Length_cortex_endo_nospace+=l_membrane

for i in range(62): #Finalizing distance averaging
    if nLayer[i]>0:
        Layer_dist[i]=Layer_dist[i]/nLayer[i]

#Discretization based on effective cell layering
r_discret=array([0])
j=0 #Counts cell layers
k=0 #Counts tissue types
rank2row=zeros((62,1))*nan
Layer_dist2=array([0])
for i in range(61, 49, -1): #rank2row list, Stele
    if nLayer[i]>0:
        rank2row[i]=j
        Layer_dist2=vstack((Layer_dist2,Layer_dist[i]))
        j+=1
k+=1
r_discret=vstack((r_discret,j))
if nLayer[16]>0:
    for i in range(16,19):
        if nLayer[i]>0:
            rank2row[i]=j #Pericycle
            Layer_dist2=vstack((Layer_dist2,Layer_dist[i]))
            j+=1
    k+=1
    r_discret=vstack((r_discret,j-sum(r_discret[1:k])))
rank2row[3]=j #2 rows for endodermis (inner and outer) + 2 rows for passage cells in between
Layer_dist2=vstack((Layer_dist2,Layer_dist[3]))
Layer_dist2=vstack((Layer_dist2,Layer_dist[3]))
Layer_dist2=vstack((Layer_dist2,Layer_dist[3]))
Layer_dist2=vstack((Layer_dist2,Layer_dist[3]))
j+=4
k+=1
r_discret=vstack((r_discret,j-sum(r_discret[1:k])))

i1=25
nLayer_ref=nLayer[i1]
Layer_dist_ref=Layer_dist[i1]
rank2row[i1]=j
Layer_dist2=vstack((Layer_dist2,Layer_dist[i1]))
j+=1 #number of the row that was just added
i1=26
if len(aerenchyma_range)>3:
    ratio_complete=0.40
else:
    ratio_complete=0.75
while i1<50: #Cortex
    if nLayer[i1]>ratio_complete*nLayer_ref: #Likely complete rank/layer
        rank2row[i1]=j
        Layer_dist2=vstack((Layer_dist2,Layer_dist[i1]))
        j+=1
        nLayer_ref=nLayer[i1]
        Layer_dist_ref=Layer_dist[i1]
        i1+=1
    elif nLayer[i1]>0: #Likely incomplete rank/layer i1
        #Find next non-empty rank (i1+i2)
        i2=i1+1
        while nLayer[i2]==0:
            i2+=1
        #Check if incomplete layers/ranks i1 and i2 would constitute a full layer
        if nLayer[i1]+nLayer[i2]>ratio_complete*nLayer_ref: #They do
            #Does layer/rank i2 constitute a full layer?
            if nLayer[i2]>ratio_complete*nLayer_ref: #Yes
                #Rank i1 is added to the closest layer
                if abs(Layer_dist[i1]-Layer_dist[i2])<abs(Layer_dist[i1]-Layer_dist_ref): #i2 is the closest layer
                    rank2row[i1]=j
                    rank2row[i2]=j
                    Layer_dist_avg=(nLayer[i1]*Layer_dist[i1]+nLayer[i2]*Layer_dist[i2])/(nLayer[i1]+nLayer[i2])
                    Layer_dist2=vstack((Layer_dist2,Layer_dist_avg))
                    j+=1
                else: #i1-1 is the closest layer to rank i1
                    rank2row[i1]=j-1
                    Layer_dist_avg=(nLayer[i1]*Layer_dist[i1]+nLayer_ref*Layer_dist_ref)/(nLayer[i1]+nLayer_ref)
                    Layer_dist2[j-1]=Layer_dist_avg
                    rank2row[i2]=j
                    Layer_dist2=vstack((Layer_dist2,Layer_dist[i2]))
                    j+=1
            else: #Rank i2 does not constitute a full layer, then i1 and i2 form a layer together
                rank2row[i1]=j
                rank2row[i2]=j
                Layer_dist_avg=(nLayer[i1]*Layer_dist[i1]+nLayer[i2]*Layer_dist[i2])/(nLayer[i1]+nLayer[i2])
                Layer_dist2=vstack((Layer_dist2,Layer_dist_avg))
                j+=1
            i1=i2+1
        else: #Ranks i1 and i2 don't constitute a full layer, then i1 and i2 can either be both added to the same layer, or to separate layers
            #Ranks i1 and i2 are added to the closest layer
            if nLayer[i2+1]>0 and i2+1<50:#the next full layer needs to be part of the cortex
                if abs(Layer_dist[i1]-Layer_dist[i2+1])<abs(Layer_dist[i1]-Layer_dist_ref): #i2+1 is the closest layer, so i1 and i2 are both added with i2+1
                    rank2row[i1]=j
                    rank2row[i2]=j
                    rank2row[i2+1]=j
                    Layer_dist_avg=(nLayer[i1]*Layer_dist[i1]+nLayer[i2]*Layer_dist[i2]+nLayer[i2+1]*Layer_dist[i2+1])/(nLayer[i1]+nLayer[i2]+nLayer[i2+1])
                    Layer_dist2=vstack((Layer_dist2,Layer_dist_avg))
                    j+=1
                else:#i1 is closer to i1-1 than i2+1
                    rank2row[i1]=j #Correction 18/04: j replaces j-1
                    if abs(Layer_dist[i2]-Layer_dist[i2+1])<abs(Layer_dist[i2]-Layer_dist_ref): #i2 is closer to i2+1 than to i1-1
                        Layer_dist_avg=(nLayer[i1]*Layer_dist[i1]+nLayer_ref*Layer_dist_ref)/(nLayer[i1]+nLayer_ref)
                        Layer_dist2[j]=Layer_dist_avg #Correction 18/04: j replaces j-1
                        rank2row[i2]=j
                        rank2row[i2+1]=j
                        Layer_dist_avg=(nLayer[i2]*Layer_dist[i2]+nLayer[i2+1]*Layer_dist[i2+1])/(nLayer[i2]+nLayer[i2+1])
                        Layer_dist2=vstack((Layer_dist2,Layer_dist_avg))
                        j+=1
                    else: #i2 is closer to i1-1 than to i2+1
                        rank2row[i2]=j-1
                        Layer_dist_avg=(nLayer[i1]*Layer_dist[i1]+nLayer[i2]*Layer_dist[i2]+nLayer_ref*Layer_dist_ref)/(nLayer[i1]+nLayer[i2]+nLayer_ref)
                        Layer_dist2[j]=Layer_dist_avg #Correction 18/04: j replaces j-1
                        rank2row[i2+1]=j
                        Layer_dist2=vstack((Layer_dist2,Layer_dist[i2+1]))
                        j+=1
            else: #We just merge i1 and i2
                rank2row[i1]=j
                rank2row[i2]=j
                Layer_dist_avg=(nLayer[i1]*Layer_dist[i1]+nLayer[i2]*Layer_dist[i2])/(nLayer[i1]+nLayer[i2])
                Layer_dist2=vstack((Layer_dist2,Layer_dist_avg))
                j+=1
            i1=i2+2
    elif nLayer[i1]==0: #No more cortex layer in the middle (likely there was only one or two of them)
        i1+=1
k+=1
r_discret=vstack((r_discret,j-sum(r_discret[1:k])))
if nLayer[1]>0: #if there is an exodermal layer
    rank2row[1]=j #2 rows for exodermis (inner and outer)
    Layer_dist2=vstack((Layer_dist2,Layer_dist[1]))
    Layer_dist2=vstack((Layer_dist2,Layer_dist[1]))
    j+=2
    k+=1
    r_discret=vstack((r_discret,j-sum(r_discret[1:k])))
rank2row[2]=j #Epidermis
Layer_dist2=vstack((Layer_dist2,Layer_dist[2]))
j+=1
k+=1
r_discret=vstack((r_discret,j-sum(r_discret[1:k])))
r_discret[0]=j
Layer_dist2=Layer_dist2[1:]
row_outercortex=rank2row[1]-1
if isnan(row_outercortex): #No exodermal layer
    row_outercortex=rank2row[2]-1
row_outercortex=int(row_outercortex) 

#Then min and max distances to cross-section centre
ratio_cortex=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('ratio_cortex')[0].get("value"))
if not ratio_cortex==1:
    dmax_cortex=0.0 #For gradient of relative AQP distribution
    dmin_cortex=inf #For gradient of relative AQP distribution
davg_epi=0.0 #avg distance of between exodermis and centre of gravity (micrometers)
dist_grav=zeros((Nwalls,1)) #distance between node and cross-section center of gravity (micrometers)
temp=0.0 #Counting the number of membranes from epidermis
for w in Cell2Wall_loop: #Loop on cells. Cell2Wall_loop contains cell wall groups info (one group by cell)
	cellnumber1 = int(w.getparent().get("id")) #Cell ID number
	for r in w: #Loop for wall elements around the cell
		wid= int(r.get("id")) #Cell wall ID
		dist_grav[wid]=sqrt(square(position[wid][0]-x_grav)+square(position[wid][1]-y_grav)) #Distance between membrane and cross section centre of gravity (micrometers)
		if G.node[NwallsJun + cellnumber1]['cgroup']==4 and not ratio_cortex==1: #Cortex
		    dmax_cortex=max(dmax_cortex,dist_grav[wid])
		    dmin_cortex=min(dmin_cortex,dist_grav[wid])
		elif G.node[NwallsJun + cellnumber1]['cgroup']==2: #Epidermis
		    davg_epi+=dist_grav[wid]
		    temp+=1.0
if temp>0:
    davg_epi/=temp #Last step of averaging (note that we take both inner and outer membranes into account in the averaging)
    perimeter=2*pi*davg_epi*1.0E-04 #(cm)
else:
    davg_epi=nan
    perimeter=nan
t1 = time.perf_counter()
print(t1-t0, "seconds process time")

#Calculating cell borders accounting for wall thickness
if Paraview==1 or ParTrack==1:# or Apo_Contagion>0 or Sym_Contagion>0:
    print('Preparing geometrical properties for Paraview')
    t0 = time.perf_counter()
    ThickWalls=[] #Was created to display 3D membranes
    nThickWalls=zeros((2*Nwalls,1)) #This will save how many new junction wall ID were already saved in connections to new wall ID (the vector is a little too long)
    ThickWallsX=[] #Same as Thickwalls except that it includes extra info about borderlink nodes, was created to display walls
    Wall2NewWall=empty((Nwalls,2))
    Wall2NewWall[:]=NAN
    nWall2NewWall=zeros((Nwalls,1))
    Wall2NewWallX=empty((NwallsJun,8)) #This one also includes "junction" to "new junctions" ID
    Wall2NewWallX[:]=NAN
    nWall2NewWallX=zeros((NwallsJun,1))
    ThickWallPolygonX=empty((Nwalls*2,4))
    ThickWallPolygonX[:]=NAN
    nThickWallPolygonX=zeros((Nwalls*2,1))
    Wall2Cell=empty((Nwalls,2))
    Wall2Cell[:]=NAN
    nWall2Cell=zeros((Nwalls,1))
    Junction2Wall2Cell=empty((NwallsJun-Nwalls,12))
    Junction2Wall2Cell[:]=NAN
    nJunction2Wall2Cell=zeros((NwallsJun-Nwalls,1))
    Junction2Wall2=empty((NwallsJun-Nwalls,12)) #This Junction2Wall will be incomplete because it needs to have the same dimension as Junction2Wall2Cell
    Junction2Wall2[:]=NAN
    nJunction2Wall2=zeros((NwallsJun-Nwalls,1))
    Wall2Junction=empty((Nwalls,2))
    Wall2Junction[:]=NAN
    nWall2Junction=zeros((Nwalls,1))
    Cell2ThickWalls=empty((Ncells,maxCell2ThickWalls))
    Cell2ThickWalls[:]=NAN
    nCell2ThickWalls=zeros((Ncells,1))
    r_rel=empty((Nwalls,1))
    r_rel[:]=NAN
    x_rel=empty((NwallsJun+Ncells,1))
    x_rel[:]=NAN
    L_diff=((abs(float(Layer_dist[3]-Layer_dist[2])*1.0E-04),abs(float(Layer_dist[3] - xyl80_dist)*1.0E-04))) #Diffusion lengths (cm)
    twpid=0 #Thick wall point ID
    twpidX=0
    for node, edges in G.adjacency_iter():
        wid=indice[node]
        if wid<Nwalls: #wall that is not a junction (connected to a cell)
            for neighbour, eattr in edges.items(): #Loop on connections (edges)
                cid=int(indice[neighbour])
                if G.node[cid]['type']=='cell':
                    if not position[cid][0]==position[wid][0]:
                        slopeCG=(position[cid][1]-position[wid][1])/(position[cid][0]-position[wid][0]) #slope of the line connecting the wall node to the center of gravity of the cell
                    else:
                        slopeCG=inf
                    x=position[wid][0]+cos(arctan(slopeCG))*thickness_disp/2*sign(position[cid][0]-position[wid][0]) #-> new position of the wall node on the current cell side
                    y=position[wid][1]+sin(arctan(slopeCG))*thickness_disp/2*sign(position[cid][0]-position[wid][0])
                    ThickWalls.append(array((twpid,wid,cid,x,y,inf,inf,Borderlink[wid]))) #Adds the thick wall node ID, its parent node ID, the associated cell ID, the new coordinates in X and Y, and the two neighbouring new junction walls IDs (not known at this point in the loop)     G.node[i]['borderlink']
                    Cell2ThickWalls[cid-NwallsJun][int(nCell2ThickWalls[cid-NwallsJun])]=twpid
                    nCell2ThickWalls[cid-NwallsJun]+=1
                    ThickWallsX.append((twpidX,x,y,wid,cid)) #Adds the thick wall node ID, the new coordinates in X and Y, its parent node ID, the associated cell ID
                    ThickWallPolygonX[2*wid][int(nThickWallPolygonX[2*wid])]=twpidX #This one is for "polygon 1"
                    ThickWallPolygonX[2*wid+1][int(nThickWallPolygonX[2*wid+1])]=twpidX #This one for "polygon 2" (wid new thick nodes are included in two polygons)
                    nThickWallPolygonX[2*wid]+=1
                    nThickWallPolygonX[2*wid+1]+=1
                    Wall2Cell[wid][int(nWall2Cell[wid])]=cid
                    nWall2Cell[wid]+=1
                    Wall2NewWall[wid][int(nWall2NewWall[wid])]=twpid
                    nWall2NewWall[wid]+=1  #the count in nWall2NewWall is actually the same as nWall2Cell 
                    Wall2NewWallX[wid][int(nWall2NewWallX[wid])]=twpidX
                    nWall2NewWallX[wid]+=1
                    twpid+=1
                    twpidX+=1
                    if Borderlink[wid]==1: #G.node[wid]['borderlink']==1:
                        x=position[wid][0]-cos(arctan(slopeCG))*thickness_disp/2*sign(position[cid][0]-position[wid][0]) #-> new position of the wall node opposite to the current cell side
                        y=position[wid][1]-sin(arctan(slopeCG))*thickness_disp/2*sign(position[cid][0]-position[wid][0])
                        ThickWallsX.append((twpidX,x,y,wid,inf))
                        ThickWallPolygonX[2*wid][int(nThickWallPolygonX[2*wid])]=twpidX #This one is for "polygon 1"
                        ThickWallPolygonX[2*wid+1][int(nThickWallPolygonX[2*wid+1])]=twpidX #This one for "polygon 2" (wid new thick nodes are included in two polygons)
                        nThickWallPolygonX[2*wid]+=1
                        nThickWallPolygonX[2*wid+1]+=1
                        Wall2NewWallX[wid][int(nWall2NewWallX[wid])]=twpidX
                        nWall2NewWallX[wid]+=1
                        twpidX+=1
                elif G.node[neighbour]['type']=='apo': #Node j is a junction
                    Wall2Junction[wid][int(nWall2Junction[wid])]=indice[neighbour]
                    nWall2Junction[wid]+=1
            cid1=Wall2Cell[wid,0]
            cid2=Wall2Cell[wid,1]
            rank1=Cell_rank[int(cid1-NwallsJun)]
            row1=rank2row[int(rank1)]
            if not isnan(cid2):
                rank2=Cell_rank[int(cid2-NwallsJun)]
                row2=rank2row[int(rank2)]
            if not isnan(xyl80_dist):
                if row1 <= rank2row[3] and row2 <= rank2row[3]: #Row on the stelar side
                    if not isnan(cid2):
                        rad_pos=-((Layer_dist2[int(row1)]+Layer_dist2[int(row2)])/2 - xyl80_dist)/(Layer_dist[3] - xyl80_dist) #Radial position of the wall layer relative to the cortical layer closest to tthe xylem vessels (zero) and to the endodermis (1)
                    else:
                        rad_pos=-(Layer_dist2[int(row1)] - xyl80_dist)/(Layer_dist[3] - xyl80_dist) #Radial position of the wall layer relative to the cortical layer closest to the xylem vessels (zero) and to the endodermis (1)
                    r_rel[wid]=max(min(rad_pos,-0.00001),-1)
                else: #Positive r_rel means a wall on the cortical side of the endodermis
                    if not isnan(cid2):
                        rad_pos=(Layer_dist[2]-(Layer_dist2[int(row1)]+Layer_dist2[int(row2)])/2)/(Layer_dist[2]-Layer_dist[3]) #Radial position of the wall layer relative to the cortical layer closest to the epidermis (zero) and to the endodermis (1)
                    else:
                        rad_pos=(Layer_dist[2]-Layer_dist2[int(row1)])/(Layer_dist[2]-Layer_dist[3]) #Radial position of the wall layer relative to the cortical layer closest to the epidermis (zero) and to the endodermis (1)
                    r_rel[wid]=min(max(rad_pos,0.00001),1)
        x_rel[wid]=(position[wid][0]-min_x_wall)/(max_x_wall-min_x_wall)
else:
    print('Preparing geometrical properties')
    t0 = time.perf_counter()
    Wall2Cell=empty((Nwalls,2))
    Wall2Cell[:]=NAN
    nWall2Cell=zeros((Nwalls,1))
    Junction2Wall2Cell=empty((NwallsJun-Nwalls,12))
    Junction2Wall2Cell[:]=NAN
    nJunction2Wall2Cell=zeros((NwallsJun-Nwalls,1))
    Wall2Junction=empty((Nwalls,2))
    Wall2Junction[:]=NAN
    nWall2Junction=zeros((Nwalls,1))
    r_rel=empty((Nwalls,1))
    r_rel[:]=NAN
    x_rel=empty((NwallsJun+Ncells,1))
    x_rel[:]=NAN
    L_diff=((abs(float(Layer_dist[2]-Layer_dist[3]+0.5*(Layer_dist[2]-Layer_dist[4]))*1.0E-04),abs(float(Layer_dist[3] - xyl80_dist)*1.0E-04))) #Diffusion lengths (cm)
    for node, edges in G.adjacency_iter():
        wid=indice[node]
        if wid<Nwalls: #wall that is not a junction (connected to a cell)
            for neighbour, eattr in edges.items(): #Loop on connections (edges)
                cid=int(indice[neighbour])
                if G.node[cid]['type']=='cell':
                    Wall2Cell[wid][int(nWall2Cell[wid])]=cid
                    nWall2Cell[wid]+=1
                elif G.node[neighbour]['type']=='apo': #Node j is a junction
                    Wall2Junction[wid][int(nWall2Junction[wid])]=indice[neighbour]
                    nWall2Junction[wid]+=1
            cid1=Wall2Cell[wid][0]
            cid2=Wall2Cell[wid][1]
            rank1=Cell_rank[int(cid1-NwallsJun)]
            row1=rank2row[int(rank1)]
            if not isnan(cid2):
                rank2=Cell_rank[int(cid2-NwallsJun)]
                row2=rank2row[int(rank2)]
            if not isnan(xyl80_dist):
                if row1 <= rank2row[3] and row2 <= rank2row[3]: #Row on the stelar side
                    if not isnan(cid2):
                        rad_pos=-((Layer_dist2[int(row1)]+Layer_dist2[int(row2)])/2 - xyl80_dist)/(Layer_dist[3] - xyl80_dist) #Radial position of the wall layer relative to the cortical layer closest to tthe xylem vessels (zero) and to the endodermis (1)
                    else:
                        rad_pos=-(Layer_dist2[int(row1)] - xyl80_dist)/(Layer_dist[3] - xyl80_dist) #Radial position of the wall layer relative to the cortical layer closest to the xylem vessels (zero) and to the endodermis (1)
                    r_rel[wid]=max(min(rad_pos,-0.00001),-1)
                else: #Positive r_rel means a wall on the cortical side of the endodermis
                    if not isnan(cid2):
                        rad_pos=(Layer_dist[2]-(Layer_dist2[int(row1)]+Layer_dist2[int(row2)])/2)/(Layer_dist[2]-Layer_dist[3]) #Radial position of the wall layer relative to the cortical layer closest to the epidermis (zero) and to the endodermis (1)
                    else:
                        rad_pos=(Layer_dist[2]-Layer_dist2[int(row1)])/(Layer_dist[2]-Layer_dist[3]) #Radial position of the wall layer relative to the cortical layer closest to the epidermis (zero) and to the endodermis (1)
                    r_rel[wid]=min(max(rad_pos,0.00001),1)
        x_rel[wid]=(position[wid][0]-min_x_wall)/(max_x_wall-min_x_wall)

wallXarea=np.linspace(0,0,NwallsJun)
for i in range(Nwalls):
    wallXarea[i]=thickness*lengths[i] #Cross sectionnal area of the wall (micron²)
    for temp in range(int(nWall2Junction[i])):
        wallXarea[int(Wall2Junction[i][temp])]+=wallXarea[i]/4.0
    wallXarea[i]/=2.0

#At this point, nThickWallPolygonX equals 2 (two thick wall nodes for each wall node)
#We still have to add two more thick junction nodes for each row in nThickWallPolygonX (each row is a "wall polygon")
if Paraview==1 or ParTrack==1:# or Apo_Contagion>0 or Sym_Contagion>0:
    for node, edges in G.adjacency_iter():
        i=indice[node]
        if i<Nwalls: #wall that is not a junction (connected to a cell)
            for neighbour, eattr in edges.items(): #Loop on connections (edges) separated to make sure that Wall2Cell[i] is complete
                j=indice[neighbour]
                if G.node[j]['type']=='apo': #then j is a junction node 
                    for cid in Wall2Cell[i][0:int(nWall2Cell[i])]: #Wall2Cell[i][k] are the cell node ID associated to Wall i
                        if cid not in Junction2Wall2Cell[j-Nwalls]: 
                            Junction2Wall2Cell[j-Nwalls][int(nJunction2Wall2Cell[j-Nwalls])]=cid #Writes the cells indirectly associated to each junction
                            nJunction2Wall2Cell[j-Nwalls]+=1
                            Junction2Wall2[j-Nwalls][int(nJunction2Wall2[j-Nwalls])]=i #Writes the walls directly associated to each junction
                            nJunction2Wall2[j-Nwalls]+=1  #the count in nJunction2Wall is actually the same as nJunction2Wall2Cell 
                        else: #Cell already associated to junction j through another wall => we can consider that the junction j is directly associated to the cell Wall2Cell[i][k] geometrically
                            #Junction j is connected to a cell "cid" from two walls "i" and "wid1"
                            #Find wid1
                            for id1, val in enumerate(Junction2Wall2Cell[j-Nwalls]): 
                                if val==cid: #Finding the position of the current cell in the list of cells related to the junction 
                                    wid1=int(Junction2Wall2[j-Nwalls][id1]) #At the same position in Junction2Wall2 we can find the "other" wall ID (wid1) that was already associated to the same junction and cell
                                    break
                            #Find the thick wall node associated to wid1
                            for id1, val in enumerate(Wall2Cell[wid1]): 
                                if val==cid: #Finding the position of the current cell in the list of cells related to the "other wall"
                                    twpid1=int(Wall2NewWall[wid1][id1]) #At the same position in Wall2NewWall we can find the "new wall" ID that was already associated to the same "other wall" and cell
                                    break
                            #Find the thick wall node associated to i
                            for id1, val in enumerate(Wall2Cell[i]):
                                if val==cid: #Finding the position of the current cell in the list of cells related to the "current wall"
                                    twpid2=int(Wall2NewWall[i][id1]) #At the same position in Wall2NewWall we can find the "new wall" ID that was already associated to the same "current wall" and cell
                                    break
                            #Calculating the slope of the line connecting the junction j to the cell
                            if not position[cid][0]==position[j][0]:
                                slopeCG=(position[cid][1]-position[j][1])/(position[cid][0]-position[j][0]) #slope of the line connecting the junction node j to the center of gravity of the cell Wall2Cell[i][k]
                            else:
                                slopeCG=inf
                            #Calculating the position of the "thick junction node" on the side of cell cid
                            x=position[j][0]+cos(arctan(slopeCG))*thicknessJunction_disp/2*sign(position[cid][0]-position[j][0])
                            y=position[j][1]+sin(arctan(slopeCG))*thicknessJunction_disp/2*sign(position[cid][0]-position[j][0])
                            ThickWalls.append(array((twpid,j,int(cid),x,y,twpid1,twpid2,Borderlink[j]))) #Adds the thick wall node ID, its parent node ID, the associated cell ID, the new coordinates in X and Y, and in this case, the two neighbouring walls that will consitute 2 neighbouring "cells" in the sense of pvtk     G.node[j]['borderlink']
                            ThickWalls[twpid1][int(5+nThickWalls[twpid1])]=twpid
                            ThickWalls[twpid2][int(5+nThickWalls[twpid2])]=twpid
                            nThickWalls[twpid1]+=1
                            nThickWalls[twpid2]+=1
                            Cell2ThickWalls[int(cid-NwallsJun)][int(nCell2ThickWalls[int(cid-NwallsJun)])]=twpid
                            nCell2ThickWalls[int(cid-NwallsJun)]+=1
                            ThickWallsX.append((twpidX,x,y,j,cid,i,wid1)) #Adds the thick wall node ID, the new coordinates in X and Y, its parent (junction) node ID, the associated cell ID, the two original neighbouring walls
                            #Which of the 2 polygons associated to wall node "i" do we add the thick junction node to?
                            #Only one of the polygons is related to junction node "j". Need to look at position into Wall2Junction
                            for id1, val in enumerate(Wall2Junction[i]): 
                                if val==j: #Finding the position "id1" (index equal to 0 or 1) of the current junction in the list of junctions related to the wall "i" 
                                    break #2*i+id1 is the row of the polygon connecting wall "i" to junction "j", because apoplastic connections to wall "i" come in the same order when creating the matrices Wall2Junction and ThickWallsX
                            if ThickWallsX[int(ThickWallPolygonX[int(2*i+id1)][1])][4]==cid: #The polygon nodes need to be consecutive ("no zigzag"), hence the ordering of the 3rd and 4th twpid matter (the 4th twpid should be on the same side as the 1st)
                                ThickWallPolygonX[int(2*i+id1)][2]=twpidX
                            else:
                                ThickWallPolygonX[int(2*i+id1)][3]=twpidX
                            nThickWallPolygonX[int(2*i+id1)]+=1
                            #We do the same for the other polygon that includes twpid (connects wall "wid1" to node "j")
                            for id1, val in enumerate(Wall2Junction[wid1]): 
                                if val==j: #Finding the position "id1" (index equal to 0 or 1) of the current junction in the list of junctions related to the wall "i" 
                                    break #2*i+id1 is the row of the polygon connecting wall "i" to junction "j"
                            if ThickWallsX[int(ThickWallPolygonX[int(2*wid1+id1)][1])][4]==cid: #The polygon nodes need to be consecutive ("no zigzag"), hence the ordering of the 3rd and 4th twpid matter
                                ThickWallPolygonX[int(2*wid1+id1)][2]=twpidX
                            else:
                                ThickWallPolygonX[int(2*wid1+id1)][3]=twpidX
                            nThickWallPolygonX[int(2*wid1+id1)]+=1
                            #New walls ID for each original wall
                            Wall2NewWallX[j][int(nWall2NewWallX[j])]=twpidX
                            nWall2NewWallX[j]+=1
                            #Increase the wall count
                            twpid+=1
                            twpidX+=1
                            if Borderlink[j]==1: #G.node[j]['borderlink']==1: #If the wall is at a border, there is only one cell => need to add the opposite new wall node independently
                                if not slopeCG==0.0:
                                    slopeCG=-1/slopeCG
                                else:
                                    slopeCG=inf
                                x=position[j][0]-cos(arctan(slopeCG))*thickness_disp/2*sign(position[cid][0]-position[j][0]) #-> new position of the wall node opposite to the current cell side
                                y=position[j][1]-sin(arctan(slopeCG))*thickness_disp/2*sign(position[cid][0]-position[j][0])
                                ThickWallsX.append((twpidX,x,y,j,inf,i))
                                #This new wall node is added to the border polygon. Is the border polygon associated to wall node "i" or "wid1"?
                                if Borderlink[i]==1: #G.node[i]['borderlink']==1:
                                    wid=i
                                else:
                                    wid=wid1
                                #Which of the 2 polygons associated to wall node "wid" do we add the twpid to?
                                #Only one of them is related to junction node "j". Need to look at position into Wall2Junction
                                for id1, val in enumerate(Wall2Junction[wid]): 
                                    if val==j: #Finding the position "id1" (index equal to 0 or 1) of the current junction in the list of junctions related to the wall "i" 
                                        break #2*i+id1 is the row of the polygon connecting wall "i" to junction "j"
                                if ThickWallsX[int(ThickWallPolygonX[int(2*wid+id1)][1])][4]==cid: #The polygon nodes need to be consecutive ("no zigzag"), hence the ordering of the 3rd and 4th twpid matter
                                    ThickWallPolygonX[int(2*wid+id1)][3]=twpidX
                                else:
                                    ThickWallPolygonX[int(2*wid+id1)][2]=twpidX
                                nThickWallPolygonX[int(2*wid+id1)]+=1
                                #New walls ID for each original wall
                                Wall2NewWallX[j][int(nWall2NewWallX[j])]=twpidX
                                nWall2NewWallX[j]+=1
                                twpidX+=1
else:
    for node, edges in G.adjacency_iter():
        i=indice[node]
        if i<Nwalls: #wall that is not a junction (connected to a cell)
            for neighbour, eattr in edges.items(): #Loop on connections (edges) separated to make sure that Wall2Cell[i] is complete
                j=indice[neighbour]
                if G.node[j]['type']=='apo': #then j is a junction node 
                    for cid in Wall2Cell[i][0:int(nWall2Cell[i])]: #Wall2Cell[i][k] are the cell node ID associated to Wall i
                        if cid not in Junction2Wall2Cell[j-Nwalls]: 
                            Junction2Wall2Cell[j-Nwalls][int(nJunction2Wall2Cell[j-Nwalls])]=cid #Writes the cells indirectly associated to each junction
                            nJunction2Wall2Cell[j-Nwalls]+=1

t1 = time.perf_counter()
print(t1-t0, "seconds process time")

#######################################
##Variables for potential calculation##
#######################################

#Import Hydraulic data
print('Importing hydraulic data')
t0 = time.perf_counter()
kwrange=etree.parse(dir + Project + inputs + Hydr).getroot().xpath('kwrange/kw')
kw_barrier_range=etree.parse(dir + Project + inputs + Hydr).getroot().xpath('kw_barrier_range/kw_barrier')
kmb=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('km')[0].get("value"))
kAQPrange=etree.parse(dir + Project + inputs + Hydr).getroot().xpath('kAQPrange/kAQP')
Kplrange=etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Kplrange/Kpl')  #cm^3/hPa/d/plasmodesmata
PD_freq_source=int(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Freq_source')[0].get("value"))
PD_freq_factor=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Freq_factor')[0].get("value"))
if PD_freq_source==2:
    newpath_geom=False
    Wall_PD=0.0*ones((Nwalls,1))
    for w in Walls_PD: #Loop on walls, by cell - wall association, hence a wall can be repeated if associated to two cells
        wid= int(w.get("id")) #Wall id number
        try:
            Wall_PD[wid]=float(w.get("PD_freq"))*PD_freq_factor  #Number of plasmodesmata per micron² on this wall
        except:
            if wid>0:
                print('Missing PD frequency data at wid =',wid)
                error
            else:
                newpath_geom=True
                break
    if newpath_geom:
        Wall_PD = np.loadtxt(fname = dir + 'cellsetdata/' + path_geom[0:-4] + '_PDfreq.txt')*PD_freq_factor/100*10 #Contains    freq (number of PD per 100 micron wall length) for a cross section of 0.1 micron thickness -> number of PD per micron² of wall
        newtree = ET.parse(dir + 'cellsetdata/' + path_geom)
        newroot = newtree.getroot()
        for wall in newroot[2].iter('wall'):
            wid=int(wall.get("id"))
            #slices=ET.SubElement(wall, 'slices')
            #z_prox=height #position of the proximal end of the root segment (initialization)
            #for iLayer in range(AxialLayers):
            #    slice_iL=ET.SubElement(slices, 'slice')
            #    slice_iL.set('z_prox', str(z_prox))
            wall.set('PD_freq', str(float(Wall_PD[wid]))) #freq (number of PD per micron² of wall)
            #    z_prox+=height
        path_geom=path_geom[0:-4]+'_PDfreq.xml'
        newtree.write(dir + 'cellsetdata/' + path_geom)
Fplxheight=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Fplxheight')[0].get("value"))*PD_freq_factor #Product of plasmodesmatal frequency (per square cm) by cell axial length when the frequency was counted (cm) => units: per cm of cell membrane perimeter, to be multiplied by cell perimeter in cm to obtain a quantity of plasmodesmata
Fplxheight_epi_exo=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Fplxheight_epi_exo')[0].get("value"))*PD_freq_factor #Because they are "cell axial length independent", these values conserve the quantity of plasmodesmata when cells elongate
Fplxheight_outer_cortex=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Fplxheight_outer_cortex')[0].get("value"))*PD_freq_factor
Fplxheight_cortex_cortex=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Fplxheight_cortex_cortex')[0].get("value"))*PD_freq_factor
Fplxheight_cortex_endo=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Fplxheight_cortex_endo')[0].get("value"))*PD_freq_factor
Fplxheight_endo_endo=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Fplxheight_endo_endo')[0].get("value"))*PD_freq_factor
Fplxheight_endo_peri=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Fplxheight_endo_peri')[0].get("value"))*PD_freq_factor
Fplxheight_peri_peri=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Fplxheight_peri_peri')[0].get("value"))*PD_freq_factor
Fplxheight_peri_stele=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Fplxheight_peri_stele')[0].get("value"))*PD_freq_factor
Fplxheight_stele_stele=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Fplxheight_stele_stele')[0].get("value"))*PD_freq_factor
Fplxheight_stele_comp=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Fplxheight_stele_comp')[0].get("value"))*PD_freq_factor
Fplxheight_peri_comp=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Fplxheight_peri_comp')[0].get("value"))*PD_freq_factor
Fplxheight_comp_comp=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Fplxheight_comp_comp')[0].get("value"))*PD_freq_factor
Fplxheight_comp_sieve=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Fplxheight_comp_sieve')[0].get("value"))*PD_freq_factor
Fplxheight_peri_sieve=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Fplxheight_peri_sieve')[0].get("value"))*PD_freq_factor
Fplxheight_stele_sieve=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Fplxheight_stele_sieve')[0].get("value"))*PD_freq_factor
Defheight=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Defheight')[0].get("value"))
Kax_source=int(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Kax_source')[0].get("value")) #Source of axial conductance values (1=Poiseuille;2=prescribed in Hydraulics input file)
#K_sieve=float(etree.parse(dir + Project + inputs + Hydr).getroot().xpath('K_sieve')[0].get("value")) #Sieve tube hydraulic conductance
if Kax_source==2:
    K_xyl_range=etree.parse(dir + Project + inputs + Hydr).getroot().xpath('K_xyl_range/K_xyl') #Xylem vessel axial hydraulic conductance
    K_sieve_range=etree.parse(dir + Project + inputs + Hydr).getroot().xpath('K_sieve_range/K_sieve') #Phloem sieve tubes axial hydraulic conductance
Xcontactrange=etree.parse(dir + Project + inputs + Hydr).getroot().xpath('Xcontactrange/Xcontact')
path_hydraulics=etree.parse(dir + Project + inputs + Hydr).getroot().xpath('path_hydraulics/Output')
t1 = time.perf_counter()
print(t1-t0, "seconds process time")

#Import boundary conditions
print('Importing boundary conditions')
t0 = time.perf_counter()
Psi_soil_range=etree.parse(dir + Project + inputs + BC).getroot().xpath('Psi_soil_range/Psi_soil')
BC_xyl_range=etree.parse(dir + Project + inputs + BC).getroot().xpath('BC_xyl_range/BC_xyl')
BC_sieve_range=etree.parse(dir + Project + inputs + BC).getroot().xpath('BC_sieve_range/BC_sieve')
Psi_cell_range=etree.parse(dir + Project + inputs + BC).getroot().xpath('Psi_cell_range/Psi_cell')
Elong_cell_range=etree.parse(dir + Project + inputs + BC).getroot().xpath('Elong_cell_range/Elong_cell')
#Elong_cell_kappa=float(etree.parse(dir + Project + inputs + BC).getroot().xpath('Elong_cell')[0].get("kappa_dot")) #Rate change of curvature (radian/micron/d)
Water_fraction_apo=float(etree.parse(dir + Project + inputs + BC).getroot().xpath('Water_fractions')[0].get("Apoplast")) #Relative volumetric fraction of water in the apoplast (dimensionless)
Water_fraction_sym=float(etree.parse(dir + Project + inputs + BC).getroot().xpath('Water_fractions')[0].get("Symplast")) #Relative volumetric fraction of water in the symplast (dimensionless)
path_BC=etree.parse(dir + Project + inputs + BC).getroot().xpath('path_scenarios/Output')[0].get("path")
Nhydraulics=len(path_hydraulics) #Total number of hydraulic parameters sets
Nkw=len(kwrange)
NKpl=len(Kplrange)
NkAQP=len(kAQPrange)
NXcontact=len(Xcontactrange)
Nkw_barrier=len(kw_barrier_range)
Nscenarios=len(Psi_soil_range) #Including the "forced scenario"
Nmaturity=len(Maturityrange)
if PileUp==2:
    Nmaturity_loop=1
else:
    Nmaturity_loop=Nmaturity
Nprint_layer=len(Printrange)
Psi_soil=zeros((2,Nscenarios))
Os_soil=zeros((6,Nscenarios))
Os_xyl=zeros((6,Nscenarios))
Ana_Os=False #Do we calculate solute stationnary fluxes analytical solution?
for count in range(1,Nscenarios):
    Psi_soil[0][count]=float(Psi_soil_range[count].get("pressure_left")) #Soil pressure potential (hPa)
    Psi_soil[1][count]=float(Psi_soil_range[count].get("pressure_right"))
    Os_soil[0][count]=float(Psi_soil_range[count].get("osmotic_left"))
    Os_soil[1][count]=float(Psi_soil_range[count].get("osmotic_right"))
    Os_soil[2][count]=float(Psi_soil_range[count].get("osmotic_symmetry"))
    Os_soil[3][count]=float(Psi_soil_range[count].get("osmotic_shape")) #1 for linear, >1 for outer slope flat, <1 for inner slope flat
    Os_soil[4][count]=float(etree.parse(dir + Project + inputs + BC).getroot().xpath('Psi_soil_range/osmotic_diffusivity')[0].get("value"))
    #Os_soil[5][count]=float(etree.parse(dir + Project + inputs + BC).getroot().xpath('Psi_soil_range/osmotic_convection')[0].get("flag"))
    Os_xyl[0][count]=float(BC_xyl_range[count].get("osmotic_xyl"))
    Os_xyl[1][count]=float(BC_xyl_range[count].get("osmotic_endo"))
    Os_xyl[2][count]=float(BC_xyl_range[count].get("osmotic_symmetry"))
    Os_xyl[3][count]=float(BC_xyl_range[count].get("osmotic_shape")) #The symmetry is central. 1 for linear, >1 for outer slope flat, <1 for inner slope flat
    Os_xyl[4][count]=float(etree.parse(dir + Project + inputs + BC).getroot().xpath('BC_xyl_range/osmotic_diffusivity')[0].get("value"))
    #Os_xyl[5][count]=float(etree.parse(dir + Project + inputs + BC).getroot().xpath('BC_xyl_range/osmotic_convection')[0].get("flag"))
    if not Os_xyl[4][count]==0 and not Os_soil[4][count]==0:
        Ana_Os=True
        print('Calculation of analytical solution for radial solute transport in cell walls in scenario ',count)

#Get solute BC conditions
Transient_tracer=False
Dif_BC=False
if Sym_Contagion==2 and Apo_Contagion==2:
    N_Dirichlet_ini=[]
    cc_Dirichlet_ini=[]
    N_AxialNeumann_ini=[]
    cc_AxialNeumann_ini=[]
    N_Prod_ini=[]
    cc_Prod_ini=[]
    if int(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Time_information/Transient_sim')[0].get("value"))==1:
        Transient_tracer=True
    if int(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Hormone_movement/Dif_BC')[0].get("flag"))==1:
        Dif_BC=True
    for source in Sym_source_ini_range:
        cellnumber1=int(source.get("id"))
        if cellnumber1>-1 and cellnumber1<Ncells:
            if float(source.get("flow_BC"))==0.0:
                for iLayer in range(int(source.get("layer_range_min")),int(source.get("layer_range_max"))+1):
                    if iLayer*Ntot+cellnumber1+NwallsJun not in N_Dirichlet_ini:
                        N_Dirichlet_ini.append(iLayer*Ntot+cellnumber1+NwallsJun)
                        cc=float(source.get("cc_dist"))
                        if isnan(cc):
                            error('Dirichlet solute BC is set (flow_BC=0) but no concentration is prescribed in cc_dist')
                        else:
                            cc_Dirichlet_ini.append(cc)
            elif float(source.get("flow_BC"))==1.0:
                if cellnumber1 not in N_AxialNeumann_ini:
                    N_AxialNeumann_ini.append(cellnumber1)
                    cc_AxialNeumann_ini.append(array((float(source.get("cc_dist")),float(source.get("cc_prox")))))
            elif float(source.get("flow_BC"))==2.0:
                if int(source.get("layer_range_min"))>int(source.get("layer_range_max")):
                    error('layer_range_min must be smaller than layer_range_max')
                else:
                    for iLayer in range(int(source.get("layer_range_min")),int(source.get("layer_range_max"))+1):
                        cid=iLayer*Ntot+cellnumber1+NwallsJun
                        if cid not in N_Prod_ini:
                            N_Prod_ini.append(cid)
                            cc=float(source.get("cc_prod"))
                            if isnan(cc):
                                error('Solute production BC is set (flow_BC=2) but no concentration is prescribed in cc_prod')
                            else:
                                cc_Prod_ini.append(cc)
    if Transient_tracer==True:
        N_Dirichlet_transi=[]
        cc_Dirichlet_transi=[]
        N_AxialNeumann_transi=[]
        cc_AxialNeumann_transi=[]
        N_Prod_transi=[]
        cc_Prod_transi=[]
        for source in Sym_source_transi_range:
            cellnumber1=int(source.get("id"))
            if cellnumber1>-1 and cellnumber1<Ncells:
                if float(source.get("flow_BC"))==0.0:
                    if int(source.get("layer_range_min"))>int(source.get("layer_range_max")):
                        error('layer_range_min must be smaller than layer_range_max')
                    else:
                        for iLayer in range(int(source.get("layer_range_min")),int(source.get("layer_range_max"))+1):
                            if iLayer*Ntot+cellnumber1+NwallsJun not in N_Dirichlet_transi:
                                N_Dirichlet_transi.append(iLayer*Ntot+cellnumber1+NwallsJun)
                                cc=float(source.get("cc_dist"))
                                if isnan(cc):
                                    error('Dirichlet solute BC is set (flow_BC=0) but no concentration is prescribed in cc_dist')
                                else:
                                    cc_Dirichlet_transi.append(cc)
                elif float(source.get("flow_BC"))==1.0:
                    if cellnumber1 not in N_AxialNeumann_transi:
                        N_AxialNeumann_transi.append(cellnumber1)
                        cc_AxialNeumann_transi.append(array((float(source.get("cc_dist")),float(source.get("cc_prox")))))
                elif float(source.get("flow_BC"))==2.0:
                    if int(source.get("layer_range_min"))>int(source.get("layer_range_max")):
                        error('layer_range_min must be smaller than layer_range_max')
                    else:
                        for iLayer in range(int(source.get("layer_range_min")),int(source.get("layer_range_max"))+1):
                            cid=iLayer*Ntot+cellnumber1+NwallsJun
                            if cid not in N_Prod_transi:
                                N_Prod_transi.append(cid)
                                cc=float(source.get("cc_prod"))
                                if isnan(cc):
                                    error('Solute production BC is set (flow_BC=2) but no concentration is prescribed in cc_prod')
                                else:
                                    cc_Prod_transi.append(cc)
        dt_tracer=float(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Time_information/Time_step_transi')[0].get("value")) #Time step for transient tracer simulation (day)
        dt_adapt=False
        if int(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Time_information/Time_step_adapt')[0].get("value"))==1:
            if PileUp==1:
                error('Adaptive time steps not compatible with PileUp==1 option (see Hormone_carriers* and Geometry* files respectively')
            dt_adapt=True
            dcc_low=float(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Time_information/Time_step_adapt')[0].get("cc_change_low"))
            dcc_high=float(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Time_information/Time_step_adapt')[0].get("cc_change_high"))
            dt_factor=float(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Time_information/Time_step_adapt')[0].get("time_step_factor"))
            dt_min=float(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Time_information/Time_step_adapt')[0].get("time_step_min"))
            dt_max=float(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Time_information/Time_step_adapt')[0].get("time_step_max"))
        switch_dur=float(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Time_information/Switch_duration')[0].get("value")) #Transition time during which the boundary condition progressively changes from ini to transi values
        switch_nt=float(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Time_information/Switch_fractionning')[0].get("value")) #Fractionning of the transition time into small time steps of duration "Switch_duration/Switch_fractionning"
        Freq_obs=float(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Time_information/Observation_frequency')[0].get("value"))
        Print_times=etree.parse(dir + Project + inputs + Horm).getroot().xpath('Time_information/Print_times/Print_time') #Print times of transient tracer simulation outputs (day)
        Obs_node_range=etree.parse(dir + Project + inputs + Horm).getroot().xpath('Observation_range/Obs') #Cell id where the tracer is observed
        if int(Obs_node_range[0].get("id"))==-2: #All cells of the bottom layer
            Obs_nodes=[]
            for cid in range(Ncells):
                Obs_nodes.append(int(cid)+NwallsJun) #+int(node.get("layer"))*Ntot
        elif int(Obs_node_range[0].get("id"))==-1:
            Obs_nodes=[]
        else:
            Obs_nodes=[]
            for node in Obs_node_range:
                Obs_nodes.append(int(node.get("id"))+NwallsJun+int(node.get("layer"))*Ntot)
        Ini_cond=int(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Initial_conditions/Type')[0].get("value"))
        if Ini_cond==2:
            Ini_path=etree.parse(dir + Project + inputs + Horm).getroot().xpath('Initial_conditions/Type')[0].get("path")
            tree_C = etree.parse(dir + Project + Ini_path)
            rootelt_C = tree_C.getroot()
            t_resume=float(rootelt_C[0][2].get('value'))
            AxialLayers=int(rootelt_C[0][3].get('value'))
            while float(Print_times[0].get('value'))<=t_resume:
                print('Remove',float(Print_times[0].get('value')))
                Print_times.remove(Print_times[0])
        t_end=float(Print_times[-1].get('value'))
        Nprint=len(Print_times)
        Num_meth=int(etree.parse(dir + Project + inputs + Horm).getroot().xpath('Numerical_method/Method')[0].get("value"))
    N_SoilNeumann_ini=[]
    cc_SoilNeumann_ini=[]
    for source in Apo_source_ini_range:
        cellnumber1=int(source.get("id"))
        tissue=int(source.get("tissue"))
        if cellnumber1>-1 and cellnumber1<Ncells:
            for r in Cell2Wall_loop[cellnumber1]:
                for iLayer in range(int(source.get("layer_range_min")),int(source.get("layer_range_max"))+1):
                    if float(source.get("flow_BC"))==0.0:
                        if iLayer*Ntot+int(r.get("id")) not in N_Dirichlet_ini:
                            N_Dirichlet_ini.append(iLayer*Ntot+int(r.get("id")))
                            cc=float(source.get("cc_dist"))
                            if isnan(cc):
                                error('Dirichlet solute BC is set (flow_BC=0) but no concentration is prescribed in cc_dist')
                            else:
                                cc_Dirichlet_ini.append(cc)
                    if float(source.get("flow_BC"))==1.0:
                        wid=int(r.get("id"))
                        if Borderwall[wid]==1:
                            cc_soil=float(source.get("cc_soil"))
                            if Ntot*iLayer+wid not in N_SoilNeumann_ini:
                                N_SoilNeumann_ini.append(Ntot*iLayer+wid)
                                cc_SoilNeumann_ini.append(cc_soil)
                            jid1=Ntot*iLayer+Wall2Junction[wid][0]
                            if jid1 not in N_SoilNeumann_ini:
                                N_SoilNeumann_ini.append(jid1)
                                cc_SoilNeumann_ini.append(cc_soil)
                            jid2=Ntot*iLayer+Wall2Junction[wid][1]
                            if jid2 not in N_SoilNeumann_ini:
                                N_SoilNeumann_ini.append(jid2)
                                cc_SoilNeumann_ini.append(cc_soil)
        if tissue==2 and float(source.get("flow_BC"))==1.0:
            cc_soil=float(source.get("cc_soil"))
            for iLayer in range(int(source.get("layer_range_min")),int(source.get("layer_range_max"))+1):
                for wid in Borderwall:
                    if Ntot*iLayer+wid not in N_SoilNeumann_ini:
                        N_SoilNeumann_ini.append(Ntot*iLayer+wid)
                        cc_SoilNeumann_ini.append(cc_soil)
                for jid in Borderjunction:
                    if Ntot*iLayer+jid not in N_SoilNeumann_ini:
                        N_SoilNeumann_ini.append(Ntot*iLayer+jid)
                        cc_SoilNeumann_ini.append(cc_soil)
        elif tissue==2 and float(source.get("flow_BC"))==0.0 and not isnan(float(source.get("cc_soil"))): #Dirichlet at root surface only
            cc_soil=float(source.get("cc_soil"))
            for iLayer in range(int(source.get("layer_range_min")),int(source.get("layer_range_max"))+1):
                for wid in Borderwall:
                    if Ntot*iLayer+wid not in N_Dirichlet_ini:
                        N_Dirichlet_ini.append(Ntot*iLayer+wid)
                        cc_Dirichlet_ini.append(cc_soil)
                for jid in Borderjunction:
                    if Ntot*iLayer+jid not in N_Dirichlet_ini:
                        N_Dirichlet_ini.append(Ntot*iLayer+jid)
                        cc_Dirichlet_ini.append(cc_soil)
        elif tissue>-1:
            for cid in range(NwallsJun,Ntot):
                if tissue==G.node[cid]['cgroup']:
                    for r in Cell2Wall_loop[cid-NwallsJun]:
                        for iLayer in range(int(source.get("layer_range_min")),int(source.get("layer_range_max"))+1):
                            if float(source.get("flow_BC"))==0.0:
                                if iLayer*Ntot+int(r.get("id")) not in N_Dirichlet_ini:
                                    N_Dirichlet_ini.append(iLayer*Ntot+int(r.get("id")))
                                    cc=float(source.get("cc_dist"))
                                    if isnan(cc):
                                        error('Dirichlet solute BC is set (flow_BC=0) but no concentration is prescribed in cc_dist')
                                    else:
                                        cc_Dirichlet_ini.append(cc)
    if Transient_tracer==True:
        N_SoilNeumann_transi=[]
        cc_SoilNeumann_transi=[]
        for source in Apo_source_transi_range:
            cellnumber1=int(source.get("id"))
            tissue=int(source.get("tissue"))
            if cellnumber1>-1 and cellnumber1<Ncells:
                for iLayer in range(int(source.get("layer_range_min")),int(source.get("layer_range_max"))+1):
                    if float(source.get("flow_BC"))==0.0:
                        for r in Cell2Wall_loop[cellnumber1]:
                            if iLayer*Ntot+int(r.get("id")) not in N_Dirichlet_transi:
                                N_Dirichlet_transi.append(iLayer*Ntot+int(r.get("id")))
                                cc=float(source.get("cc_dist"))
                                if isnan(cc):
                                    error('Dirichlet solute BC is set (flow_BC=0) but no concentration is prescribed in cc_dist')
                                else:
                                    cc_Dirichlet_transi.append(cc)
                    else:
                        for r in Cell2Wall_loop[cellnumber1]:
                            wid=int(r.get("id"))
                            if Borderlink[wid]==1:
                                cc_soil=float(source.get("cc_soil"))
                                if Ntot*iLayer+wid not in N_SoilNeumann_transi:
                                    N_SoilNeumann_transi.append(iLayer*Ntot+wid)
                                    cc_SoilNeumann_transi.append(cc_soil)
                                jid1=Ntot*iLayer+Wall2Junction[wid][0]
                                if jid1 not in N_SoilNeumann_transi:
                                    N_SoilNeumann_transi.append(jid1)
                                    cc_SoilNeumann_transi.append(cc_soil)
                                jid2=Ntot*iLayer+Wall2Junction[wid][1]
                                if jid2 not in N_SoilNeumann_transi:
                                    N_SoilNeumann_transi.append(jid2)
                                    cc_SoilNeumann_transi.append(cc_soil)
            if tissue==2 and float(source.get("flow_BC"))==1.0:
                cc_soil=float(source.get("cc_soil"))
                for iLayer in range(int(source.get("layer_range_min")),int(source.get("layer_range_max"))+1):
                    for wid in Borderwall:
                        if Ntot*iLayer+wid not in N_SoilNeumann_transi:
                            N_SoilNeumann_transi.append(Ntot*iLayer+wid)
                            cc_SoilNeumann_transi.append(cc_soil)
                    for jid in Borderjunction:
                        if Ntot*iLayer+jid not in N_SoilNeumann_transi:
                            N_SoilNeumann_transi.append(Ntot*iLayer+jid)
                            cc_SoilNeumann_transi.append(cc_soil)
            elif tissue==2 and float(source.get("flow_BC"))==0.0 and not isnan(float(source.get("cc_soil"))):
                cc_soil=float(source.get("cc_soil"))
                for iLayer in range(int(source.get("layer_range_min")),int(source.get("layer_range_max"))+1):
                    for wid in Borderwall:
                        if Ntot*iLayer+wid not in N_Dirichlet_transi:
                            N_Dirichlet_transi.append(Ntot*iLayer+wid)
                            cc_Dirichlet_transi.append(cc_soil)
                    for jid in Borderjunction:
                        if Ntot*iLayer+jid not in N_Dirichlet_transi:
                            N_Dirichlet_transi.append(Ntot*iLayer+jid)
                            cc_Dirichlet_transi.append(cc_soil)
            elif tissue>-1:
                for cid in range(NwallsJun,Ntot):
                    if tissue==G.node[cid]['cgroup']:
                        for r in Cell2Wall_loop[cid-NwallsJun]:
                            for iLayer in range(int(source.get("layer_range_min")),int(source.get("layer_range_max"))+1):
                                if float(source.get("flow_BC"))==0.0:
                                    if iLayer*Ntot+int(r.get("id")) not in N_Dirichlet_transi:
                                        N_Dirichlet_transi.append(iLayer*Ntot+int(r.get("id")))
                                        cc=float(source.get("cc_dist"))
                                        if isnan(cc):
                                            error('Dirichlet solute BC is set (flow_BC=0) but no concentration is prescribed in cc_dist')
                                        else:
                                            cc_Dirichlet_transi.append(cc)
    if PileUp==1:
        Sym_cc_flows_ini=zeros((Nhydraulics,Nmaturity,Nscenarios,len(N_AxialNeumann_ini),2))
        for ih in range(Nhydraulics):
            for iMaturity in range(Nmaturity):
                for i in range(len(N_AxialNeumann_ini)):
                    Sym_cc_flows_ini[ih][iMaturity][0][i][:]=cc_AxialNeumann_ini[i][:]
        if size(cc_AxialNeumann_transi)>0:
            Sym_cc_flows_transi=zeros((Nhydraulics,Nmaturity,Nscenarios,len(N_AxialNeumann_transi),2,Nprint+1))
            for ih in range(Nhydraulics):
                for iMaturity in range(Nmaturity):
                    for it in range(Nprint+1):
                        for i in range(len(N_AxialNeumann_transi)):
                            if not cc_AxialNeumann_transi[i][1] == 0.0:
                                error('The current version of the code cannot deal with the proximal solute flux BC in PileUp mode yet')
                            Sym_cc_flows_transi[ih][iMaturity][0][i][:][it]=cc_AxialNeumann_transi[i][:] #The first slice currently has the same BC at all times after time 0 (though the BC differs before time 0, see Sym_cc_flows_ini)
    if Transient_tracer:
        if Ini_cond==2:
            #Read initial condition from xml file
            soln_C = empty((Ntot*AxialLayers))
            for wall in rootelt_C.xpath('walls/wall/slices'):
                wid=int(wall.getparent().get('id'))
                for iLayer in range(AxialLayers):
                    soln_C[iLayer*Ntot+wid]=float(wall[iLayer].get('cc1'))
                    #print(iLayer*Ntot+wid,float(wall[iLayer].get('cc1')))
            for cell in rootelt_C.xpath('cells/cell/slices'):
                cid=int(cell.getparent().get('id'))+NwallsJun
                for iLayer in range(AxialLayers):
                    soln_C[iLayer*Ntot+cid]=float(cell[iLayer].get('cc1'))
            for junction in rootelt_C.xpath('junctions/junction/slices'):
                jid=int(junction.getparent().get('id'))
                for iLayer in range(AxialLayers):
                    soln_C[iLayer*Ntot+jid]=float(junction[iLayer].get('cc1'))
t1 = time.perf_counter()
print(t1-t0, "seconds process time")

#Unit changes
sperd=24.0*3600.0 #(seconds per day)
cmperm=100.0 #(cm per metre)

#Start the loop of hydraulic properties
print('Hydraulic properties loop')
t0 = time.perf_counter()
for h in range(Nhydraulics):
    print('   ')
    print('Hydraulic network #'+str(h))
    newpath=dir+Project+Output_path+Plant+'/'+path_BC+path_hydraulics[h].get("path")
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    
    text_file = open(newpath+"position.txt", "w")
    with open(newpath+"position.txt", "a") as myfile:
        for i in range(Ntot):
            myfile.write(str(position[i][:])+" \n")
    myfile.close()
    text_file.close()
    
    #System solving
    if PileUp==2:
        Psi_xyl=empty((2,1,Nscenarios))
        dPsi_xyl=empty((1,Nscenarios))
        Psi_sieve=empty((2,1,Nscenarios))
        dPsi_sieve=empty((1,Nscenarios))
        Q_tot=zeros((1,Nscenarios)) #(cm^3/d) Total flow rate at root surface
        kr_tot=zeros((1,1))
        OsCellLayer=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
        nOsCellLayer=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
        OsWallLayer=zeros((int(r_discret[0]),1,Nscenarios))
        nOsWallLayer=zeros((int(r_discret[0]),1,Nscenarios)) #Used for averaging OsWallLayer
        NWallLayer=zeros((int(r_discret[0]),1,Nscenarios))
        Q_xyl_layer=zeros((int(r_discret[0]),1,Nscenarios))
        Q_sieve_layer=zeros((int(r_discret[0]),1,Nscenarios))
        Q_elong_layer=zeros((int(r_discret[0]),1,Nscenarios))
        Hydropatterning=empty((1,Nscenarios))
        Hydrotropism=empty((1,Nscenarios))
    else:
        Psi_xyl=empty((2,Nmaturity,Nscenarios))
        dPsi_xyl=empty((Nmaturity,Nscenarios))
        Psi_sieve=empty((2,Nmaturity,Nscenarios))
        dPsi_sieve=empty((Nmaturity,Nscenarios))
        Q_tot=zeros((Nmaturity,Nscenarios)) #(cm^3/d) Total flow rate at root surface
        kr_tot=zeros((Nmaturity,1))
        OsCellLayer=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
        nOsCellLayer=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
        OsWallLayer=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
        nOsWallLayer=zeros((int(r_discret[0]),Nmaturity,Nscenarios)) #Used for averaging OsWallLayer
        NWallLayer=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
        Q_xyl_layer=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
        Q_sieve_layer=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
        Q_elong_layer=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
        Hydropatterning=empty((Nmaturity,Nscenarios))
        Hydrotropism=empty((Nmaturity,Nscenarios))
    Psi_xyl[:]=NAN
    dPsi_xyl[:]=NAN
    iEquil_xyl=NAN #index of the equilibrium root xylem pressure scenario
    Flow_xyl=empty((2,Nxyl+1,Nscenarios))
    Flow_xyl[:]=NAN
    Psi_sieve[:]=NAN
    dPsi_sieve[:]=NAN
    iEquil_sieve=NAN #index of the equilibrium root phloem pressure scenario
    Flow_sieve=empty((2,Nsieve+1,Nscenarios))
    Flow_sieve[:]=NAN
    Os_sieve=zeros((1,Nscenarios))
    Os_cortex=zeros((1,Nscenarios))
    Os_hetero=zeros((1,Nscenarios))
    s_factor=zeros((1,Nscenarios))
    s_hetero=zeros((1,Nscenarios))
    Elong_cell=zeros((1,Nscenarios))
    Elong_cell_side_diff=zeros((1,Nscenarios))
    TotLayers=0
    iM=0
    Nlayers=empty((Nmaturity))
    for Maturity in Maturityrange:
        Nlayers[iM]=int(Maturity.get("Nlayers"))
        iM+=1
    TotLayers=int(sum(Nlayers))
    Os_apo_eq=zeros((Nmaturity,2,Nscenarios))
    Os_sym_eq=zeros((Nmaturity,2,Nscenarios))
    if PileUp==2 or Nlayers[-1]>1:
        Os_apo_eq_lay=zeros((TotLayers,2,Nscenarios))
        Os_sym_eq_lay=zeros((TotLayers,2,Nscenarios))
    UptakeLayer_plus=zeros((int(r_discret[0]),TotLayers,Nscenarios))
    UptakeLayer_minus=zeros((int(r_discret[0]),TotLayers,Nscenarios))
    UptakeLayer_plus_InterCid=zeros((int(r_discret[0]),TotLayers,Nscenarios))
    UptakeLayer_minus_InterCid=zeros((int(r_discret[0]),TotLayers,Nscenarios))
    STFmb=zeros((Nmb,TotLayers))
    STFcell_plus=zeros((Ncells,TotLayers))
    STFcell_minus=zeros((Ncells,TotLayers))
    STFlayer_plus=zeros((int(r_discret[0]),TotLayers))
    STFlayer_minus=zeros((int(r_discret[0]),TotLayers))
    PsiCellLayer=zeros((int(r_discret[0]),TotLayers,Nscenarios))
    PsiWallLayer=zeros((int(r_discret[0]),TotLayers,Nscenarios))
    #UptakeDistri_plus=zeros((40,3,8))#the size will be adjusted, but won't be more than 40. Dimension 1: radial position, 2: compartment, 3: scenario
    #UptakeDistri_minus=zeros((40,3,8))
    Hydropatterning[:]=NAN
    Hydrotropism[:]=NAN
    
    #iMaturity=-1 #Iteration index in the Barriers loop
    TopLayer=0
    for iMaturity in range(Nmaturity_loop):
        if PileUp==2: #Will only take one lap (including all maturity levels simultaneously) in the Maturity loop
            TopLayer=TotLayers
            AxialLayers=TotLayers
            Barriers=empty((Nmaturity))
            heights=empty((AxialLayers,1))
            Cumheights=zeros((AxialLayers,1))
            MaturL=empty((AxialLayers))
            BarriersL=empty((AxialLayers))
            iM=0
            iLayer=0
            for tempM in Maturityrange:
                Barriers[iM]=int(tempM.get("Barrier"))
                BarriersL[iLayer:iLayer+int(Nlayers[iM])]=int(tempM.get("Barrier"))
                MaturL[iLayer:iLayer+int(Nlayers[iM])]=iM
                heights[iLayer:iLayer+int(Nlayers[iM])]=int(tempM.get("height"))
                for iL in range(iLayer,iLayer+int(Nlayers[iM])):
                    try:
                        Cumheights[iL+1]=Cumheights[iL]+heights[iL]
                    except:
                        print('Cumheights done')
                iLayer+=int(Nlayers[iM])
                iM+=1
        else:
            AxialLayers=int(Nlayers[iMaturity])
            MaturL=zeros((AxialLayers,1))
            Barrier=int(Maturityrange[iMaturity].get("Barrier")) #Apoplastic barriers (0: No apoplastic barrier, 1:Endodermis radial walls, 2:Endodermis with passage cells, 3: Endodermis full, 4: Endodermis full and exodermis radial walls)
            BarriersL=Barrier*ones((AxialLayers,1))
            height=int(Maturityrange[iMaturity].get("height")) #Cell length in the axial direction (microns)
            AxialLayers=int(Maturityrange[iMaturity].get("Nlayers")) #Total number of stacked cross-sections at this maturity level
            TopLayer+=AxialLayers
            Cumheights=zeros((int(Nlayers[iMaturity]),1))
            for iL in range(1,int(Nlayers[iMaturity])):
                Cumheights[iL]=Cumheights[iL-1]+height
        
        #Index for barriers loop
        if PileUp==2: #Will only take one lap (including all maturity levels simultaneously) in the Maturity loop
            print('All maturity levels explicitly coupled')
            #iMaturity=0
        else:
            #iMaturity+=1
            print('Maturity #'+str(iMaturity)+' with apoplastic barrier type #'+str(Barrier))
        
        #Scenarios concern boundary conditions only
        count=0
        print('Scenario #'+str(count))
        
        #Soil, xylem, and phloem pressure potentials
        Psi_xyl[1][iMaturity][count]=float(BC_xyl_range[count].get("pressure_prox")) #Xylem pressure potential at proximal end (hPa)
        Psi_xyl[0][iMaturity][count]=float(BC_xyl_range[count].get("pressure_dist")) #Xylem pressure potential at distal end (hPa)
        dPsi_xyl[iMaturity][count]=float(BC_xyl_range[count].get("deltaP_prox")) #Xylem pressure potential change as compared to equilibrium pressure (hPa)
        Flow_xyl[1][0][count]=float(BC_xyl_range[count].get("flowrate_prox")) #Xylem flow rate on proximal end of the segment (cm^3/d)
        Flow_xyl[0][0][count]=float(BC_xyl_range[count].get("flowrate_dist")) #Xylem flow rate on distal end of the segment (cm^3/d)
        i_dist=0
        no_dist_BCx=False
        if not isnan(Psi_xyl[0][iMaturity][count]):
            i_dist+=1
        if not isnan(Flow_xyl[0][0][count]): #dist&prox, Nxyl+1, Nscenarios
            i_dist+=1
        if i_dist==0:
            no_dist_BCx=True
        elif i_dist>1:
            error('May not have more than 1 water BC type at the distal end of xylem vessels')
        i_prox=0
        no_prox_BCx=False
        if not isnan(Psi_xyl[1][iMaturity][count]):
            i_prox+=1
        if not isnan(dPsi_xyl[iMaturity][count]):
            i_prox+=1
        if not isnan(Flow_xyl[1][0][count]):
            i_prox+=1
        if i_prox==0:
            no_prox_BCx=True
        elif i_prox>1:
            error('May not have more than 1 water BC type at the proximal end of xylem vessels')
        if not isnan(Flow_xyl[1][0][count]):
            tot_flow_prox=Flow_xyl[1][0][count]
            sum_area=0
            i=0
            for cid in listxyl:
                area=cellarea[cid-NwallsJun]
                Flow_xyl[1][i+1][count]=tot_flow_prox*area
                sum_area+=area
                i+=1
            i=0
            for cid in listxyl:
                Flow_xyl[1][i+1][count]/=sum_area #Total xylem flow rate partitioned proportionnally to xylem cross-section area
                i+=1
            if Flow_xyl[1][0][count]==0.0 and Flow_xyl[0][0][count]==0.0:
                iEquil_xyl=count
        elif not isnan(dPsi_xyl[iMaturity][count]):
            if not isnan(iEquil_xyl):
                Psi_xyl[1][iMaturity][count]=Psi_xyl[1][iMaturity][iEquil_xyl]+dPsi_xyl[iMaturity][count]
            else:
                print('Error: Cannot have xylem pressure change relative to equilibrium without having a prior scenario with equilibrium xylem boundary condition')
        if not isnan(Flow_xyl[0][0][count]):
            tot_flow_dist=Flow_xyl[0][0][count]
            sum_area=0
            i=0
            for cid in listxyl:
                area=cellarea[cid-NwallsJun]
                Flow_xyl[0][i+1][count]=tot_flow_dist*area
                sum_area+=area
                i+=1
            i=0
            for cid in listxyl: 
                Flow_xyl[0][i+1][count]/=sum_area #Total xylem flow rate partitioned proportionnally to xylem cross-section area
                i+=1
        
        Psi_sieve[1][iMaturity][count]=float(BC_sieve_range[count].get("pressure_prox")) #Phloem sieve element pressure potential (hPa) at proximal end
        Psi_sieve[0][iMaturity][count]=float(BC_sieve_range[count].get("pressure_dist")) #Phloem sieve element pressure potential (hPa) at distal end
        dPsi_sieve[iMaturity][count]=float(BC_sieve_range[count].get("deltaP_prox")) #Phloem pressure potential change as compared to equilibrium pressure (hPa)
        Flow_sieve[1][0][count]=float(BC_sieve_range[count].get("flowrate_prox")) #Phloem flow rate (cm^3/d) at proximal end
        Flow_sieve[0][0][count]=float(BC_sieve_range[count].get("flowrate_dist")) #Phloem flow rate (cm^3/d) at distal end
        i_dist=0
        no_dist_BCp=False
        if not isnan(Psi_sieve[0][iMaturity][count]):
            i_dist+=1
        if not isnan(Flow_sieve[0][0][count]):
            i_dist+=1
        if i_dist==0:
            no_dist_BCp=True
        elif i_dist>1:
            error('May not have more than 1 water BC type at the distal end of sieve tubes')
        i_prox=0
        no_prox_BCp=False
        if not isnan(Psi_sieve[1][iMaturity][count]):
            i_prox+=1
        if not isnan(dPsi_sieve[iMaturity][count]):
            i_prox+=1
        if not isnan(Flow_sieve[1][0][count]):
            i_prox+=1
        if i_prox==0:
            no_prox_BCp=True
        elif i_prox>1:
            error('May not have more than 1 water BC type at the proximal end of sieve tubes')
        XylOK=0
        if PileUp==2:
            if Barriers[-1]>0:
                XylOK=1
        else:
            if Barrier>0:
                XylOK=1
        if not isnan(Flow_sieve[1][0][count]):
            tot_flow=Flow_sieve[1][0][count]
            sum_area=0
            i=1
            if XylOK:
                for cid in listsieve:
                    area=cellarea[cid-NwallsJun]
                    Flow_sieve[1][i][count]=tot_flow*area
                    sum_area+=area
                    i+=1
                i=1
                for cid in listsieve:
                    Flow_sieve[1][i][count]/=sum_area #Total phloem flow rate partitioned proportionnally to phloem cross-section area
                    i+=1
            else:
                for cid in listsieve:
                    if cid in listprotosieve:
                        area=cellarea[cid-NwallsJun]
                        Flow_sieve[1][i][count]=tot_flow*area
                        sum_area+=area
                    else:
                        Flow_sieve[1][i][count]=0
                    i+=1
                i=1
                for cid in listsieve:
                    if cid in listprotosieve:
                        Flow_sieve[1][i][count]/=sum_area #Total phloem flow rate partitioned proportionnally to phloem cross-section area
                    i+=1
            if Flow_sieve[1][0][count]==0.0:
                iEquil_sieve=count
        elif not isnan(dPsi_sieve[iMaturity][count]):
            if not isnan(iEquil_sieve):
                Psi_sieve[1][iMaturity][count]=Psi_sieve[1][iMaturity][iEquil_sieve]+dPsi_sieve[iMaturity][count]
            else:
                print('Error: Cannot have phloem pressure change relative to equilibrium without having a prior scenario with equilibrium phloem boundary condition')
        if not isnan(Flow_sieve[0][0][count]):
            tot_flow=Flow_sieve[0][0][count]
            sum_area=0
            if PileUp==2:
                if Barriers[0]==0:
                    i=1
                    for cid in listsieve:
                        if cid in listprotosieve:
                            area=cellarea[cid-NwallsJun]
                            Flow_sieve[0][i][count]=tot_flow*area
                            sum_area+=area
                        else:
                            Flow_sieve[0][i][count]=0
                        i+=1
                    i=1
                    for cid in listsieve:
                        if cid in listprotosieve:
                            Flow_sieve[0][i][count]/=sum_area #Total phloem flow rate partitioned proportionnally to phloem cross-section area
                        i+=1
                else:
                    i=1
                    for cid in listsieve:
                        area=cellarea[cid-NwallsJun]
                        Flow_sieve[0][i][count]=tot_flow*area
                        sum_area+=area
                        i+=1
                    i=1
                    for cid in listsieve:
                        Flow_sieve[0][i][count]/=sum_area #Total phloem flow rate partitioned proportionnally to phloem cross-section area
                        i+=1
            else:
                if Barrier==0:
                    i=1
                    for cid in listsieve:
                        if cid in listprotosieve:
                            area=cellarea[cid-NwallsJun]
                            Flow_sieve[0][i][count]=tot_flow*area
                            sum_area+=area
                        i+=1
                    i=1
                    for cid in listsieve:
                        if cid in listprotosieve:
                            Flow_sieve[0][i][count]/=sum_area #Total phloem flow rate partitioned proportionnally to phloem cross-section area
                        i+=1
                else:
                    i=1
                    for cid in listsieve:
                        area=cellarea[cid-NwallsJun]
                        Flow_sieve[0][i][count]=tot_flow*area
                        sum_area+=area
                        i+=1
                    i=1
                    for cid in listsieve:
                        Flow_sieve[0][i][count]/=sum_area #Total phloem flow rate partitioned proportionnally to phloem cross-section area
                        i+=1
        
        #Soil - root contact limit
        if PileUp==2:
            Xcontacts=empty((Nmaturity))
            if NXcontact is Nmaturity:
                iM=0
                for tempM in Maturityrange:
                    Xcontacts[iM]=float(Xcontactrange[iM].get("value"))
                    iM+=1
            elif NXcontact is Nhydraulics:
                Xcontacts[:]=float(Xcontactrange[h].get("value")) #(micrometers) X threshold coordinate of contact between soil and root (lower X not in contact with soil)
            elif NXcontact == 1:
                Xcontacts[:]=float(Xcontactrange[0].get("value"))
            else:
                error('PileUp=2 not compatible with NXcontact other than equal to 1, Nhydraulics, or Nmaturity')
        else:
            if NXcontact is Nhydraulics:
                Xcontact=float(Xcontactrange[h].get("value")) #(micrometers) X threshold coordinate of contact between soil and root (lower X not in contact with soil)
            elif NXcontact == 1:
                Xcontact=float(Xcontactrange[0].get("value"))
            else:
                Xcontact=float(Xcontactrange[int(h/(NkAQP*NKpl*Nkw*Nkw_barrier))].get("value")) #OK
        
        #Cell wall hydraulic conductivity
        if Nkw is Nhydraulics:
            kw = float(kwrange[h].get("value"))
        elif Nkw == 1:
            kw = float(kwrange[0].get("value"))
        else:
            kw = float(kwrange[int(h/(NkAQP*NKpl))%Nkw].get("value"))
        if PileUp==2:
            kw_barriers=empty((Nmaturity,2))
            if Nkw_barrier is Nmaturity:
                iM=0
                for tempM in Maturityrange:
                    kw_barriers[iM,0]=float(kw_barrier_range[iM].get("Casp"))
                    kw_barriers[iM,1]=float(kw_barrier_range[iM].get("Sub"))
                    iM+=1
            elif Nkw_barrier is Nhydraulics:
                kw_barriers[:,0] = float(kw_barrier_range[h].get("Casp"))
                kw_barriers[:,1] = float(kw_barrier_range[h].get("Sub"))
            elif Nkw_barrier == 1:
                kw_barriers[:,0] = float(kw_barrier_range[0].get("Casp"))
                kw_barriers[:,1] = float(kw_barrier_range[0].get("Sub"))
            else:
                kw_barriers[:,0] = float(kw_barrier_range[int(h/(NkAQP*NKpl*Nkw))%Nkw_barrier].get("Casp"))
                kw_barriers[:,1] = float(kw_barrier_range[int(h/(NkAQP*NKpl*Nkw))%Nkw_barrier].get("Sub"))
        else:
            kw_barrier=empty((1,2))
            if Nkw_barrier is Nhydraulics:
                kw_barrier[0,0] = float(kw_barrier_range[h].get("Casp"))
                kw_barrier[0,1] = float(kw_barrier_range[h].get("Sub"))
            elif Nkw_barrier == 1:
                kw_barrier[0,0] = float(kw_barrier_range[0].get("Casp"))
                kw_barrier[0,1] = float(kw_barrier_range[0].get("Sub"))
            else:
                kw_barrier[0,0] = float(kw_barrier_range[int(h/(NkAQP*NKpl*Nkw))%Nkw_barrier].get("Casp"))
                kw_barrier[0,1] = float(kw_barrier_range[int(h/(NkAQP*NKpl*Nkw))%Nkw_barrier].get("Sub"))
        #kw_barrier = kw/10.0
        if PileUp==2:
            kws_endo_endo=ones((Nmaturity))*kw
            kws_endo_endo[Barriers>0]=kw_barriers[Barriers>0,0]
            kws_puncture=ones((Nmaturity))*kw
            kws_exo_exo=ones((Nmaturity))*kw
            kws_exo_exo[Barriers>3]=kw_barriers[Barriers>3,0]
            kws_cortex_cortex=ones((Nmaturity))*kw
            kws_endo_peri=ones((Nmaturity))*kw
            kws_endo_peri[np.logical_and(Barriers>1,Barriers<5)]=kw_barriers[np.logical_and(Barriers>1,Barriers<5),1]
            kws_endo_cortex=kws_endo_peri
            kws_passage=ones((Nmaturity))*kw
            kws_passage[np.logical_or(Barriers==3,Barriers==4)]=kw_barriers[np.logical_or(Barriers==3,Barriers==4),1]
        else:
            if Barrier==0: #No Casparian strip ###Yet to come: Punctured Casparian strip as in Steudle et al. (1993)
                kw_endo_endo=kw
                kw_puncture=kw
                kw_exo_exo=kw #(cm^2/hPa/d) hydraulic conductivity of the suberised walls between exodermis cells
                kw_cortex_cortex=kw
                kw_endo_peri=kw #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
                kw_endo_cortex=kw #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
                kw_passage=kw #(cm^2/hPa/d) hydraulic conductivity of passage cells tangential walls
            elif Barrier==1: #Endodermis radial walls
                kw_endo_endo=kw_barrier[0,0]
                kw_exo_exo=kw #(cm^2/hPa/d) hydraulic conductivity of the suberised walls between exodermis cells
                kw_cortex_cortex=kw
                kw_endo_peri=kw #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
                kw_endo_cortex=kw #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
                kw_passage=kw #(cm^2/hPa/d) hydraulic conductivity of passage cells tangential walls
            elif Barrier==2: #Endodermis with passage cells
                kw_endo_endo=kw_barrier[0,0]
                kw_exo_exo=kw #(cm^2/hPa/d) hydraulic conductivity of the suberised walls between exodermis cells
                kw_cortex_cortex=kw
                kw_endo_peri=kw_barrier[0,1] #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
                kw_endo_cortex=kw_barrier[0,1] #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
                kw_passage=kw #(cm^2/hPa/d) hydraulic conductivity of passage cells tangential walls
            elif Barrier==3: #Endodermis full
                kw_endo_endo=kw_barrier[0,0]
                kw_exo_exo=kw #(cm^2/hPa/d) hydraulic conductivity of the suberised walls between exodermis cells
                kw_cortex_cortex=kw
                kw_endo_peri=kw_barrier[0,1] #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
                kw_endo_cortex=kw_barrier[0,1] #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
                kw_passage=kw_barrier[0,1] #(cm^2/hPa/d) hydraulic conductivity of passage cells tangential walls
            elif Barrier==4: #Endodermis full and exodermis radial walls
                kw_endo_endo=kw_barrier[0,0]
                kw_exo_exo=kw_barrier[0,0] #(cm^2/hPa/d) hydraulic conductivity of the suberised walls between exodermis cells
                kw_cortex_cortex=kw
                kw_endo_peri=kw_barrier[0,1] #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
                kw_endo_cortex=kw_barrier[0,1] #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
                kw_passage=kw_barrier[0,1] #(cm^2/hPa/d) hydraulic conductivity of passage cells tangential walls
            elif Barrier==5:
                kw_endo_endo=kw_barrier[0,0]
                kw_exo_exo=kw_barrier[0,0] #(cm^2/hPa/d) hydraulic conductivity of the suberised walls between exodermis cells
                kw_cortex_cortex=kw
                kw_endo_peri=kw #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
                kw_endo_cortex=kw #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
                kw_passage=kw #(cm^2/hPa/d) hydraulic conductivity of passage cells tangential walls
        
        #Plasmodesmatal hydraulic conductance
        if NKpl is Nhydraulics:
            iPD=h
        elif NKpl == 1:
            iPD=0
        else:
            iPD=int(h/NkAQP)%NKpl
        Kpl = float(Kplrange[iPD].get("value")) #cm^3/hPa/d/plasmodesmata
        
        #Contribution of aquaporins to membrane hydraulic conductivity
        if NkAQP is Nhydraulics:
            iAQP=h
        elif NkAQP == 1:
            iAQP=0
        else:
            iAQP=h%NkAQP
        kaqp = float(kAQPrange[iAQP].get("value"))
        kaqp_stele= kaqp*float(kAQPrange[iAQP].get("stele_factor"))
        kaqp_endo= kaqp*float(kAQPrange[iAQP].get("endo_factor"))
        kaqp_exo= kaqp*float(kAQPrange[iAQP].get("exo_factor"))
        kaqp_epi= kaqp*float(kAQPrange[iAQP].get("epi_factor"))
        kaqp_cortex= kaqp*float(kAQPrange[iAQP].get("cortex_factor"))
        
        #Calculate parameter a
        if ratio_cortex==1 or PileUp==2: #Uniform AQP activity in all cortex membranes
            a_cortex=0.0  #(1/hPa/d)
            b_cortex=kaqp_cortex #(cm/hPa/d)
            if not ratio_cortex==1:
                print('PileUp==2 not compatible with non-uniform aqp distribution in cortex; kaqp_cortex set to ',kaqp_cortex,' cm/hPa/d')
        else:
            tot_surf_cortex=0.0 #Total membrane exchange surface in cortical cells (square centimeters)
            temp=0.0 #Term for summation (cm3)
            for w in Cell2Wall_loop: #Loop on cells. Cell2Wall_loop contains cell wall groups info (one group by cell)
                cellnumber1 = int(w.getparent().get("id")) #Cell ID number
                for r in w: #Loop for wall elements around the cell
                    wid= int(r.get("id")) #Cell wall ID
                    if G.node[NwallsJun + cellnumber1]['cgroup']==4: #Cortex
                        dist_cell=sqrt(square(position[wid][0]-position[NwallsJun+cellnumber1][0])+square(position[wid][1]-position[NwallsJun+cellnumber1][1])) #distance between wall node and cell node (micrometers)
                        surf=(height+dist_cell)*lengths[wid]*1.0E-08 #(square centimeters)
                        temp+=surf*1.0E-04*(dist_grav[wid]+(ratio_cortex*dmax_cortex-dmin_cortex)/(1-ratio_cortex))
                        tot_surf_cortex+=surf
            a_cortex=kaqp_cortex*tot_surf_cortex/temp  #(1/hPa/d)
            b_cortex=a_cortex*1.0E-04*(ratio_cortex*dmax_cortex-dmin_cortex)/(1-ratio_cortex) #(cm/hPa/d)
        
        ######################
        ##Filling the matrix##
        ######################
        print('Filling matrix_W.2 & matrix_C.2')
        t00 = time.perf_counter()  #Lasts about 65 sec in the 24mm root
        if sparseM==1:
            if PileUp==2:
                matrix_W2 = lil_matrix((AxialLayers*Ntot, AxialLayers*Ntot),dtype=dp)
            else:
                matrix_W = lil_matrix((Ntot, Ntot),dtype=dp) #Initializes the Doussan matrix
                matrix_W2 = lil_matrix((AxialLayers*Ntot, AxialLayers*Ntot),dtype=dp)
        else:
            if AxialLayers>1:
                error('If there are several layers (3D problem, see *Geometry.in), sparse matrices have to be used (see *General.in)')
            else:
                matrix_W = np.zeros(((Ntot),Ntot)) #Initializes the Doussan matrix
        if PileUp==2:
            kwi=zeros((Ntot,Nmaturity)) #Stores wall conductivities
            kaqpi=zeros((Ncells,Nmaturity)) #Stores cell membrane conductivities
            Fpli=zeros((Ncells,Nmaturity)) #FRequncy of plasmodesmata in the axial direction (qtity per cm^2), used to calculate Kaxial
        else:
            kwi=zeros((Ntot,1)) #Stores wall conductivities
            kaqpi=zeros((Ncells,1)) #Stores cell membrane conductivities
            Fpli=zeros((Ncells,1)) #FRequncy of plasmodesmata in the axial direction (qtity per cm^2), used to calculate Kaxial
        for i in range(NwallsJun,Ntot): #Loop over cells
            if G.node[i]['cgroup']==1: #Exodermis
                if PileUp==2:
                    iLayer=0
                    for iM in range(Nmaturity):
                        kwi[i,iM] = kws_exo_exo[iM]
                        iLayer+=int(Nlayers[iM])
                    kaqpi[i-NwallsJun,:]=kaqp_exo
                    Fpli[i-NwallsJun,:]=Fplxheight/(Defheight*1.0E-04) #Quantity per cm^2
                else:
                    kwi[i] = kw_exo_exo
                    kaqpi[i-NwallsJun]=kaqp_exo
                    Fpli[i-NwallsJun]=Fplxheight/(Defheight*1.0E-04) #Quantity per cm^2
            elif G.node[i]['cgroup']==2: #Epidermis
                if PileUp==2:
                    kwi[i,:] = kw
                    kaqpi[i-NwallsJun,:]=kaqp_epi
                    Fpli[i-NwallsJun,:]=Fplxheight/(Defheight*1.0E-04)
                else:
                    kwi[i] = kw
                    kaqpi[i-NwallsJun]=kaqp_epi
                    Fpli[i-NwallsJun]=Fplxheight/(Defheight*1.0E-04)
            elif G.node[i]['cgroup']==3: #Endodermis
                if PileUp==2:
                    iLayer=0
                    for iM in range(Nmaturity):
                        kwi[i,iM] = kws_endo_endo[iM]
                        iLayer+=int(Nlayers[iM])
                    kaqpi[i-NwallsJun,:]=kaqp_endo
                    Fpli[i-NwallsJun,:]=Fplxheight_endo_endo/(Defheight*1.0E-04)
                else:
                    kwi[i] = kw_endo_endo
                    kaqpi[i-NwallsJun]=kaqp_endo
                    Fpli[i-NwallsJun]=Fplxheight_endo_endo/(Defheight*1.0E-04)
            elif G.node[i]['cgroup']==13 or G.node[i]['cgroup']==19 or G.node[i]['cgroup']==20: #xylem cell or vessel
                if PileUp==2:
                    iLayer=0
                    for iM in range(Nmaturity):
                        if Barriers[iM]>0: #Xylem vessel
                            kaqpi[i-NwallsJun,iM]=1.0 #No membrane resistance (very small) because no membrane
                            #kwi[iM] & Fpli[iM] remain null...
                        else:
                            kwi[i,iM] = kw
                            kaqpi[i-NwallsJun,iM]=kaqp_stele
                            Fpli[i-NwallsJun,iM]=Fplxheight_stele_stele/(Defheight*1.0E-04) #Quantity per cm^2
                        iLayer+=int(Nlayers[iM])
                else:
                    if Barrier>0: #Xylem vessel
                        kaqpi[i-NwallsJun]=1.0 #No membrane resistance (very small) because no membrane
                        #kwi & Fpli remain null...
                    else:
                        kwi[i] = kw
                        kaqpi[i-NwallsJun]=kaqp_stele
                        Fpli[i-NwallsJun]=Fplxheight_stele_stele/(Defheight*1.0E-04) #Quantity per cm^2
            elif G.node[i]['cgroup']>4: #Stele and pericycle but not xylem
                if PileUp==2:
                    kwi[i,:] = kw
                    kaqpi[i-NwallsJun,:]=kaqp_stele
                    if G.node[i]['cgroup']==26 or G.node[i]['cgroup']==24 or G.node[i]['cgroup']==12: #Companion Cell
                        Fpli[i-NwallsJun,:]=Fplxheight_stele_stele/(Defheight*1.0E-04)*float(Kplrange[iPD].get("PCC_factor")) #Quantity per cm^2
                    elif G.node[i]['cgroup']==23 or G.node[i]['cgroup']==11: #Phloem sieve tube
                        Fpli[i-NwallsJun,:]=Fplxheight_stele_stele/(Defheight*1.0E-04)*float(Kplrange[iPD].get("PST_factor")) #Quantity per cm^2
                    elif i-NwallsJun in PPP: #Phloem pole pericyle
                        Fpli[i-NwallsJun,:]=Fplxheight_stele_stele/(Defheight*1.0E-04)*float(Kplrange[iPD].get("PPP_factor"))
                    else:
                        Fpli[i-NwallsJun,:]=Fplxheight_stele_stele/(Defheight*1.0E-04)
                else:
                    kwi[i] = kw
                    kaqpi[i-NwallsJun]=kaqp_stele
                    if G.node[i]['cgroup']==26 or G.node[i]['cgroup']==24 or G.node[i]['cgroup']==12: #Companion Cell
                        Fpli[i-NwallsJun]=Fplxheight_stele_stele/(Defheight*1.0E-04)*float(Kplrange[iPD].get("PCC_factor")) #Quantity per cm^2
                    elif G.node[i]['cgroup']==23 or G.node[i]['cgroup']==11: #Phloem sieve tube
                        Fpli[i-NwallsJun]=Fplxheight_stele_stele/(Defheight*1.0E-04)*float(Kplrange[iPD].get("PST_factor")) #Quantity per cm^2
                    elif i-NwallsJun in PPP: #Phloem pole pericyle
                        Fpli[i-NwallsJun]=Fplxheight_stele_stele/(Defheight*1.0E-04)*float(Kplrange[iPD].get("PPP_factor"))
                    else:
                        Fpli[i-NwallsJun]=Fplxheight_stele_stele/(Defheight*1.0E-04)
            elif G.node[i]['cgroup']==4: #Cortex
                if PileUp==2:
                    iLayer=0
                    for iM in range(Nmaturity):
                        if Barriers[iM]>0 and (i-NwallsJun in InterCid): #Intercellular space "cell":
                            if Xylem_pieces:
                                kwi[i,iM] = kInterC
                            else:
                                kwi[i,iM] = kws_cortex_cortex[iM]
                            kaqpi[i-NwallsJun,iM]=kInterC
                            Fpli[i-NwallsJun,iM]=0 #Quantity per cm^2
                        else:
                            kwi[i,iM] = kws_cortex_cortex[iM]
                            if not isnan(x_grav):
                                dist_grav_cell=sqrt(square(position[i][0]-x_grav)+square(position[i][1]-y_grav))
                                kaqpi[i-NwallsJun,iM]=float(a_cortex*dist_grav_cell*1.0E-04+b_cortex) #AQP activity (cm/hPa/d)
                                if kaqpi[i-NwallsJun,iM] < 0:    
                                    error('Error, negative kaqp in cortical cell, adjust Paqp_cortex')
                            else:
                                kaqpi[i-NwallsJun,iM]=kaqp_cortex
                            Fpli[i-NwallsJun,iM]=Fplxheight_cortex_cortex/(Defheight*1.0E-04)*float(Kplrange[iPD].get("cortex_factor")) #Quantity per cm^2 including correction for specific cell-type PD aperture
                        iLayer+=int(Nlayers[iM])
                else:
                    if Barrier>0 and (i-NwallsJun in InterCid):
                        if Xylem_pieces:
                            kwi[i] = kInterC #At some point we might want to calculate it based on poiseuille to avoid having water-filled IS behaving like superxylem
                        else:
                            kwi[i] = kw_cortex_cortex
                        kaqpi[i-NwallsJun]=kInterC
                        Fpli[i-NwallsJun]=0 #Quantity per cm^2
                    else:
                        kwi[i] = kw_cortex_cortex
                        if not isnan(x_grav):
                            dist_grav_cell=sqrt(square(position[i][0]-x_grav)+square(position[i][1]-y_grav))
                            kaqpi[i-NwallsJun]=float(a_cortex*dist_grav_cell*1.0E-04+b_cortex) #AQP activity (cm/hPa/d)
                            if kaqpi[i-NwallsJun] < 0:    
                                error('Error, negative kaqp in cortical cell, adjust Paqp_cortex')
                        else:
                            kaqpi[i-NwallsJun]=kaqp_cortex
                        Fpli[i-NwallsJun]=Fplxheight_cortex_cortex/(Defheight*1.0E-04)*float(Kplrange[iPD].get("cortex_factor")) #Quantity per cm^2 including correction for specific cell-type PD aperture
        if Apo_Contagion==2 and Sym_Contagion==2:
            if sparseM==1:
                if PileUp==2:
                    matrix_C2 = lil_matrix((AxialLayers*Ntot, AxialLayers*Ntot),dtype=dp) #Initializes the matrix of convection diffusion
                else:
                    matrix_C = lil_matrix((Ntot, Ntot),dtype=dp)
                    matrix_C2 = lil_matrix((AxialLayers*Ntot, AxialLayers*Ntot),dtype=dp) #Initializes the matrix of convection diffusion
            else:
                if AxialLayers>1:
                    error('If there are several layers (3D problem, see *Geometry.in), sparse matrices have to be used (see *General.in)')
                else:
                    matrix_C = np.zeros(((Ntot),Ntot)) #Initializes the matrix of convection diffusion
            Volume_W2=zeros((AxialLayers*Ntot),dtype=dp)#empty(AxialLayers*Ntot,1)#
            if PileUp==2:
                for i in range(Nwalls):
                    Volume_W2[i::Ntot]=Water_fraction_apo*1.0E-12*(lat_dists[i][0]*thickness*lengths[i]+heights[:].T*thickness*lengths[i]/2-square(thickness)*lengths[i])*V_modifier
                    for temp in range(int(nWall2Junction[i])):
                        Volume_W2[int(Wall2Junction[i][temp])::Ntot]+=Volume_W2[i::Ntot]/4.0
                    Volume_W2[i::Ntot]/=2.0
                    #Decomposition rate (mol decomp/mol-day * cm^3)
                    for iLayer in range(AxialLayers):
                        matrix_C2[Ntot*iLayer+i,Ntot*iLayer+i]-=Degrad1*Volume_W2[i]
                for j in range(Nwalls,NwallsJun):
                    #Decomposition rate (mol decomp/mol-day * cm^3)
                    for iLayer in range(AxialLayers):
                        matrix_C2[Ntot*iLayer+j,Ntot*iLayer+j]-=Degrad1*Volume_W2[j]
                for cellnumber1 in range(Ncells):
                    iLayer=0
                    for iM in range(Nmaturity):
                        if Barriers[iM]>0 and (cellnumber1+NwallsJun in listxyl): #Xylem vessels are assumed full of water
                            Volume_W2[iLayer*Ntot+NwallsJun+cellnumber1:int(iLayer+Nlayers[iM])*Ntot:Ntot]=1.0E-12*cellarea[cellnumber1]*heights[iLayer:iLayer+int(Nlayers[iM])].T*V_modifier
                        else:
                            Volume_W2[iLayer*Ntot+NwallsJun+cellnumber1:int(iLayer+Nlayers[iM])*Ntot:Ntot]=Water_fraction_sym*1.0E-12*cellarea[cellnumber1]*heights[iLayer:iLayer+int(Nlayers[iM])].T*V_modifier
                        iLayer+=int(Nlayers[iM])
                    #Decomposition rate (mol decomp/mol-day * cm^3)
                    for iLayer in range(AxialLayers):
                        matrix_C2[Ntot*iLayer+NwallsJun+cellnumber1,Ntot*iLayer+NwallsJun+cellnumber1]-=Degrad1*Volume_W2[NwallsJun+cellnumber1]
            else:
                Volume_W=zeros((Ntot),dtype=dp)
                for i in range(Nwalls):
                    Volume_W[i]=Water_fraction_apo*1.0E-12*(lat_dists[i][0]*thickness*lengths[i]+height*thickness*lengths[i]/2-square(thickness)*lengths[i])*V_modifier
                    for temp in range(int(nWall2Junction[i])):
                        Volume_W[int(Wall2Junction[i][temp])]+=Volume_W[i]/4.0
                    Volume_W[i]/=2.0
                    #Decomposition rate (mol decomp/mol-day * cm^3)
                    if sparseM==1:
                        matrix_C[i,i]-=Degrad1*Volume_W[i]
                    else:
                        matrix_C[i][i]-=Degrad1*Volume_W[i]
                for j in range(Nwalls,NwallsJun):
                        #Decomposition rate (mol decomp/mol-day * cm^3)
                        if sparseM==1:
                            matrix_C[j,j]-=Degrad1*Volume_W[j]
                        else:
                            matrix_C[j][j]-=Degrad1*Volume_W[j]
                for cellnumber1 in range(Ncells):
                    if Barrier>0 and (cellnumber1+NwallsJun in listxyl): #Xylem vessels are assumed full of water
                        Volume_W[NwallsJun+cellnumber1]=1.0E-12*cellarea[cellnumber1]*height*V_modifier
                    else:
                        Volume_W[NwallsJun+cellnumber1]=Water_fraction_sym*1.0E-12*cellarea[cellnumber1]*height*V_modifier
                    #Decomposition rate (mol decomp/mol-day * cm^3)
                    if sparseM==1:
                        matrix_C[NwallsJun+cellnumber1,NwallsJun+cellnumber1]-=Degrad1*Volume_W[NwallsJun+cellnumber1]
                    else:
                        matrix_C[NwallsJun+cellnumber1][NwallsJun+cellnumber1]-=Degrad1*Volume_W[NwallsJun+cellnumber1]
                for iLayer in range(AxialLayers):
                    Volume_W2[iLayer*Ntot:(iLayer+1)*Ntot]=Volume_W
        elif Apo_Contagion==2:
            matrix_ApoC = np.zeros(((NwallsJun),NwallsJun)) #Initializes the matrix of convection
            rhs_ApoC = np.zeros((NwallsJun,1)) #Initializing the right-hand side matrix of solute apoplastic concentrations
            Volume_W=zeros((NwallsJun,1))
            for i in range(Nwalls):
                Volume_W[i]=Water_fraction_apo*1.0E-12*(lat_dists[i][0]*thickness*lengths[i]+height*thickness*lengths[i]/2-square(thickness)*lengths[i])*V_modifier
                for temp in range(int(nWall2Junction[i])):
                    Volume_W[int(Wall2Junction[i][temp])]+=Volume_W[i]/4.0
                    Volume_W[i]/=2.0
                if i in N_Dirichlet_ini:
                    matrix_ApoC[i][i]=1.0
                    rhs_ApoC[i][0]=cc_Dirichlet_ini[N_Dirichlet_ini.index(i)] #1 #Concentration in source wall i equals 1 by default
                else: #Decomposition rate (mol decomp/mol-day * cm^3)
                    matrix_ApoC[i][i]-=Degrad1*Volume_W[i]
            for j in range(Nwalls,NwallsJun):
                if j in N_Dirichlet_ini:
                    matrix_ApoC[j][j]=1.0
                    rhs_ApoC[j][0]=cc_Dirichlet_ini[N_Dirichlet_ini.index(j)] #1 #Concentration in source junction j equals 1 by default
                else: #Decomposition rate (mol decomp/mol-day * cm^3)
                    matrix_ApoC[j][j]-=Degrad1*Volume_W[j]
        elif Sym_Contagion==2:
            matrix_SymC = np.zeros(((Ncells),Ncells)) #Initializes the matrix of convection
            rhs_SymC = np.zeros((Ncells,1)) #Initializing the right-hand side matrix of solute symplastic concentrations
            for cellnumber1 in range(Ncells):
                if cellnumber1+NwallsJun in N_Dirichlet_ini:
                    matrix_SymC[cellnumber1][cellnumber1]=1.0
                    rhs_SymC[cellnumber1][0]=cc_Dirichlet_ini[N_Dirichlet_ini.index(cellnumber1+NwallsJun)] #1 #Concentration in source protoplasts equals 1 by default
                else: #Decomposition rate (mol decomp/mol-day * cm^3)
                    matrix_SymC[cellnumber1][cellnumber1]-=Water_fraction_sym*Degrad1*1.0E-12*cellarea[cellnumber1]*height
        
        jmb=0 #Index of membrane in Kmb
        if PileUp==2:
            Kmb=zeros((Nmb,Nmaturity)) #Stores membranes conductances for the second K loop
            K_axial=zeros((Ntot*AxialLayers,1)) #Vector of apoplastic and plasmodesmatal axial conductances
            Fraction_notTM=ones((Ncells,Nmaturity)) #Records the fraction of axial flow not across membranes
            DF_axial=zeros((Ntot*AxialLayers,1)) #Vector of apoplastic and plasmodesmatal axial diffusion coefficients
            if Kax_source==2:
                for K_xyl in K_xyl_range:
                    cellnumber=int(K_xyl.get("id"))
                    iLayer=0
                    for iM in range(Nmaturity):
                        if Barriers[iM]>0:
                            K_axial[int(iLayer*Ntot+NwallsJun+cellnumber):int(iLayer+Nlayers[iM])*Ntot:Ntot]=float(K_xyl.get("value"))
                        iLayer+=Nlayers[iM]
                K_xyl_spec=sum(K_axial[-Ntot:])*heights[-1,0]/1.0E04
                for temp in K_sieve_range:
                    K_sieve=float(temp.get("value"))
                    cellnumber=int(temp.get("id"))
                    iLayer=0
                    for iM in range(Nmaturity):
                        if Barriers[iM]>0 or (cellnumber+NwallsJun in listprotosieve):
                            K_axial[int(iLayer*Ntot+NwallsJun+cellnumber):int(iLayer+Nlayers[iM])*Ntot:Ntot]=K_sieve
                        iLayer+=Nlayers[iM]
            else: #Calculated from Poiseuille law (cm^3/hPa/d)
                for cid in listxyl:
                    iLayer=0
                    for iM in range(Nmaturity):
                        if Barriers[iM]>0:
                            K_axial[int(iLayer*Ntot+cid):int(iLayer+Nlayers[iM])*Ntot:Ntot]=cellarea[cid-NwallsJun]**2/(8*math.pi*heights[int(iLayer):int(iLayer+Nlayers[iM])]*1.0E-05/3600/24)*1.0E-12 #(micron^4/micron)->(cm^3) & (1.0E-3 Pa.s)->(1.0E-05/3600/24 hPa.d) 
                        iLayer+=Nlayers[iM]
                K_xyl_spec=sum(K_axial[-Ntot:])*heights[-1,0]/1.0E04
                for cid in listsieve:
                    iLayer=0
                    for iM in range(Nmaturity):
                        if Barriers[iM]>0 or (cid in listprotosieve):
                            K_axial[int(iLayer*Ntot+cid):int(iLayer+Nlayers[iM])*Ntot:Ntot]=cellarea[cid-NwallsJun]**2/(8*math.pi*heights[int(iLayer):int(iLayer+Nlayers[iM])]*1.0E-05/3600/24)*1.0E-12 #(micron^4/micron)->(cm^3) & (1.0E-3 Pa.s)->(1.0E-05/3600/24 hPa.d)
                        iLayer+=Nlayers[iM]
        else:
            Kmb=zeros((Nmb,1)) #Stores membranes conductances for the second K loop
            K_axial=zeros((Ntot,1)) #Vector of apoplastic and plasmodesmatal axial conductances
            Fraction_notTM=ones((Ncells,1)) #Records the fraction of axial flow not across membranes
            DF_axial=zeros((Ntot,1)) #Vector of apoplastic and plasmodesmatal axial diffusion coefficients
            if Barrier>0:
                if Kax_source==2:
                    for K_xyl in K_xyl_range:
                        cellnumber=int(K_xyl.get("id"))
                        K_axial[cellnumber+NwallsJun]=float(K_xyl.get("value"))
                    K_xyl_spec=sum(K_axial)*height/1.0E04
                    for K_sieve in K_sieve_range:
                        cellnumber=int(K_sieve.get("id"))
                        K_axial[cellnumber+NwallsJun]=float(K_sieve.get("value"))
                else: #Calculated from Poiseuille law (cm^3/hPa/d)
                    for cid in listxyl:
                        K_axial[cid]=cellarea[cid-NwallsJun]**2/(8*math.pi*height*1.0E-05/3600/24)*1.0E-12 #(micron^4/micron)->(cm^3) & (1.0E-3 Pa.s)->(1.0E-05/3600/24 hPa.d) 
                    K_xyl_spec=sum(K_axial)*height/1.0E04
                    for cid in listsieve:
                        K_axial[cid]=cellarea[cid-NwallsJun]**2/(8*math.pi*height*1.0E-05/3600/24)*1.0E-12 #(micron^4/micron)->(cm^3) & (1.0E-3 Pa.s)->(1.0E-05/3600/24 hPa.d) 
            else:
                if Kax_source==2:
                    for K_sieve in K_sieve_range:
                        cellnumber=int(K_sieve.get("id"))
                        if cellnumber+NwallsJun in listprotosieve:
                            K_axial[cellnumber+NwallsJun]=float(K_sieve.get("value"))
                else: #Calculated from Poiseuille law (cm^3/hPa/d)
                    for cid in listprotosieve:
                        K_axial[cid]=cellarea[cid-NwallsJun]**2/(8*math.pi*height*1.0E-05/3600/24)*1.0E-12 #(micron^4/micron)->(cm^3) & (1.0E-3 Pa.s)->(1.0E-05/3600/24 hPa.d) 
            
        
        #K_axial_MB=zeros((Ncells,1)) #Vector of transmembrane axial conductances -> We consider that there is not a single transmembrane conductance across two membranes as their respective flow rates could very well have opposite directions (ex: both from the symplast to the apoplast)
        #DF_axial_MB=zeros((Ncells,1)) #Vector of transmembrane axial diffusion coefficients
        list_ghostwalls=[] #"Fake walls" not to be displayed
        list_ghostjunctions=[] #"Fake junctions" not to be displayed
        nGhostJunction2Wall=0
        #Adding matrix components at cell-cell, cell-wall, and wall-junction connections
        for node, edges in G.adjacency_iter() : #adjacency_iter returns an iterator of (node, adjacency dict) tuples for all nodes. This is the fastest way to look at every edge. For directed graphs, only outgoing adjacencies are included.
            i=indice[node] #Node ID number
            #Here we count surrounding cell types in order to position apoplastic barriers
            count_endo=0 #total number of endodermis cells around the wall
            count_xyl=0 #total number of xylem cells around the wall
            count_stele_overall=0 #total number of stelar cells around the wall
            count_exo=0 #total number of exodermis cells around the wall
            count_epi=0 #total number of epidermis cells around the wall
            count_cortex=0 #total number of cortical cells around the wall
            count_passage=0 #total number of passage cells around the wall
            count_interC=0 #total number of intercellular spaces around the wall
            if i<Nwalls: #wall ID
                for neighbour, eattr in edges.items(): #Loop on connections (edges)
                    if eattr['path'] == 'membrane': #Wall connection
                        if any(passage_cell_ID==array((indice[neighbour])-NwallsJun)):
                            count_passage+=1
                        if any(InterCid==array((indice[neighbour])-NwallsJun)):
                            count_interC+=1
                            if Xylem_pieces:
                                if count_interC==2 and i not in list_ghostwalls:
                                    list_ghostwalls.append(i)
                        if G.node[neighbour]['cgroup']==3:#Endodermis
                            count_endo+=1
                        elif G.node[neighbour]['cgroup']==13 or G.node[neighbour]['cgroup']==19 or G.node[neighbour]['cgroup']==20:#Xylem cell or vessel
                            count_xyl+=1
                            if Xylem_pieces:
                                if count_xyl==2 and i not in list_ghostwalls:
                                    list_ghostwalls.append(i)
                        elif G.node[neighbour]['cgroup']>4:#Pericycle or stele but not xylem
                            count_stele_overall+=1
                        elif G.node[neighbour]['cgroup']==4:#Cortex
                            count_cortex+=1
                        elif G.node[neighbour]['cgroup']==1:#Exodermis
                            count_exo+=1
                        elif G.node[neighbour]['cgroup']==2:#Epidermis
                            count_epi+=1
            
            #Calculating axial conductances for each sub-element (still within the nodes loop)
            if i<Nwalls:
                if PileUp==2:
                    if (count_xyl==2 and Xylem_pieces): #"Fake wall" splitting a xylem cell in two
                        kwi[i,:] = 1.0E-16 #Non conductive
                        iLayer=0
                        for iM in range(Nmaturity):
                            K_axial[i+Ntot*iLayer:Ntot*(iLayer+int(Nlayers[iM])):Ntot]=wallXarea[i]*kwi[i,iM]#*K_axial_factor
                            if Apo_Contagion==2 and Sym_Contagion==2:
                                DF_axial[i+Ntot*iLayer:Ntot*(iLayer+int(Nlayers[iM])):Ntot]+=Diff_PW1*kwi[i,iM]/kw*1.0E-04*wallXarea[i]/heights[iLayer]#*DF_axial_factor
                            for j in Wall2Junction[i]:
                                j=int(j)
                                K_axial[j+Ntot*iLayer:Ntot*(iLayer+int(Nlayers[iM])):Ntot]+=1.0E-4*wallXarea[i]/2.0*kwi[i,iM]/heights[iLayer]#*K_axial_factor #(cm^3/hPa/d) WallXarea accounts for half of the total wall length. 2 quarters the wall length contributes to the 2 junctions, whose Xareas are thus half of the wallXarea
                                if Apo_Contagion==2 and Sym_Contagion==2:
                                    DF_axial[j+Ntot*iLayer:Ntot*(iLayer+int(Nlayers[iM])):Ntot]+=Diff_PW1*kwi[i,iM]/kw*1.0E-04*wallXarea[i]/2.0/heights[iLayer]#*DF_axial_factor
                            iLayer+=int(Nlayers[iM])
                        for j in Wall2Junction[i]:
                            if j not in list_ghostjunctions:
                                fakeJ=True
                                for ind in range(int(nJunction2Wall[j-Nwalls])):
                                    if Junction2Wall[j-Nwalls][ind] not in list_ghostwalls:
                                        fakeJ=False #If any of the surrounding walls is real, the junction is real
                                if fakeJ:
                                    list_ghostjunctions.append(j)
                                    nGhostJunction2Wall+=int(nJunction2Wall[j-Nwalls])+2 #The first and second thick junction nodes each appear twice in the text file for Paraview
                    elif (count_interC>=2): #"Fake wall" splitting an intercellular space if Xylem_pieces is true (due to segmentation of aerenchyma in CellSet)
                        iLayer=0
                        for iM in range(Nmaturity):
                            if Barriers[iM]>0:
                                if Xylem_pieces:
                                    kwi[i,iM] = 1.0E-16 #Non conductive
                                elif count_cortex>=2: #wall between two cortical cells:
                                    kwi[i,iM] = kws_cortex_cortex[iM]
                                else:
                                    error('Intercellular space wall not inside cortex??')
                                K_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot]=wallXarea[i]*kwi[i,iM]#*K_axial_factor
                                if Apo_Contagion==2 and Sym_Contagion==2:
                                    DF_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot]+=Diff_PW1*kwi[i,iM]/kw*1.0E-04*wallXarea[i]/heights[iLayer]#*DF_axial_factor
                                for j in Wall2Junction[i]:
                                    j=int(j)
                                    K_axial[iLayer*Ntot+j:int(iLayer+Nlayers[iM])*Ntot:Ntot]+=1.0E-4*wallXarea[i]/2.0*kwi[i,iM]/heights[iLayer]#*K_axial_factor #(cm^3/hPa/d) WallXarea accounts for half of the total wall length. 2 quarters the wall length contributes to the 2 junctions, whose Xareas are thus half of the wallXarea
                                    if Apo_Contagion==2 and Sym_Contagion==2:
                                        DF_axial[iLayer*Ntot+j:int(iLayer+Nlayers[iM])*Ntot:Ntot]+=Diff_PW1*kwi[i,iM]/kw*1.0E-04*wallXarea[i]/2.0/heights[iLayer]#*DF_axial_factor
                                    if Xylem_pieces and j not in list_ghostjunctions:
                                        fakeJ=True
                                        for ind in range(int(nJunction2Wall[j-Nwalls])):
                                            if Junction2Wall[j-Nwalls][ind] not in list_ghostwalls:
                                                fakeJ=False #If any of the surrounding walls is real, the junction is real
                                        if fakeJ:
                                            list_ghostjunctions.append(j)
                                            nGhostJunction2Wall+=int(nJunction2Wall[j-Nwalls])+2 #The first and second thick junction nodes each appear twice in the text file for Paraview
                            else:
                                if count_cortex>=2: #wall between two cortical cells
                                    kwi[i,iM] = kws_cortex_cortex[iM]
                                else:
                                    error('Intercellular space wall not inside cortex??')
                                K_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot]=1.0E-4*wallXarea[i]*kwi[i,iM]/heights[iLayer]#*K_axial_factor
                                if Apo_Contagion==2 and Sym_Contagion==2:
                                    DF_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot]+=Diff_PW1*kwi[i,iM]/kw*1.0E-04*wallXarea[i]/heights[iLayer]#*DF_axial_factor
                                for j in Wall2Junction[i]:
                                    j=int(j)
                                    K_axial[iLayer*Ntot+j:int(iLayer+Nlayers[iM])*Ntot:Ntot]+=1.0E-4*wallXarea[i]/2.0*kwi[i,iM]/heights[iLayer]#*K_axial_factor
                                    if Apo_Contagion==2 and Sym_Contagion==2:
                                        DF_axial[iLayer*Ntot+j:int(iLayer+Nlayers[iM])*Ntot:Ntot]+=Diff_PW1*kwi[i,iM]/kw*1.0E-04*wallXarea[i]/2.0/heights[iLayer]#*DF_axial_factor
                            iLayer+=int(Nlayers[iM])
                    else:
                        iLayer=0
                        for iM in range(Nmaturity):
                            if count_cortex>=2: #wall between two cortical cells
                                kwi[i,iM] = kws_cortex_cortex[iM]
                            elif count_endo>=2: #wall between two endodermis cells
                                kwi[i,iM] = kws_endo_endo[iM]
                            elif count_exo>=2: #wall between two exodermis cells
                                kwi[i,iM] = kws_exo_exo[iM]
                            else: #other walls
                                kwi[i,iM] = kw
                            K_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot]=1.0E-4*wallXarea[i]*kwi[i,iM]/heights[iLayer]#*K_axial_factor
                            if Apo_Contagion==2 and Sym_Contagion==2:
                                DF_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot]+=Diff_PW1*kwi[i,iM]/kw*1.0E-04*wallXarea[i]/heights[iLayer]#*DF_axial_factor
                            for j in Wall2Junction[i]:
                                j=int(j)
                                K_axial[iLayer*Ntot+j:int(iLayer+Nlayers[iM])*Ntot:Ntot]+=1.0E-4*wallXarea[i]/2.0*kwi[i,iM]/heights[iLayer]#*K_axial_factor
                                if Apo_Contagion==2 and Sym_Contagion==2:
                                    DF_axial[iLayer*Ntot+j:int(iLayer+Nlayers[iM])*Ntot:Ntot]+=Diff_PW1*kwi[i,iM]/kw*1.0E-04*wallXarea[i]/2.0/heights[iLayer]#*DF_axial_factor
                            iLayer+=int(Nlayers[iM])
                else:
                    if Xylem_pieces and ((count_interC>=2 and Barrier>0) or count_xyl==2): #"Fake wall" splitting an intercellular space or a xylem cell in two
                        kwi[i] = 1.0E-16 #Non conductive
                        K_axial[i]=wallXarea[i]*kwi[i]#*K_axial_factor
                        if Apo_Contagion==2 and Sym_Contagion==2:
                            temp_factor=kwi[i]/kw
                            DF_axial[i]+=Diff_PW1*temp_factor*1.0E-04*wallXarea[i]/height#*DF_axial_factor
                        for j in Wall2Junction[i]:
                            j=int(j)
                            K_axial[j]+=1.0E-4*wallXarea[i]/2.0*kwi[i]/height#*K_axial_factor #(cm^3/hPa/d) WallXarea accounts for half of the total wall length. 2 quarters the wall length contributes to the 2 junctions, whose Xareas are thus half of the wallXarea
                            if Apo_Contagion==2 and Sym_Contagion==2:
                                DF_axial[j]+=Diff_PW1*temp_factor*1.0E-04*wallXarea[i]/2.0/height#*DF_axial_factor
                            if j not in list_ghostjunctions:
                                fakeJ=True
                                for ind in range(int(nJunction2Wall[j-Nwalls])):
                                    if Junction2Wall[j-Nwalls][ind] not in list_ghostwalls:
                                        fakeJ=False #If any of the surrounding walls is real, the junction is real
                                if fakeJ:
                                    list_ghostjunctions.append(j)
                                    nGhostJunction2Wall+=int(nJunction2Wall[j-Nwalls])+2 #The first and second thick junction nodes each appear twice in the text file for Paraview
                    else:
                        if count_cortex>=2: #wall between two cortical cells
                            kwi[i] = kw_cortex_cortex
                        elif count_endo>=2: #wall between two endodermis cells
                            kwi[i] = kw_endo_endo
                        elif count_exo>=2: #wall between two exodermis cells
                            kwi[i] = kw_exo_exo
                        else: #other walls
                            kwi[i] = kw
                        K_axial[i]=1.0E-4*wallXarea[i]*kwi[i]/height#*K_axial_factor
                        if Apo_Contagion==2 and Sym_Contagion==2:
                            temp_factor=kwi[i]/kw
                            DF_axial[i]+=Diff_PW1*temp_factor*1.0E-04*wallXarea[i]/height#*DF_axial_factor
                        for j in Wall2Junction[i]:
                            j=int(j)
                            K_axial[j]+=1.0E-4*wallXarea[i]/2.0*kwi[i]/height#*K_axial_factor
                            if Apo_Contagion==2 and Sym_Contagion==2:
                                DF_axial[j]+=Diff_PW1*temp_factor*1.0E-04*wallXarea[i]/2.0/height#*DF_axial_factor
            
            elif i>=NwallsJun:
                if PileUp==2:
                    iLayer=0
                    for iM in range(Nmaturity):
                        if not ((Barriers[iM]>0 and (G.node[i]['cgroup']==13 or G.node[i]['cgroup']==19 or G.node[i]['cgroup']==20 or i in listsieve)) or (Barriers[iM]==0 and i in listprotosieve)): #neither mature xylem vessel nor mature sieve tube
                            if G.node[i]['cgroup']==3:#Endodermis
                                K_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot] = (Kpl*Fpli[i-NwallsJun,iM] + 1/(1/(kws_endo_peri[iM]/(thickness*1.0E-04))+1/(kmb+kaqpi[i-NwallsJun,iM]))) * 1.0E-08*cellarea[i-NwallsJun]#*K_axial_factor #[cm^3/hPa/d]
                                Fraction_notTM[i-NwallsJun,iM]=Kpl*Fpli[i-NwallsJun,iM]/(Kpl*Fpli[i-NwallsJun,iM] + 1/(1/(kws_endo_peri[iM]/(thickness*1.0E-04))+1/(kmb+kaqpi[i-NwallsJun,iM])))
                            else:
                                K_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot] = (Kpl*Fpli[i-NwallsJun,iM] + 1/(1/(kwi[i,iM]/(thickness*1.0E-04))+1/(kmb+kaqpi[i-NwallsJun,iM]))) * 1.0E-08*cellarea[i-NwallsJun]#*K_axial_factor #[cm^3/hPa/d]
                                Fraction_notTM[i-NwallsJun,iM]=Kpl*Fpli[i-NwallsJun,iM]/(Kpl*Fpli[i-NwallsJun,iM] + 1/(1/(kwi[i,iM]/(thickness*1.0E-04))+1/(kmb+kaqpi[i-NwallsJun,iM])))
                            if Apo_Contagion==2 and Sym_Contagion==2:
                                DF_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot]+=PD_section*1.0E-08*Fpli[i-NwallsJun,iM]*cellarea[i-NwallsJun]/thickness*1.0E-04*Diff_PD1#*DF_axial_factor #"Diffusive flux": Total PD cross-section area (micron^2) per unit PD length (micron) (tunred into cm) multiplied by solute diffusivity (cm^2/d) (yields cm^3/d)
                        else:#if not (Barriers[iM]>0 and (i-NwallsJun in InterCid)): changed 17/03/2020 -> Intercid are all accepted in first "if". Need to make sure that their Fpli is null
                            if (Barriers[iM]>0 and Barriers[iM-1]==0):
                                if G.node[i]['cgroup']==13 or G.node[i]['cgroup']==19 or G.node[i]['cgroup']==20: #First mature xylem vessel (its K_axial was already written, and the former Kaxial needs to be corrected to account for the single membrane)
                                    K_axial[(iLayer-1)*Ntot+i] = 1/(1/(kwi[i,iM-1]/(thickness*1.0E-04))+1/(kmb+kaqpi[i-NwallsJun,iM]))*1.0E-08*cellarea[i-NwallsJun]
                            if Apo_Contagion==2 and Sym_Contagion==2:
                                #Longitudinal diffusion term in xylem and phloem, assuming the same diffusion coefficient
                                DF_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot]+=cellarea[i-NwallsJun]/heights[iLayer]*1.0E-04*Diff_X1#*DF_axial_factor
                        iLayer+=int(Nlayers[iM])
                else:
                    if not ((Barrier>0 and (G.node[i]['cgroup']==13 or G.node[i]['cgroup']==19 or G.node[i]['cgroup']==20 or i in listsieve)) or Barrier==0 and i in listprotosieve): #xylem cell or vessel
                        if G.node[i]['cgroup']==3:#Endodermis
                            K_axial[i] = (Kpl*Fpli[i-NwallsJun] + 1/(1/(kw_endo_peri/(thickness*1.0E-04))+1/(kmb+kaqpi[i-NwallsJun]))) * 1.0E-08*cellarea[i-NwallsJun]#*K_axial_factor #[cm^3/hPa/d]
                            Fraction_notTM[i-NwallsJun]=Kpl*Fpli[i-NwallsJun]/(Kpl*Fpli[i-NwallsJun] + 1/(1/(kw_endo_peri/(thickness*1.0E-04))+1/(kmb+kaqpi[i-NwallsJun])))
                        else:
                            if kwi[i]==0.0: #air-filled intercellular space or aerenchyma
                                K_axial[i] = 0.0
                                Fraction_notTM[i-NwallsJun] = 1.0
                            else:
                                K_axial[i] = (Kpl*Fpli[i-NwallsJun] + 1/(1/(kwi[i]/(thickness*1.0E-04))+1/(kmb+kaqpi[i-NwallsJun]))) * 1.0E-08*cellarea[i-NwallsJun]#*K_axial_factor #[cm^3/hPa/d]
                                Fraction_notTM[i-NwallsJun]=Kpl*Fpli[i-NwallsJun]/(Kpl*Fpli[i-NwallsJun] + 1/(1/(kwi[i]/(thickness*1.0E-04))+1/(kmb+kaqpi[i-NwallsJun])))
                            #R_axial_MB = 2/(1.0E-08*cellarea[i-NwallsJun]*(kaqpi[i-NwallsJun]+kmb)) #[hPa*d/cm^3] Two membranes have to be crossed
                            #R_axial_PW = thickness*1.0E-04/(kwi[i]*cellarea[i-NwallsJun]) #[hPa*d/cm^3] A cell wall thickness has to be crossed
                            #K_axial_MB[i-NwallsJun] = 1/(R_axial_MB+R_axial_PW) #[cm^3/hPa/d]
                            #K_axial[i]+=K_axial_MB[i-NwallsJun]
                        if Apo_Contagion==2 and Sym_Contagion==2:
                            DF_axial[i]+=PD_section*1.0E-08*Fpli[i-NwallsJun]*cellarea[i-NwallsJun]/thickness*1.0E-04*Diff_PD1#*DF_axial_factor #"Diffusive flux": Total PD cross-section area (micron^2) per unit PD length (micron) (tunred into cm) multiplied by solute diffusivity (cm^2/d) (yields cm^3/d)
                            #DF_axial_MB[i-NwallsJun]=1/(1/(cellarea[i-NwallsJun]/thickness*1.0E-04*Diff_PW1)+1/(Diff_MB1*K_MB1*1.0E-04*cellarea[i-NwallsJun]/(2*dx_MB1))) #Membrane area to thickness ratio (cm) times partition coefficient of membrane K_MB1
                    else:#elif not (Barrier>0 and (i-NwallsJun in InterCid)): changed 17/03/2020 -> Intercid are all accepted in first "if". Need to make sure that their Fpli is null
                        if Apo_Contagion==2 and Sym_Contagion==2:
                            #Longitudinal diffusion term in xylem and phloem, assuming the same diffusion coefficient
                            DF_axial[i]+=cellarea[i-NwallsJun]/height*1.0E-04*Diff_X1#*DF_axial_factor
            
            for neighbour, eattr in edges.items(): #Loop on connections (edges)
                j = (indice[neighbour]) #neighbouring node number
                if j > i: #Only treating the information one way to save time
                    path = eattr['path'] #eattr is the edge attribute (i.e. connection type)
                    K=empty((AxialLayers,1))
                    if path == 'wall': #Wall connection
                        if PileUp==2:
                            iLayer=0
                            temp=empty((Nmaturity))
                            for iM in range(Nmaturity):
                                temp[iM]=1.0E-04*((eattr['lat_dist']*Xwalls+heights[iLayer]-thickness)*thickness)/eattr['length'] #Wall section to length ratio (cm)
                                K[iLayer:iLayer+int(Nlayers[iM])] = kwi[i,iM]*temp[iM] #Junction-Wall conductance (cm^3/hPa/d)
                                iLayer+=int(Nlayers[iM])
                            ########Solute fluxes (diffusion across walls and junctions)
                            if Apo_Contagion==2 and Sym_Contagion==2: #Sym & Apo contagion
                                DF=empty((AxialLayers))
                                iLayer=0
                                for iM in range(Nmaturity):
                                    DF[iLayer:iLayer+int(Nlayers[iM])]=temp[iM]*kwi[i,iM]/kw*Diff_PW1 #"Diffusive flux" (cm^3/d) temp is the section to length ratio of the wall to junction path
                                    iLayer+=int(Nlayers[iM])
                                for iLayer in range(AxialLayers):
                                    matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] -= DF[iLayer]
                                    matrix_C2[iLayer*Ntot+i,iLayer*Ntot+j] += DF[iLayer] #Convection will be dealt with further down
                                    matrix_C2[iLayer*Ntot+j,iLayer*Ntot+j] -= DF[iLayer] #temp_factor is the factor for reduced diffusion across impermeable walls
                                    matrix_C2[iLayer*Ntot+j,iLayer*Ntot+i] += DF[iLayer]
                        else:
                            temp=1.0E-04*((eattr['lat_dist']*Xwalls+height-thickness)*thickness)/eattr['length'] #Wall section to length ratio (cm)
                            K = kwi[i]*temp #Junction-Wall conductance (cm^3/hPa/d)
                            ########Solute fluxes (diffusion across walls and junctions)
                            if Apo_Contagion==2:
                                temp_factor=kwi[i]/kw #Factor for reduced diffusion across impermeable walls
                                DF=temp*temp_factor*Diff_PW1 #"Diffusive flux" (cm^3/d) temp is the section to length ratio of the wall to junction path
                                if Sym_Contagion==2: #Sym & Apo contagion
                                        if sparseM==1:
                                            matrix_C[i,i] -= DF
                                            matrix_C[i,j] += DF #Convection will be dealt with further down
                                        else:
                                            matrix_C[i][i] -= DF
                                            matrix_C[i][j] += DF #Convection will be dealt with further down
                                        if sparseM==1:
                                            matrix_C[j,j] -= DF #temp_factor is the factor for reduced diffusion across impermeable walls
                                            matrix_C[j,i] += DF
                                        else:
                                            matrix_C[j][j] -= DF #temp_factor is the factor for reduced diffusion across impermeable walls
                                            matrix_C[j][i] += DF
                                else: #Only Apo contagion
                                    if i not in N_Dirichlet_ini:
                                        matrix_ApoC[i][i] -= DF
                                        matrix_ApoC[i][j] += DF
                                    if j not in N_Dirichlet_ini:
                                        matrix_ApoC[j][j] -= DF #Convection will be dealt with further down
                                        matrix_ApoC[j][i] += DF
                    elif path == "membrane": #Membrane connection
                        if Apo_Contagion==2 and Sym_Contagion==2:
                            if PileUp==2:
                                iLayer=0
                                for iM in range(Nmaturity):
                                    if Barriers[iM]>0 and (G.node[j]['cgroup']==13 or G.node[j]['cgroup']==19 or G.node[j]['cgroup']==20 or ((j-NwallsJun in InterCid) and kInterC>0.0)) and (i not in list_ghostwalls): #
                                        #Diffusion between mature xylem vessels and their walls
                                        temp=1.0E-04*(lengths[i]*heights[iLayer:iLayer+int(Nlayers[iM])])/thickness #Section to length ratio (cm) for the xylem wall
                                        iK=0
                                        for iL in range(iLayer,int(iLayer+Nlayers[iM])):
                                            matrix_C2[iL*Ntot+i,iL*Ntot+i] -= temp[iK]*Diff_PW1
                                            matrix_C2[iL*Ntot+i,iL*Ntot+j] += temp[iK]*Diff_PW1
                                            matrix_C2[iL*Ntot+j,iL*Ntot+j] -= temp[iK]*Diff_PW1
                                            matrix_C2[iL*Ntot+j,iL*Ntot+i] += temp[iK]*Diff_PW1
                                            iK+=1
                                    elif D2O1==1 and (j-NwallsJun not in InterCid) and (i not in list_ghostwalls):
                                        #Diffusion of D2O across phospholipid bilayer
                                        temp=K_MB1*1.0E-04*(lengths[i]*heights[iLayer:iLayer+int(Nlayers[iM])])/dx_MB1 #Membrane area to thickness ratio (cm) times partition coefficient of membrane
                                        iK=0
                                        for iL in range(iLayer,int(iLayer+Nlayers[iM])):
                                            matrix_C2[iL*Ntot+i,iL*Ntot+i] -= temp[iK]*Diff_MB1
                                            matrix_C2[iL*Ntot+i,iL*Ntot+j] += temp[iK]*Diff_MB1
                                            matrix_C2[iL*Ntot+j,iL*Ntot+j] -= temp[iK]*Diff_MB1
                                            matrix_C2[iL*Ntot+j,iL*Ntot+i] += temp[iK]*Diff_MB1
                                            iK+=1
                                    for carrier in Active_transport_range:
                                        if int(carrier.get("tissue"))==G.node[j]['cgroup']:
                                            #Condition is that the protoplast (j) is an actual protoplast with membranes
                                            if j-NwallsJun not in InterCid and not (Barriers[iM]>0 and (G.node[j]['cgroup']==13 or G.node[j]['cgroup']==19 or G.node[j]['cgroup']==20)):
                                                temp=float(carrier.get("constant"))*(heights[iLayer:iLayer+int(Nlayers[iM])]+eattr['dist'])*eattr['length'] #Linear transport constant (Vmax/KM) [liter/day^-1/micron^-2] * membrane surface [micron²]
                                                iK=0
                                                if int(carrier.get("direction"))==1: #Influx transporter
                                                    for iL in range(iLayer,int(iLayer+Nlayers[iM])):
                                                        matrix_C2[iL*Ntot+j,iL*Ntot+i] += temp[iK] #Increase of concentration in protoplast (j) depends on concentration in cell wall (i)
                                                        matrix_C2[iL*Ntot+i,iL*Ntot+i] -= temp[iK] #Decrease of concentration in apoplast (i) depends on concentration in apoplast (i)
                                                        iK+=1
                                                elif int(carrier.get("direction"))==int(-1): #Efflux transporter
                                                    for iL in range(iLayer,int(iLayer+Nlayers[iM])):
                                                        matrix_C2[iL*Ntot+j,iL*Ntot+j] -= temp[iK]
                                                        matrix_C2[iL*Ntot+i,iL*Ntot+j] += temp[iK]
                                                        iK+=1
                                                else:
                                                    error('Error, carrier direction is either 1 (influx) or -1 (efflux), please correct in *_Hormones_Carriers_*.xml')
                                    iLayer+=int(Nlayers[iM])
                            else:
                                if Barrier>0 and (G.node[j]['cgroup']==13 or G.node[j]['cgroup']==19 or G.node[j]['cgroup']==20 or ((j-NwallsJun in InterCid) and kInterC>0.0)) and (i not in list_ghostwalls): #
                                    #Diffusion between mature xylem vessels and their walls
                                        temp=1.0E-04*(lengths[i]*height)/thickness #Section to length ratio (cm) for the xylem wall
                                        if sparseM==1:
                                            matrix_C[i,i] -= temp*Diff_PW1
                                            matrix_C[i,j] += temp*Diff_PW1
                                        else:
                                            matrix_C[i][i] -= temp*Diff_PW1
                                            matrix_C[i][j] += temp*Diff_PW1
                                        if sparseM==1:
                                            matrix_C[j,j] -= temp*Diff_PW1
                                            matrix_C[j,i] += temp*Diff_PW1
                                        else:
                                            matrix_C[j][j] -= temp*Diff_PW1
                                            matrix_C[j][i] += temp*Diff_PW1
                                elif D2O1==1 and (j-NwallsJun not in InterCid) and (i not in list_ghostwalls):
                                    #Diffusion of D2O across phospholipid bilayer
                                        temp=K_MB1*1.0E-04*(lengths[i]*height)/dx_MB1 #Membrane area to thickness ratio (cm) times partition coefficient of membrane
                                        if sparseM==1:
                                            matrix_C[i,i] -= temp*Diff_MB1
                                            matrix_C[i,j] += temp*Diff_MB1
                                        else:
                                            matrix_C[i][i] -= temp*Diff_MB1
                                            matrix_C[i][j] += temp*Diff_MB1
                                        if sparseM==1:
                                            matrix_C[j,j] -= temp*Diff_MB1
                                            matrix_C[j,i] += temp*Diff_MB1
                                        else:
                                            matrix_C[j][j] -= temp*Diff_MB1
                                            matrix_C[j][i] += temp*Diff_MB1
                                for carrier in Active_transport_range:
                                    if int(carrier.get("tissue"))==G.node[j]['cgroup']:
                                        #Condition is that the protoplast (j) is an actual protoplast with membranes
                                        if j-NwallsJun not in InterCid and not (Barrier>0 and (G.node[j]['cgroup']==13 or G.node[j]['cgroup']==19 or G.node[j]['cgroup']==20)):
                                            temp=float(carrier.get("constant"))*(height+eattr['dist'])*eattr['length'] #Linear transport constant (Vmax/KM) [liter/day^-1/micron^-2] * membrane surface [micron²]
                                            if int(carrier.get("direction"))==1: #Influx transporter
                                                    if sparseM==1:
                                                        matrix_C[j,i] += temp #Increase of concentration in protoplast (j) depends on concentration in cell wall (i)
                                                    else:
                                                        matrix_C[j][i] += temp #Increase of concentration in protoplast (j) depends on concentration in cell wall (i)
                                                    if sparseM==1:
                                                        matrix_C[i,i] -= temp #Decrease of concentration in apoplast (i) depends on concentration in apoplast (i)
                                                    else:
                                                        matrix_C[i][i] -= temp #Decrease of concentration in apoplast (i) depends on concentration in apoplast (i)
                                            elif int(carrier.get("direction"))==int(-1): #Efflux transporter
                                                    if sparseM==1:
                                                        matrix_C[j,j] -= temp
                                                    else:
                                                        matrix_C[j][j] -= temp #Increase of concentration in protoplast (j) depends on concentration in protoplast (j)
                                                    if sparseM==1:
                                                        matrix_C[i,j] += temp
                                                    else:
                                                        matrix_C[i][j] += temp #Decrease of concentration in apoplast (i) depends on concentration in protoplast (j)
                                            else:
                                                error('Error, carrier direction is either 1 (influx) or -1 (efflux), please correct in *_Hormones_Carriers_*.xml')
                        if G.node[j]['cgroup']==4: #Cortex (used to be elif)
                            if PileUp==2:
                                iLayer=0
                                for iM in range(Nmaturity):
                                    if Barriers[iM]>0 and (j-NwallsJun in InterCid): #Intercellular space "cell":
                                        kaqp=kInterC
                                    else:
                                        if not isnan(dist_grav[i]):
                                            kaqp=float(a_cortex*dist_grav[i]*1.0E-04+b_cortex) #AQP activity (cm/hPa/d) recalculated to allow for different kaqp on either sides of the same cortex cell
                                            if kaqp < 0:    
                                                error('Error, negative kaqp in cortical cell, adjust Paqp_cortex')
                                        else:
                                            kaqp=kaqp_cortex
                                    #Calculating each conductance
                                    if kwi[i,iM]==0.0 or kaqp==0.0: # kaqp==0 in case of air-filled intercellular space
                                        K[iLayer:iLayer+int(Nlayers[iM])] = 1.00E-16 #0?
                                    else:
                                        K[iLayer:iLayer+int(Nlayers[iM])] = 1/(1/(kwi[i,iM]/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(heights[iLayer]+eattr['dist'])*eattr['length']
                                    Kmb[jmb,iM]=K[iLayer]
                                    iLayer+=int(Nlayers[iM])
                                #if jmb<=10:
                                #    print(jmb,'Kmb',Kmb[jmb,iM],'wid',i,'cid',j-NwallsJun)
                                jmb+=1
                            else:
                                if Barrier>0 and (j-NwallsJun in InterCid):
                                    kaqp=kInterC
                                else:
                                    if not isnan(dist_grav[i]):
                                        kaqp=float(a_cortex*dist_grav[i]*1.0E-04+b_cortex) #AQP activity (cm/hPa/d) recalculated to allow for different kaqp on either sides of the same cortex cell
                                        if kaqp < 0:    
                                            error('Error, negative kaqp in cortical cell, adjust Paqp_cortex')
                                    else:
                                        kaqp=kaqp_cortex
                                #Calculating each conductance
                                if kwi[i]==0.0 or kaqp==0.0: # kaqp==0 in case of air-filled intercellular space
                                    K = 1.00E-16 #0?
                                else:
                                    K = 1/(1/(kwi[i]/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(height+eattr['dist'])*eattr['length']
                                Kmb[jmb]=K
                                #if jmb<=10:
                                #    print(jmb,'Kmb',K,'wid',i,'cid',j-NwallsJun)
                                jmb+=1
                        elif not G.node[j]['cgroup']==3: #Not endodermis
                            if PileUp==2:
                                iLayer=0
                                for iM in range(Nmaturity):
                                    kaqp=kaqpi[j-NwallsJun,iM]
                                    #Calculating each conductance
                                    if kwi[i,iM]==0.0 or kaqp==0.0: # kaqp==0 in case of air-filled intercellular space
                                        K[iLayer:iLayer+int(Nlayers[iM])] = 1.00E-16 #0?
                                    else:
                                        K[iLayer:iLayer+int(Nlayers[iM])] = 1/(1/(kwi[i,iM]/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(heights[iLayer]+eattr['dist'])*eattr['length']
                                    Kmb[jmb,iM]=K[iLayer]
                                    iLayer+=int(Nlayers[iM])
                                jmb+=1
                            else:
                                kaqp=kaqpi[j-NwallsJun]
                                #Calculating each conductance
                                if kwi[i]==0.0 or kaqp==0.0: # kaqp==0 in case of air-filled intercellular space
                                    K = 1.00E-16 #0?
                                else:
                                    K = 1/(1/(kwi[i]/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(height+eattr['dist'])*eattr['length']
                                Kmb[jmb]=K
                                #if jmb<=10:
                                #    print(jmb,'K init',K,'wid',i,'cid',j-NwallsJun)
                                jmb+=1
                        else: #Endodermis
                            if PileUp==2:
                                iLayer=0
                                for iM in range(Nmaturity):
                                    if count_stele_overall>0: #wall between endodermis and pericycle
                                        if count_passage>0:
                                            kw_temp = kws_passage[iM]
                                        else:
                                            kw_temp = kws_endo_peri[iM]
                                    elif count_endo>=2: #radial wall between endodermis cells
                                        kw_temp = min(kwi[i,iM],kws_endo_peri[iM])
                                    else:#elif count_stele_overall==0: #wall between endodermis and cortex
                                        if count_passage>0:
                                            kw_temp = kws_passage[iM]
                                        else:
                                            kw_temp = kws_endo_cortex[iM]
                                    kaqp=kaqpi[j-NwallsJun,iM]
                                    #Calculating each conductance
                                    if kw_temp==0.0 or kaqp==0.0: # kaqp==0 in case of air-filled intercellular space
                                        K[iLayer:iLayer+int(Nlayers[iM])] = 1.00E-16 #0?
                                    else:
                                        K[iLayer:iLayer+int(Nlayers[iM])] = 1/(1/(kw_temp/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(heights[iLayer]+eattr['dist'])*eattr['length']
                                    Kmb[jmb,iM]=K[iLayer]
                                    iLayer+=int(Nlayers[iM])
                                jmb+=1
                            else:
                                if count_stele_overall>0: #wall between endodermis and pericycle
                                    if count_passage>0:
                                        kw_temp = kw_passage
                                    else:
                                        kw_temp = kw_endo_peri
                                elif count_endo>=2: #radial wall between endodermis cells
                                    kw_temp = min(kwi[i],kw_endo_peri)
                                else:#elif count_stele_overall==0 and count_endo==1: #wall between endodermis and cortex
                                    if count_passage>0:
                                        kw_temp = kw_passage
                                    else:
                                        kw_temp = kw_endo_cortex 
                                kaqp=kaqpi[j-NwallsJun]
                                #Calculating each conductance
                                if kw_temp==0.0 or kaqp==0.0: # kaqp==0 in case of air-filled intercellular space
                                    K = 1.00E-16 #0?
                                else:
                                    K = 1/(1/(kw_temp/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(height+eattr['dist'])*eattr['length']
                                Kmb[jmb]=K
                                #if jmb<=10:
                                #    print(jmb,'K init',K,'wid',i,'cid',j-NwallsJun)
                                jmb+=1
                    elif path == "plasmodesmata": #Plasmodesmata connection
                        if PileUp==2:
                            if PD_freq_source==2:
                                #We have to find back the wall crossed by the PD from the 2 cells ID
                                wids=[]
                                for r in Cell2Wall_loop[i-NwallsJun]:
                                    wids.append(r.get("id"))
                                for r in Cell2Wall_loop[j-NwallsJun]:
                                    wid1=r.get("id")
                                    if wid1 in wids:
                                        break
                                temp_factor=Wall_PD[int(wid1)]*heights[:]*lengths[int(wid1)] #Total number of plasmodesmata on this wall (Wall_PD is now in PD per micron²)
                            else:
                                iLayer=0
                                temp_factor=ones((AxialLayers,1))#Quantity of plasmodesmata (adjusted by relative aperture)
                                for iM in range(Nmaturity):
                                    cgroupi=G.node[i]['cgroup']
                                    cgroupj=G.node[j]['cgroup']
                                    if cgroupi==19 or cgroupi==20:  #Xylem in new Cellset version
                                        cgroupi=13
                                    elif cgroupi==21: #Xylem Pole Pericyle in new Cellset version
                                        cgroupi=16
                                    elif cgroupi==23: #Phloem in new Cellset version
                                        cgroupi==11
                                    elif cgroupi==26 or cgroupi==24: #Companion Cell in new Cellset version
                                        cgroupi==12
                                    if cgroupj==19 or cgroupj==20:  #Xylem in new Cellset version
                                        cgroupj=13
                                    elif cgroupj==21: #Xylem Pole Pericyle in new Cellset version
                                        cgroupj=16
                                    elif cgroupj==23: #Phloem in new Cellset version
                                        cgroupj==11
                                    elif cgroupj==26 or cgroupj==24: #Companion Cell in new Cellset version
                                        cgroupj==12
                                    if ((j-NwallsJun in InterCid) or (i-NwallsJun in InterCid)) and Barriers[iM]>0: #one of the connected cells is an intercellular space "cell".
                                        temp_factor[iLayer:iLayer+int(Nlayers[iM])]=0.0
                                    elif Xylem_pieces and cgroupj==13 and cgroupi==13: #Fake wall splitting a xylem cell or vessel, high conductance in order to ensure homogeneous pressure within the splitted cell
                                        temp_factor[iLayer:iLayer+int(Nlayers[iM])]=10000*Fplxheight*1.0E-04*eattr['length'] #Quantity of PD
                                    elif Barriers[iM]>0 and (cgroupj==13 or cgroupi==13): #Mature xylem vessels, so no plasmodesmata with surrounding cells
                                        temp_factor[iLayer:iLayer+int(Nlayers[iM])]=0.0 #If Barrier==0, this case is treated like xylem is a stelar parenchyma cell
                                    elif (cgroupi==2 and cgroupj==1) or (cgroupj==2 and cgroupi==1):#Epidermis to exodermis cell or vice versa
                                        temp_factor[iLayer:iLayer+int(Nlayers[iM])]=Fplxheight_epi_exo*1.0E-04*eattr['length'] #Will not be used in case there is no exodermal layer
                                    elif (cgroupi==outercortex_connec_rank and cgroupj==4) or (cgroupj==outercortex_connec_rank and cgroupi==4):#Exodermis to cortex cell or vice versa
                                        temp=float(Kplrange[iPD].get("cortex_factor")) #Correction for specific cell-type PD aperture
                                        if Barriers[iM]>0:
                                            temp_factor[iLayer:iLayer+int(Nlayers[iM])]=2*temp/(temp+1)*Fplxheight_outer_cortex*1.0E-04*eattr['length']*Length_outer_cortex_tot /Length_outer_cortex_nospace
                                        else: #No aerenchyma
                                            temp_factor[iLayer:iLayer+int(Nlayers[iM])]=2*temp/(temp+1)*Fplxheight_outer_cortex*1.0E-04*eattr['length']
                                    elif (cgroupi==4 and cgroupj==4):#Cortex to cortex cell
                                        temp=float(Kplrange[iPD].get("cortex_factor")) #Correction for specific cell-type PD aperture
                                        if Barriers[iM]>0:
                                            temp_factor[iLayer:iLayer+int(Nlayers[iM])]=temp*Fplxheight_cortex_cortex*1.0E-04*eattr['length']*Length_cortex_cortex_tot /Length_cortex_cortex_nospace
                                        else: #No aerenchyma
                                            temp_factor[iLayer:iLayer+int(Nlayers[iM])]=temp*Fplxheight_cortex_cortex*1.0E-04*eattr['length']
                                    elif (cgroupi==3 and cgroupj==4) or (cgroupj==3 and cgroupi==4):#Cortex to endodermis cell or vice versa
                                        temp=float(Kplrange[iPD].get("cortex_factor")) #Correction for specific cell-type PD aperture
                                        if Barriers[iM]>0:
                                            temp_factor[iLayer:iLayer+int(Nlayers[iM])]=2*temp/(temp+1)*Fplxheight_cortex_endo*1.0E-04*eattr['length']*Length_cortex_endo_tot /Length_cortex_endo_nospace
                                        else: #No aerenchyma
                                            temp_factor[iLayer:iLayer+int(Nlayers[iM])]=2*temp/(temp+1)*Fplxheight_cortex_endo*1.0E-04*eattr['length']
                                    elif (cgroupi==3 and cgroupj==3):#Endodermis to endodermis cell
                                        temp_factor[iLayer:iLayer+int(Nlayers[iM])]=Fplxheight_endo_endo*1.0E-04*eattr['length']
                                    elif (cgroupi==3 and cgroupj==16) or (cgroupj==3 and cgroupi==16):#Pericycle to endodermis cell or vice versa
                                        if (i-NwallsJun in PPP) or (j-NwallsJun in PPP):
                                            temp=float(Kplrange[iPD].get("PPP_factor")) #Correction for specific cell-type PD aperture
                                        else:
                                            temp=1
                                        temp_factor[iLayer:iLayer+int(Nlayers[iM])]=2*temp/(temp+1)*Fplxheight_endo_peri*1.0E-04*eattr['length']
                                    elif (cgroupi==16 and (cgroupj==5 or cgroupj==13)) or (cgroupj==16 and (cgroupi==5 or cgroupi==13)):#Pericycle to stele cell or vice versa
                                        if (i-NwallsJun in PPP) or (j-NwallsJun in PPP):
                                            temp=float(Kplrange[iPD].get("PPP_factor")) #Correction for specific cell-type PD aperture
                                        else:
                                            temp=1
                                        temp_factor[iLayer:iLayer+int(Nlayers[iM])]=2*temp/(temp+1)*Fplxheight_peri_stele*1.0E-04*eattr['length']
                                    elif ((cgroupi==5 or cgroupi==13) and cgroupj==12) or (cgroupi==12 and (cgroupj==5 or cgroupj==13)):#Stele to companion cell
                                        temp=float(Kplrange[iPD].get("PCC_factor")) #Correction for specific cell-type PD aperture
                                        temp_factor[iLayer:iLayer+int(Nlayers[iM])]=2*temp/(temp+1)*Fplxheight_stele_comp*1.0E-04*eattr['length']
                                    elif (cgroupi==16 and cgroupj==12) or (cgroupi==12 and cgroupj==16):#Pericycle to companion cell
                                        temp1=float(Kplrange[iPD].get("PCC_factor"))
                                        if (i-NwallsJun in PPP) or (j-NwallsJun in PPP):
                                            temp2=float(Kplrange[iPD].get("PPP_factor")) #Correction for specific cell-type PD aperture
                                        else:
                                            temp2=1
                                        temp_factor[iLayer:iLayer+int(Nlayers[iM])]=2*temp1*temp2/(temp1+temp2)*Fplxheight_peri_comp*1.0E-04*eattr['length']
                                    elif (cgroupi==12 and cgroupj==12):#Companion to companion cell
                                        temp=float(Kplrange[iPD].get("PCC_factor"))
                                        temp_factor[iLayer:iLayer+int(Nlayers[iM])]=temp*Fplxheight_comp_comp*1.0E-04*eattr['length']
                                    elif (cgroupi==12 and cgroupj==11) or (cgroupi==11 and cgroupj==12):#Companion to phloem sieve tube cell
                                        temp=float(Kplrange[iPD].get("PCC_factor"))
                                        temp_factor[iLayer:iLayer+int(Nlayers[iM])]=2*temp/(temp+1)*Fplxheight_comp_sieve*1.0E-04*eattr['length']
                                    elif (cgroupi==16 and cgroupj==11) or (cgroupi==11 and cgroupj==16):#Pericycle to phloem sieve tube cell
                                        if (i-NwallsJun in PPP) or (j-NwallsJun in PPP):
                                            temp=float(Kplrange[iPD].get("PPP_factor")) #Correction for specific cell-type PD aperture
                                        else:
                                            temp=1
                                        temp_factor[iLayer:iLayer+int(Nlayers[iM])]=2*temp/(temp+1)*Fplxheight_peri_sieve*1.0E-04*eattr['length']
                                    elif ((cgroupi==5 or cgroupi==13) and cgroupj==11) or (cgroupi==11 and (cgroupj==5 or cgroupj==13)):#Stele to phloem sieve tube cell
                                        temp_factor[iLayer:iLayer+int(Nlayers[iM])]=Fplxheight_stele_sieve*1.0E-04*eattr['length']
                                    #elif cgroupi==13 and cgroupj==13: #Fake wall splitting a xylem cell or vessel, high conductance in order to ensure homogeneous pressure within the splitted cell
                                    #    temp_factor=10000*Fplxheight*1.0E-04*eattr['length']
                                    elif ((cgroupi==5 or cgroupi==13) and (cgroupj==5 or cgroupj==13)):#Stele to stele cell
                                        temp_factor[iLayer:iLayer+int(Nlayers[iM])]=Fplxheight_stele_stele*1.0E-04*eattr['length']
                                    else: #Default plasmodesmatal frequency
                                        temp_factor[iLayer:iLayer+int(Nlayers[iM])]=Fplxheight*1.0E-04*eattr['length'] #eattr['kpl']
                                    iLayer+=int(Nlayers[iM])
                            K = Kpl*temp_factor
                            ########Solute fluxes (diffusion across plasmodesmata)
                            if Sym_Contagion==2 and Apo_Contagion==2: #Sym & Apo contagion:
                                DF=PD_section*temp_factor/thickness*1.0E-04*Diff_PD1 #"Diffusive flux": Total PD cross-section area (micron^2) per unit PD length (micron) (tunred into cm) multiplied by solute diffusivity (cm^2/d) (yields cm^3/d)
                                iK=0
                                for iLayer in range(AxialLayers):
                                    matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] -= DF[iK]
                                    matrix_C2[iLayer*Ntot+i,iLayer*Ntot+j] += DF[iK]
                                    matrix_C2[iLayer*Ntot+j,iLayer*Ntot+j] -= DF[iK]
                                    matrix_C2[iLayer*Ntot+j,iLayer*Ntot+i] += DF[iK]
                                    iK+=1
                        else:
                            if PD_freq_source==2:
                                #We have to find back the wall crossed by the PD from the 2 cells ID
                                wids=[]
                                for r in Cell2Wall_loop[i-NwallsJun]:
                                    wids.append(r.get("id"))
                                for r in Cell2Wall_loop[j-NwallsJun]:
                                    wid1=r.get("id")
                                    if wid1 in wids:
                                        break
                                temp_factor=Wall_PD[int(wid1)]*height*lengths[int(wid1)] #Total number of plasmodesmata on this wall (Wall_PD is now in PD per micron²)
                            else:
                                cgroupi=G.node[i]['cgroup']
                                cgroupj=G.node[j]['cgroup']
                                if cgroupi==19 or cgroupi==20:  #Xylem in new Cellset version
                                    cgroupi=13
                                elif cgroupi==21: #Xylem Pole Pericyle in new Cellset version
                                    cgroupi=16
                                elif cgroupi==23: #Phloem in new Cellset version
                                    cgroupi==11
                                elif cgroupi==26 or cgroupi==24: #Companion Cell in new Cellset version
                                    cgroupi==12
                                if cgroupj==19 or cgroupj==20:  #Xylem in new Cellset version
                                    cgroupj=13
                                elif cgroupj==21: #Xylem Pole Pericyle in new Cellset version
                                    cgroupj=16
                                elif cgroupj==23: #Phloem in new Cellset version
                                    cgroupj==11
                                elif cgroupj==26 or cgroupj==24: #Companion Cell in new Cellset version
                                    cgroupj==12
                                temp_factor=1.0 #Quantity of plasmodesmata (adjusted by relative aperture)
                                if ((j-NwallsJun in InterCid) or (i-NwallsJun in InterCid)) and Barrier>0: #one of the connected cells is an intercellular space "cell".
                                    temp_factor=0.0
                                elif Xylem_pieces and cgroupj==13 and cgroupi==13: #Fake wall splitting a xylem cell or vessel, high conductance in order to ensure homogeneous pressure within the splitted cell
                                    temp_factor=10000*Fplxheight*1.0E-04*eattr['length'] #Quantity of PD
                                elif Barrier>0 and (cgroupj==13 or cgroupi==13): #Mature xylem vessels, so no plasmodesmata with surrounding cells
                                    temp_factor=0.0 #If Barrier==0, this case is treated like xylem is a stelar parenchyma cell
                                elif (cgroupi==2 and cgroupj==1) or (cgroupj==2 and cgroupi==1):#Epidermis to exodermis cell or vice versa
                                    temp_factor=Fplxheight_epi_exo*1.0E-04*eattr['length'] #Will not be used in case there is no exodermal layer
                                elif (cgroupi==outercortex_connec_rank and cgroupj==4) or (cgroupj==outercortex_connec_rank and cgroupi==4):#Exodermis to cortex cell or vice versa
                                    temp=float(Kplrange[iPD].get("cortex_factor")) #Correction for specific cell-type PD aperture
                                    if Barrier>0:
                                        temp_factor=2*temp/(temp+1)*Fplxheight_outer_cortex*1.0E-04*eattr['length']*Length_outer_cortex_tot /Length_outer_cortex_nospace
                                    else: #No aerenchyma
                                        temp_factor=2*temp/(temp+1)*Fplxheight_outer_cortex*1.0E-04*eattr['length']
                                elif (cgroupi==4 and cgroupj==4):#Cortex to cortex cell
                                    temp=float(Kplrange[iPD].get("cortex_factor")) #Correction for specific cell-type PD aperture
                                    if Barrier>0:
                                        temp_factor=temp*Fplxheight_cortex_cortex*1.0E-04*eattr['length']*Length_cortex_cortex_tot /Length_cortex_cortex_nospace
                                    else: #No aerenchyma
                                        temp_factor=temp*Fplxheight_cortex_cortex*1.0E-04*eattr['length']
                                elif (cgroupi==3 and cgroupj==4) or (cgroupj==3 and cgroupi==4):#Cortex to endodermis cell or vice versa
                                    temp=float(Kplrange[iPD].get("cortex_factor")) #Correction for specific cell-type PD aperture
                                    if Barrier>0:
                                        temp_factor=2*temp/(temp+1)*Fplxheight_cortex_endo*1.0E-04*eattr['length']*Length_cortex_endo_tot /Length_cortex_endo_nospace
                                    else: #No aerenchyma
                                        temp_factor=2*temp/(temp+1)*Fplxheight_cortex_endo*1.0E-04*eattr['length']
                                elif (cgroupi==3 and cgroupj==3):#Endodermis to endodermis cell
                                    temp_factor=Fplxheight_endo_endo*1.0E-04*eattr['length']
                                elif (cgroupi==3 and cgroupj==16) or (cgroupj==3 and cgroupi==16):#Pericycle to endodermis cell or vice versa
                                    if (i-NwallsJun in PPP) or (j-NwallsJun in PPP):
                                        temp=float(Kplrange[iPD].get("PPP_factor")) #Correction for specific cell-type PD aperture
                                    else:
                                        temp=1
                                    temp_factor=2*temp/(temp+1)*Fplxheight_endo_peri*1.0E-04*eattr['length']
                                elif (cgroupi==16 and (cgroupj==5 or cgroupj==13)) or (cgroupj==16 and (cgroupi==5 or cgroupi==13)):#Pericycle to stele cell or vice versa
                                    if (i-NwallsJun in PPP) or (j-NwallsJun in PPP):
                                        temp=float(Kplrange[iPD].get("PPP_factor")) #Correction for specific cell-type PD aperture
                                    else:
                                        temp=1
                                    temp_factor=2*temp/(temp+1)*Fplxheight_peri_stele*1.0E-04*eattr['length']
                                elif ((cgroupi==5 or cgroupi==13) and cgroupj==12) or (cgroupi==12 and (cgroupj==5 or cgroupj==13)):#Stele to companion cell
                                    temp=float(Kplrange[iPD].get("PCC_factor")) #Correction for specific cell-type PD aperture
                                    temp_factor=2*temp/(temp+1)*Fplxheight_stele_comp*1.0E-04*eattr['length']
                                elif (cgroupi==16 and cgroupj==12) or (cgroupi==12 and cgroupj==16):#Pericycle to companion cell
                                    temp1=float(Kplrange[iPD].get("PCC_factor"))
                                    if (i-NwallsJun in PPP) or (j-NwallsJun in PPP):
                                        temp2=float(Kplrange[iPD].get("PPP_factor")) #Correction for specific cell-type PD aperture
                                    else:
                                        temp2=1
                                    temp_factor=2*temp1*temp2/(temp1+temp2)*Fplxheight_peri_comp*1.0E-04*eattr['length']
                                elif (cgroupi==12 and cgroupj==12):#Companion to companion cell
                                    temp=float(Kplrange[iPD].get("PCC_factor"))
                                    temp_factor=temp*Fplxheight_comp_comp*1.0E-04*eattr['length']
                                elif (cgroupi==12 and cgroupj==11) or (cgroupi==11 and cgroupj==12):#Companion to phloem sieve tube cell
                                    temp=float(Kplrange[iPD].get("PCC_factor"))
                                    temp_factor=2*temp/(temp+1)*Fplxheight_comp_sieve*1.0E-04*eattr['length']
                                elif (cgroupi==16 and cgroupj==11) or (cgroupi==11 and cgroupj==16):#Pericycle to phloem sieve tube cell
                                    if (i-NwallsJun in PPP) or (j-NwallsJun in PPP):
                                        temp=float(Kplrange[iPD].get("PPP_factor")) #Correction for specific cell-type PD aperture
                                    else:
                                        temp=1
                                    temp_factor=2*temp/(temp+1)*Fplxheight_peri_sieve*1.0E-04*eattr['length']
                                elif ((cgroupi==5 or cgroupi==13) and cgroupj==11) or (cgroupi==11 and (cgroupj==5 or cgroupj==13)):#Stele to phloem sieve tube cell
                                    temp_factor=Fplxheight_stele_sieve*1.0E-04*eattr['length']
                                #elif cgroupi==13 and cgroupj==13: #Fake wall splitting a xylem cell or vessel, high conductance in order to ensure homogeneous pressure within the splitted cell
                                #    temp_factor=10000*Fplxheight*1.0E-04*eattr['length']
                                elif ((cgroupi==5 or cgroupi==13) and (cgroupj==5 or cgroupj==13)):#Stele to stele cell
                                    temp_factor=Fplxheight_stele_stele*1.0E-04*eattr['length']
                                else: #Default plasmodesmatal frequency
                                    temp_factor=Fplxheight*1.0E-04*eattr['length'] #eattr['kpl']
                            K = Kpl*temp_factor
                            ########Solute fluxes (diffusion across plasmodesmata)
                            if Sym_Contagion==2:
                                DF=PD_section*temp_factor/thickness*1.0E-04*Diff_PD1 #"Diffusive flux": Total PD cross-section area (micron^2) per unit PD length (micron) (tunred into cm) multiplied by solute diffusivity (cm^2/d) (yields cm^3/d)
                                if Apo_Contagion==2: #Sym & Apo contagion
                                        if sparseM==1:
                                            matrix_C[i,i] -= DF
                                            matrix_C[i,j] += DF
                                        else:
                                            matrix_C[i][i] -= DF
                                            matrix_C[i][j] += DF #Convection will be dealt with further down
                                        if sparseM==1:
                                            matrix_C[j,j] -= DF
                                            matrix_C[j,i] += DF
                                        else:
                                            matrix_C[j][j] -= DF
                                            matrix_C[j][i] += DF
                                else: #Only Sym contagion
                                    if i not in N_Dirichlet_ini:
                                        matrix_SymC[i-NwallsJun][i-NwallsJun] -= DF
                                        matrix_SymC[i-NwallsJun][j-NwallsJun] += DF
                                    if j not in N_Dirichlet_ini:
                                        matrix_SymC[j-NwallsJun][j-NwallsJun] -= DF #Convection will be dealt with further down
                                        matrix_SymC[j-NwallsJun][i-NwallsJun] += DF
                    if PileUp==2:
                        for iLayer in range(AxialLayers):
                            matrix_W2[iLayer*Ntot+i,iLayer*Ntot+i] -= K[iLayer] #Filling the Doussan matrix (symmetric)
                            matrix_W2[iLayer*Ntot+j,iLayer*Ntot+j] -= K[iLayer]
                            matrix_W2[iLayer*Ntot+i,iLayer*Ntot+j] += K[iLayer]
                            matrix_W2[iLayer*Ntot+j,iLayer*Ntot+i] += K[iLayer]
                    else:
                        if sparseM==1:
                            matrix_W[i,i] -= K #Filling the Doussan matrix (symmetric)
                            matrix_W[i,j] += K
                            matrix_W[j,i] += K
                            matrix_W[j,j] -= K
                        else:
                            matrix_W[i][i] -= K #Filling the Doussan matrix (symmetric)
                            matrix_W[i][j] += K
                            matrix_W[j][i] += K
                            matrix_W[j][j] -= K
                    
        if Apo_Contagion==2 and Sym_Contagion==2:
            if PileUp==2:
                    #Creating the stacked matrix_C2 and adding axial components of solute diffusivity
                    if AxialLayers>1:
                        ##The lines below cause memory errors when using very large matrices
                        #temp=len(DF_axial[:-Ntot])
                        #matrix_C2[:-Ntot,:-Ntot] -= spdiags(DF_axial[:-Ntot].T, 0, temp, temp).toarray() #A[1, 100:200] = A[0, :100]
                        #matrix_C2[:-Ntot,Ntot:] += spdiags(DF_axial[:-Ntot].T, 0, temp, temp).toarray() #Adding connections across layers
                        #matrix_C2[Ntot:,Ntot:] -= spdiags(DF_axial[:-Ntot].T, 0, temp, temp).toarray() #A[1, 100:200] = A[0, :100]
                        #matrix_C2[Ntot:,:-Ntot] += spdiags(DF_axial[:-Ntot].T, 0, temp, temp).toarray() #Adding connections across layers
                        ##This is expected not to cause such errors
                        matrix_C2[:Ntot,:Ntot] -= spdiags(DF_axial[:Ntot].T, 0, Ntot, Ntot).toarray() #A[1, 100:200] = A[0, :100]
                        matrix_C2[:Ntot,Ntot:2*Ntot] = spdiags(DF_axial[:Ntot].T, 0, Ntot, Ntot).toarray() #Adding connections across layers
                        #Central layers
                        for istack in range(1,AxialLayers-1):
                            matrix_C2[istack*Ntot:(istack+1)*Ntot,istack*Ntot:(istack+1)*Ntot] -= spdiags(DF_axial[(istack-1)*Ntot:istack*Ntot].T, 0, Ntot, Ntot).toarray()
                            matrix_C2[istack*Ntot:(istack+1)*Ntot,(istack-1)*Ntot:istack*Ntot] = spdiags(DF_axial[(istack-1)*Ntot:istack*Ntot].T, 0, Ntot, Ntot).toarray()
                            matrix_C2[istack*Ntot:(istack+1)*Ntot,istack*Ntot:(istack+1)*Ntot] -= spdiags(DF_axial[istack*Ntot:(istack+1)*Ntot].T, 0, Ntot, Ntot).toarray()
                            matrix_C2[istack*Ntot:(istack+1)*Ntot,(istack+1)*Ntot:(istack+2)*Ntot] = spdiags(DF_axial[istack*Ntot:(istack+1)*Ntot].T, 0, Ntot, Ntot).toarray()
                        istack+=1
                        matrix_C2[istack*Ntot:(istack+1)*Ntot,istack*Ntot:(istack+1)*Ntot] -= spdiags(DF_axial[(istack-1)*Ntot:istack*Ntot].T, 0, Ntot, Ntot).toarray()
                        matrix_C2[istack*Ntot:(istack+1)*Ntot,(istack-1)*Ntot:istack*Ntot] = spdiags(DF_axial[(istack-1)*Ntot:istack*Ntot].T, 0, Ntot, Ntot).toarray()
            else:
                if sparseM==1:
                    #Creating the stacked matrix_C2 and adding axial components of solute diffusivity
                    #Distal layer
                    matrix_C2[:Ntot,:Ntot] = matrix_C
                    if AxialLayers>1:
                        matrix_C2[:Ntot,:Ntot] -= spdiags(DF_axial.T, 0, Ntot, Ntot).toarray() #A[1, 100:200] = A[0, :100]
                        matrix_C2[:Ntot,Ntot:2*Ntot] = spdiags(DF_axial.T, 0, Ntot, Ntot).toarray() #Adding connections across layers
                        #Central layers
                        for istack in range(1,AxialLayers-1):
                            matrix_C2[istack*Ntot:(istack+1)*Ntot,istack*Ntot:(istack+1)*Ntot] = matrix_C
                            matrix_C2[istack*Ntot:(istack+1)*Ntot,istack*Ntot:(istack+1)*Ntot] -= spdiags(2.0*DF_axial.T, 0, Ntot, Ntot).toarray() #A[1, 100:200] = A[0, :100]
                            matrix_C2[istack*Ntot:(istack+1)*Ntot,(istack+1)*Ntot:(istack+2)*Ntot] = spdiags(DF_axial.T, 0, Ntot, Ntot).toarray() #Adding connections across layers
                            matrix_C2[istack*Ntot:(istack+1)*Ntot,(istack-1)*Ntot:istack*Ntot] = spdiags(DF_axial.T, 0, Ntot, Ntot).toarray()
                        #Proximal layer
                        matrix_C2[-Ntot:,-Ntot:] = matrix_C
                        matrix_C2[-Ntot:,-Ntot:] -= spdiags(DF_axial.T, 0, Ntot, Ntot).toarray()
                        matrix_C2[-Ntot:,-2*Ntot:-Ntot] = spdiags(DF_axial.T, 0, Ntot, Ntot).toarray()
        
        #Adding matrix components at soil-wall and wall-xylem connections & rhs terms
        rhs = np.zeros((AxialLayers*Ntot,1))
        rhs_s = np.zeros((AxialLayers*Ntot,1)) #Initializing the right-hand side matrix of soil pressure potentials
        rhs_x = np.zeros((AxialLayers*Ntot,1)) #Initializing the right-hand side matrix of xylem pressure potentials
        rhs_p = np.zeros((AxialLayers*Ntot,1)) #Initializing the right-hand side matrix of hydrostatic potentials for phloem BC
        
        #Adding matrix components at soil-wall connections
        if PileUp==2:
            for wid in Borderwall:
                iLayer=0
                for iM in range(Nmaturity):
                    if (position[wid][0]>=Xcontacts[iM]) or (Wall2Cell[wid][0]-NwallsJun in Contact): #Wall (not including junctions) connected to soil
                        K_int = kw*1.0E-04*(lengths[wid]/2*heights[iLayer:int(iLayer+Nlayers[iM])])/(thickness/2)
                        iK=0
                        for iL in range(iLayer,int(iLayer+Nlayers[iM])):
                            matrix_W2[iL*Ntot+wid,iL*Ntot+wid] -= K_int[iK] #Half the wall length is used here as the other half is attributed to the junction (Only for connection to soil)
                            iK+=1
                        rhs_s[iLayer*Ntot+wid:int(iLayer+Nlayers[iM])*Ntot:Ntot] = -K_int
                    iLayer+=int(Nlayers[iM])
        else:
            for wid in Borderwall:
                if (position[wid][0]>=Xcontact) or (Wall2Cell[wid][0]-NwallsJun in Contact): #Wall (not including junctions) connected to soil
                    temp=1.0E-04*(lengths[wid]/2*height)/(thickness/2)
                    K=kw*temp #Half the wall length is used here as the other half is attributed to the junction (Only for connection to soil)
                    if sparseM==1:
                        matrix_W[wid,wid] -= K #Doussan matrix
                        rhs_s[wid:AxialLayers*Ntot:Ntot] = -K 
                    else:
                        matrix_W[wid][wid] -= K #Doussan matrix
                        rhs_s[wid][0] = -K    #Right-hand side vector, could become Psi_soil[idwall], which could be a function of the horizontal position
                    
        #Adding matrix components at soil-junction connections
        if PileUp==2:
            for jid in Borderjunction:
                iLayer=0
                for iM in range(Nmaturity):
                    if (position[jid][0]>=Xcontacts[iM]) or (Junction2Wall2Cell[jid-Nwalls][0]-NwallsJun in Contact) or (Junction2Wall2Cell[jid-Nwalls][1]-NwallsJun in Contact) or (Junction2Wall2Cell[jid-Nwalls][2]-NwallsJun in Contact): #Junction connected to soil
                        K_int = kw*1.0E-04*(lengths[jid]*heights[iLayer:int(iLayer+Nlayers[iM])])/(thickness/2) #Doussan matrix
                        iK=0
                        for iL in range(iLayer,int(iLayer+Nlayers[iM])):
                            matrix_W2[iL*Ntot+jid,iL*Ntot+jid] -= K_int[iK]
                            iK+=1
                        rhs_s[iLayer*Ntot+jid:int(iLayer+Nlayers[iM])*Ntot:Ntot] = -K_int
                    iLayer+=int(Nlayers[iM])
        else:
            for jid in Borderjunction:
                if (position[jid][0]>=Xcontact) or (Junction2Wall2Cell[jid-Nwalls][0]-NwallsJun in Contact) or (Junction2Wall2Cell[jid-Nwalls][1]-NwallsJun in Contact) or (Junction2Wall2Cell[jid-Nwalls][2]-NwallsJun in Contact): #Junction connected to soil
                    temp=1.0E-04*(lengths[jid]*height)/(thickness/2)
                    K=kw*temp
                    if sparseM==1:
                        matrix_W[jid,jid] -= K #Doussan matrix
                        rhs_s[jid:AxialLayers*Ntot:Ntot] = -K
                    else:
                        matrix_W[jid][jid] -= K #Doussan matrix
                        rhs_s[jid][0] = -K    #Right-hand side vector, could become Psi_soil[idwall], which could be a function of the horizontal position
        
        if PileUp==2:
            if AxialLayers>1:
                ##The lines below cause memory errors when using very large matrices
                ##Construct the full matrix_W2
                #temp=len(K_axial[:-Ntot])
                #matrix_W2[:-Ntot,:-Ntot] -= spdiags(K_axial[:-Ntot].T, 0, temp, temp).toarray() #A[1, 100:200] = A[0, :100]
                #matrix_W2[:-Ntot,Ntot:] = spdiags(K_axial[:-Ntot].T, 0, temp, temp).toarray() #Adding connections across layers
                #matrix_W2[Ntot:,Ntot:] -= spdiags(K_axial[:-Ntot].T, 0, temp, temp).toarray() #A[1, 100:200] = A[0, :100]
                #matrix_W2[Ntot:,:-Ntot] = spdiags(K_axial[:-Ntot].T, 0, temp, temp).toarray() #Adding connections across layers
                ##This is expected not to cause such errors
                istack=0
                matrix_W2[:Ntot,:Ntot] -= spdiags(K_axial[:Ntot].T, 0, Ntot, Ntot).toarray() #A[1, 100:200] = A[0, :100]
                matrix_W2[:Ntot,Ntot:2*Ntot] = spdiags(K_axial[:Ntot].T, 0, Ntot, Ntot).toarray() #Adding connections across layers
                #Central layers
                for istack in range(1,AxialLayers-1):
                    matrix_W2[istack*Ntot:(istack+1)*Ntot,istack*Ntot:(istack+1)*Ntot] -= spdiags(K_axial[(istack-1)*Ntot:istack*Ntot].T, 0, Ntot, Ntot).toarray()
                    matrix_W2[istack*Ntot:(istack+1)*Ntot,(istack-1)*Ntot:istack*Ntot] = spdiags(K_axial[(istack-1)*Ntot:istack*Ntot].T, 0, Ntot, Ntot).toarray()
                    matrix_W2[istack*Ntot:(istack+1)*Ntot,istack*Ntot:(istack+1)*Ntot] -= spdiags(K_axial[istack*Ntot:(istack+1)*Ntot].T, 0, Ntot, Ntot).toarray()
                    matrix_W2[istack*Ntot:(istack+1)*Ntot,(istack+1)*Ntot:(istack+2)*Ntot] = spdiags(K_axial[istack*Ntot:(istack+1)*Ntot].T, 0, Ntot, Ntot).toarray()
                istack+=1
                matrix_W2[istack*Ntot:(istack+1)*Ntot,istack*Ntot:(istack+1)*Ntot] -= spdiags(K_axial[(istack-1)*Ntot:istack*Ntot].T, 0, Ntot, Ntot).toarray()
                matrix_W2[istack*Ntot:(istack+1)*Ntot,(istack-1)*Ntot:istack*Ntot] = spdiags(K_axial[(istack-1)*Ntot:istack*Ntot].T, 0, Ntot, Ntot).toarray()
        else:
            if AxialLayers>1:
                #Construct the full matrix_W2
                #Distal layer
                matrix_W2[:Ntot,:Ntot] = matrix_W
                matrix_W2[:Ntot,:Ntot] -= spdiags(K_axial.T, 0, Ntot, Ntot).toarray() #A[1, 100:200] = A[0, :100]
                matrix_W2[:Ntot,Ntot:2*Ntot] = spdiags(K_axial.T, 0, Ntot, Ntot).toarray() #Adding connections across layers
                #Central layers
                for istack in range(1,AxialLayers-1):
                    matrix_W2[istack*Ntot:(istack+1)*Ntot,istack*Ntot:(istack+1)*Ntot] = matrix_W
                    matrix_W2[istack*Ntot:(istack+1)*Ntot,istack*Ntot:(istack+1)*Ntot] -= spdiags(2.0*K_axial.T, 0, Ntot, Ntot).toarray() #A[1, 100:200] = A[0, :100]
                    matrix_W2[istack*Ntot:(istack+1)*Ntot,(istack+1)*Ntot:(istack+2)*Ntot] = spdiags(K_axial.T, 0, Ntot, Ntot).toarray() #Adding connections across layers
                    matrix_W2[istack*Ntot:(istack+1)*Ntot,(istack-1)*Ntot:istack*Ntot] = spdiags(K_axial.T, 0, Ntot, Ntot).toarray()
                #Proximal layer
                matrix_W2[-Ntot:,-Ntot:] = matrix_W
                matrix_W2[-Ntot:,-Ntot:] -= spdiags(K_axial.T, 0, Ntot, Ntot).toarray()
                matrix_W2[-Ntot:,-2*Ntot:-Ntot] = spdiags(K_axial.T, 0, Ntot, Ntot).toarray() #Adding connections across layers
            elif sparseM==1: #Single layered sparse matrix
                matrix_W2 = matrix_W
        
        #Creating connections to xylem & phloem BC elements for kr calculation (either xylem or phloem flow occurs depending on whether the segment is in the differentiation or elongation zone)
        XylOK=0
        if PileUp==2:
            if Barriers[-1]>0:
                XylOK=1
        else:
            if Barrier>0:
                XylOK=1
        if XylOK==1:
            if not isnan(Psi_xyl[1][iMaturity][count]): #Pressure xylem BC applied at proximal end only
                for cid in listxyl:
                    #K=K_axial[-Ntot+cid] #Axial conductance of xylem vessels
                    if sparseM==1:
                        matrix_W2[-Ntot+cid,-Ntot+cid] -= K_axial[-Ntot+cid] #Filling the Doussan matrix (symmetric)
                    else:
                        matrix_W[cid][cid] -= K_axial[cid] #Filling the Doussan matrix (symmetric)
                    rhs_x[-Ntot+cid][0] = -K_axial[-Ntot+cid]  #Axial conductance of xylem vessels (formerly Ntot+i)
                if not isnan(Psi_xyl[0][iMaturity][count]):
                    print('Distal xylem pressure BC not accounted for in kr estimation')
                rhs = rhs_s*Psi_soil[0][count] + rhs_x*Psi_xyl[1][iMaturity][count] #multiplication of rhs components delayed till this point so that rhs_s & rhs_x can be re-used to calculate Q
            elif not isnan(Flow_xyl[1][0][count]): #Xylem flux BC applied at proximal end and possibly distal end (positive flow towards the shoot)
                i=0
                for cid in listxyl:
                    rhs_x[-Ntot+cid][0] += Flow_xyl[1][i+1][count] #Flow boundary condition at proximal end of segment
                    i+=1
                if not isnan(Flow_xyl[0][0][count]):
                    i=0
                    for cid in listxyl:
                        rhs_x[cid][0] -= Flow_xyl[0][i+1][count] #Flow boundary condition at distal end of segment
                        i+=1    
                rhs = rhs_s*Psi_soil[0][count] + rhs_x #multiplication of rhs components delayed till this point so that rhs_s & rhs_x can be re-used to calculate Q
            else:
                rhs = rhs_s*Psi_soil[0][count]
        else:#in this case, the xylem BC is replaced by a protophloem BC, which is why there is no metaphloem in the loop
            if not isnan(Psi_sieve[1][iMaturity][count]): #Pressure phloem BC applied at proximal end only
                for cid in listprotosieve:
                    rhs_p[-Ntot+cid][0] = -K_axial[-Ntot+cid]  #Axial conductance of phloem sieve tube
                    if sparseM==1:
                        matrix_W2[-Ntot+cid,-Ntot+cid] -= K_axial[-Ntot+cid]
                    else:
                        matrix_W[cid][cid] -= K_axial[-Ntot+cid]
                if not isnan(Psi_sieve[0][iMaturity][count]):
                    print('Distal phloem pressure BC not accounted for in kr estimation')
                rhs = rhs_s*Psi_soil[0][count] + rhs_p*Psi_sieve[1][iMaturity][count] #multiplication of rhs components delayed till this point so that rhs_s & rhs_x can be re-used to calculate Q
            elif not isnan(Flow_sieve[1][0][count]): # (positive flow towards the shoot)
                i=1 
                for cid in listsieve:
                    if cid in listprotosieve:
                        rhs_p[-Ntot+cid][0] += Flow_sieve[1][i][count]
                    i+=1
                if not isnan(Flow_sieve[0][0][count]):
                    i=1
                    for cid in listsieve:
                        if cid in listprotosieve:
                            rhs_p[cid][0] -= Flow_sieve[0][i][count]
                        i+=1
                rhs = rhs_s*Psi_soil[0][count] + rhs_p #multiplication of rhs components delayed till this point so that rhs_s & rhs_x can be re-used to calculate Q
            else:
                rhs = rhs_s*Psi_soil[0][count]
        
        t01 = time.perf_counter()
        print(t01-t00, "seconds process time to fill matrix_W.2 & matrix_C.2")
        
        ##################################################
        ##Solve Doussan equation, results in soln matrix##
        ##################################################
        
        print("Solving water flow")
        t00 = time.perf_counter()  #Lasts about 35 sec in the 24mm root

        if sparseM==1:
            matrix_W20=matrix_W2
            soln = spsolve(matrix_W2.tocsr(),rhs)
        else:
            soln = np.linalg.solve(matrix_W,rhs) #Solving the equation to get potentials inside the network
            #Verification that computation was correct
            verif1=np.allclose(np.dot(matrix_W,soln),rhs)
        
        t01 = time.perf_counter()
        print(t01-t00, "seconds process time to solve water flow")
        
        #Removing xylem and phloem BC terms
        if XylOK==1:
            if not isnan(Psi_xyl[1][iMaturity][count]): #Pressure xylem BC
                if sparseM==1:
                    for cid in listxyl:
                        matrix_W2[-Ntot+cid,-Ntot+cid] += K_axial[-Ntot+cid]
                else:
                    for cid in listxyl:
                        matrix_W[cid][cid] += K_axial[cid]
        else:
            if not isnan(Psi_sieve[1][iMaturity][count]):
                if sparseM==1:
                    for cid in listprotosieve:
                        matrix_W2[-Ntot+cid,-Ntot+cid] += K_axial[-Ntot+cid]
                else:
                    for cid in listprotosieve:
                        matrix_W[cid][cid] += K_axial[-Ntot+cid]
        
        #Flow rates at interfaces
        Q_soil=[]
        for ind in Borderwall:
            for iLayer in range(AxialLayers):
                Q_soil.append(rhs_s[iLayer*Ntot+ind]*(soln[iLayer*Ntot+ind]-Psi_soil[0][count])) #(cm^3/d) Positive for water flowing into the root
        for ind in Borderjunction:
            for iLayer in range(AxialLayers):
                Q_soil.append(rhs_s[iLayer*Ntot+ind]*(soln[iLayer*Ntot+ind]-Psi_soil[0][count])) #(cm^3/d) Positive for water flowing into the root
        
        Q_xyl=zeros((2,Nxyl)) #Flow rate in xylem tubes (cm^3/d) positive towards the shoot
        Q_sieve=zeros((2,Nsieve)) #Flow rate in phloem sieve tubes (cm^3/d) positive towards the shoot
        if XylOK==1:
            if not isnan(Psi_xyl[1][iMaturity][count]):
                i=0
                for cid in listxyl:
                    Q_xyl[1][i]=K_axial[-Ntot+cid]*(soln[-Ntot+cid]-Psi_xyl[1][iMaturity][count])
                    if not isnan(Psi_xyl[0][iMaturity][count]):
                        print('No distal xylem flow when using xylem pressure BC at the round calculating kr')
                        #Q_xyl[0][i]=K_axial[cid]*(Psi_xyl[0][iMaturity][count]-soln[cid])
                    rank=int(Cell_rank[cid-NwallsJun])
                    row=int(rank2row[rank])
                    Q_xyl_layer[row][iMaturity][count] += -Q_xyl[1][i]+Q_xyl[0][i]
                    i+=1
            elif not isnan(Flow_xyl[1][0][count]):
                i=0
                for cid in listxyl:
                    Q_xyl[1][i]=Flow_xyl[1][i+1][count] #(cm^3/d) 
                    if not isnan(Flow_xyl[0][0][count]):
                        Q_xyl[0][i]=Flow_xyl[0][i+1][count]
                    rank=int(Cell_rank[cid-NwallsJun])
                    row=int(rank2row[rank])
                    Q_xyl_layer[row][iMaturity][count] += -Q_xyl[1][i]+Q_xyl[0][i]
                    i+=1
                Psi_xyl_prox = mean(soln[-Ntot+listxyl]-(sum(Q_xyl[1])/sum(K_axial[-Ntot+listxyl]))) #Check sign of Q_xyl
        else:#in this case, the xylem BC is replaced by a protophloem BC, which is why there is no metaphloem in the loop
            if not isnan(Psi_sieve[1][iMaturity][count]):
                i=0
                for cid in listsieve:
                    if cid in listprotosieve:
                        Q_sieve[1][i]=K_axial[-Ntot+cid]*(soln[-Ntot+cid]-Psi_sieve[1][iMaturity][count]) #(cm^3/d) Negative for water flowing into phloem tubes at the proximal end
                        #Q=rhs_p[-Ntot+cid]*(soln[-Ntot+cid]-Psi_sieve[1][iMaturity][count])
                        #Q_sieve.append(Q) #(cm^3/d) Negative for water flowing into phloem tubes at the proximal end
                        if not isnan(Psi_sieve[0][iMaturity][count]):
                            print('No distal phloem flow when using phloem pressure BC at the round calculating kr')
                            #Q_sieve[0][i]=K_axial[cid]*(Psi_sieve[0][iMaturity][count]-soln[cid])
                        rank=int(Cell_rank[cid-NwallsJun])
                        row=int(rank2row[rank])
                        Q_sieve_layer[row][iMaturity][count] += Q_sieve[1][i]
                    i+=1
            elif not isnan(Flow_sieve[1][0][count]):
                temp=0
                i=0
                for cid in listsieve: #(cm^3/d) Negative for water flowing towards the tip
                    if cid in listprotosieve:
                        Q_sieve[1][i] -= rhs_p[-Ntot+cid]
                        #Q=-rhs_p[-Ntot+cid] 
                        if AxialLayers>1:
                            Q_sieve[0][i] -= rhs_p[cid]
                            #Q-=rhs_p[cid]
                        #Q_sieve.append(Q) 
                        rank=int(Cell_rank[cid-NwallsJun])
                        row=int(rank2row[rank])
                        Q_sieve_layer[row][iMaturity][count] += -Q_sieve[1][i]+Q_sieve[0][i] # += Q
                        temp+=K_axial[-Ntot+cid]
                    i+=1
                Psi_sieve_prox = mean(soln[-Ntot+listprotosieve]-(sum(Q_sieve[1])/temp)) #Check sign of Q_sieve
        
        if PileUp==2:
            Totheight=sum(heights)
        else:
            Totheight=height*AxialLayers
        Q_tot[iMaturity][0]=sum(Q_soil) #Total flow rate at root surface
        SUF_Layer=empty((AxialLayers,1))
        nBorderwall=len(Borderwall)
        for iLayer in range(0,AxialLayers):
            SUF_Layer[iLayer][0]=(sum(Q_soil[iLayer:nBorderwall*AxialLayers:AxialLayers])+sum(Q_soil[nBorderwall*AxialLayers+iLayer::AxialLayers]))/Q_tot[iMaturity][0]
        if XylOK==1:
            if not isnan(Psi_xyl[1][iMaturity][0]):
                kr_tot[iMaturity][0]=Q_tot[iMaturity][0]/(Psi_soil[0][0]-Psi_xyl[1][iMaturity][0])/perimeter/Totheight/1.0E-04
            elif not isnan(Flow_xyl[1][0][count]):
                kr_tot[iMaturity][0]=Q_tot[iMaturity][0]/(Psi_soil[0][0]-Psi_xyl_prox)/perimeter/Totheight/1.0E-04
            else:
                print('Error: Scenario 0 should have proximal xylem pressure or flux boundary conditions beyond the elongation zone')
        else:
            if not isnan(Psi_sieve[1][iMaturity][0]):
                kr_tot[iMaturity][0]=Q_tot[iMaturity][0]/(Psi_soil[0][0]-Psi_sieve[1][iMaturity][0])/perimeter/Totheight/1.0E-04
            elif not isnan(Flow_sieve[1][0][count]):
                kr_tot[iMaturity][0]=Q_tot[iMaturity][0]/(Psi_soil[0][0]-Psi_sieve_prox)/perimeter/Totheight/1.0E-04
            else:
                print('Error: Scenario 0 should have proximal phloem pressure or flux boundary conditions in the elongation zone')
        print("Flow rates per unit root length: soil ",(Q_tot[iMaturity][0]/Totheight/1.0E-04),"cm^2/d, xylem ",((sum(Q_xyl[0])-sum(Q_xyl[1]))/Totheight/1.0E-04),"cm^2/d, phloem ",(sum(Q_sieve)/(Totheight)/1.0E-04),"cm^2/d")
        print("Flow rates: soil ",(Q_tot[iMaturity][0]),"cm^3/d, xylem ",((sum(Q_xyl[0])-sum(Q_xyl[1]))),"cm^3/d, phloem ",(sum(Q_sieve)),"cm^3/d")
        print("Mass balance error:",(Q_tot[iMaturity][0]+sum(Q_xyl[0])-sum(Q_xyl[1])+sum(Q_sieve[0])-sum(Q_sieve[1])),"cm^3/d")
        if PileUp==2:
            print("Radial conductivity:",kr_tot[iMaturity][0],"cm/hPa/d, Barriers:",Barriers,", stack height: ",Totheight," microns")
        else:
            print("Radial conductivity:",kr_tot[iMaturity][0],"cm/hPa/d, Barrier:",Barrier,", stack height: ",Totheight," microns")
        
        if count==iEquil_xyl:
            if XylOK==1 and isnan(Psi_xyl[1][iMaturity][0]):
                Psi_xyl[1][iMaturity][0]=0.0
                for cid in listxyl:
                    Psi_xyl[1][iMaturity][0]+=soln[-Ntot+cid]/Nxyl #Average of xylem water pressures
        if count==iEquil_sieve:
            if XylOK==0 and isnan(Psi_sieve[1][iMaturity][0]):
                Psi_sieve[1][iMaturity][0]=0.0
                for cid in listprotosieve:
                    Psi_sieve[1][iMaturity][0]+=soln[-Ntot+cid]/Nprotosieve #Average of protophloem water pressures
        
        #Calculation of standard transmembrane fractions
        jmb=0 #Index for membrane conductance vector
        for node, edges in G.adjacency_iter() : #adjacency_iter returns an iterator of (node, adjacency dict) tuples for all nodes. This is the fastest way to look at every edge. For directed graphs, only outgoing adjacencies are included.
            i = indice[node] #Node ID number
            psi = empty((AxialLayers))
            if i<Nwalls: #wall ID
                for iLayer in range(AxialLayers): #AxialLayers includes all layers even if there are several maturation zones
                    psi[iLayer] = float(soln[i+iLayer*Ntot])    #Node water potential
                #Here we count surrounding cell types in order to identify in which row of the endodermis or exodermis we are.
                count_endo=0 #total number of endodermis cells around the wall
                count_stele_overall=0 #total number of stelar cells around the wall
                count_exo=0 #total number of exodermis cells around the wall
                count_epi=0 #total number of epidermis cells around the wall
                #count_stele=0 #total number of epidermis cells around the wall
                count_cortex=0 #total number of epidermis cells around the wall
                count_passage=0 #total number of passage cells around the wall
                for neighbour, eattr in edges.items(): #Loop on connections (edges)
                    if eattr['path'] == 'membrane': #Wall connection
                        if any(passage_cell_ID==array((indice[neighbour])-NwallsJun)):
                            count_passage+=1
                        if G.node[neighbour]['cgroup']==3:#Endodermis
                            count_endo+=1
                        elif G.node[neighbour]['cgroup']>4:#Pericycle or stele
                            count_stele_overall+=1
                        elif G.node[neighbour]['cgroup']==1:#Exodermis
                            count_exo+=1
                        elif G.node[neighbour]['cgroup']==2:#Epidermis
                            count_epi+=1
                        elif G.node[neighbour]['cgroup']==4:#Cortex
                            count_cortex+=1
                    # if G.node[neighbour]['cgroup']==5:#Stele
                    #     count_stele+=1
                for neighbour, eattr in edges.items(): #Loop on connections (edges)
                    j = indice[neighbour] #Neighbouring node ID number
                    path = eattr['path'] #eattr is the edge attribute (i.e. connection type)
                    if path == "membrane": #Membrane connection
                        psin = empty((AxialLayers))
                        for iLayer in range(AxialLayers):
                            psin[iLayer] = float(soln[j+iLayer*Ntot]) #Neighbouring node water potential
                        if PileUp==2:
                            Flow = empty((AxialLayers))
                            iLayer=0
                            for iM in range(Nmaturity):
                                K=Kmb[jmb,iM]
                                Flow[iLayer:iLayer+int(Nlayers[iM])] = K * (psi[iLayer:iLayer+int(Nlayers[iM])] - psin[iLayer:iLayer+int(Nlayers[iM])]) #Note that this is only valid because we are in the scenario 0 with no osmotic potentials
                                iLayer+=int(Nlayers[iM])
                        else:
                            K=Kmb[jmb]
                            Flow = K * (psi - psin) #Note that this is only valid because we are in the scenario 0 with no osmotic potentials
                        jmb+=1
                        
                        #Flow densities calculation
                        #Macroscopic distributed parameter for transmembrane flow
                        #Discretization based on cell layers and apoplasmic barriers
                        rank = int(Cell_rank[j-NwallsJun])
                        row = int(rank2row[rank])
                        if rank == 1 and count_epi > 0: #Outer exodermis
                            row += 1
                        if rank == 3 and count_cortex > 0: #Outer endodermis
                            if any(passage_cell_ID==array(j-NwallsJun)):# and Barrier==2: Future/past passage cells separated like actual passage cells
                                row += 2
                            else:
                                row += 3
                        elif rank == 3 and count_stele_overall > 0: #Inner endodermis
                            if any(passage_cell_ID==array(j-NwallsJun)):# and Barrier==2: Future/past passage cells separated like actual passage cells
                                row += 1
                        
                        if PileUp==2:
                            iXylOK=sum(Nlayers[Barriers==0])
                        else:
                            if Barrier==0:
                                iXylOK=AxialLayers
                            else:
                                iXylOK=0
                        for iLayer in range(AxialLayers):
                            if ((j-NwallsJun not in InterCid) and (j not in listxyl)) or iLayer<iXylOK: #Not part of STF if crosses an intercellular space "membrane" or mature xylem "membrane" (that is no membrane though still labelled like one)
                                iL=TopLayer-AxialLayers+iLayer
                                if Flow[iLayer] > 0 :
                                    UptakeLayer_plus[row,iL,count] += Flow[iLayer] #grouping membrane flow rates in cell layers
                                else:
                                    UptakeLayer_minus[row,iL,count] += Flow[iLayer]
                                if not Q_tot[iMaturity][0]==0:
                                    if Flow[iLayer]/Q_tot[iMaturity][0] > 0 :
                                        STFlayer_plus[row,iL] += Flow[iLayer]/Q_tot[iMaturity][0] #Cell standard transmembrane fraction (positive)
                                        STFcell_plus[j-NwallsJun,iL] += Flow[iLayer]/Q_tot[iMaturity][0] #Cell standard transmembrane fraction (positive)
                                    else:
                                        STFlayer_minus[row,iL] += Flow[iLayer]/Q_tot[iMaturity][0] #Cell standard transmembrane fraction (negative)
                                        STFcell_minus[j-NwallsJun,iL] += Flow[iLayer]/Q_tot[iMaturity][0] #Cell standard transmembrane fraction (negative)
                                    STFmb[jmb-1,iL] = Flow[iLayer]/Q_tot[iMaturity][0]
                            elif (j-NwallsJun in InterCid):
                                iL=TopLayer-AxialLayers+iLayer
                                if Flow[iLayer] > 0 :
                                    UptakeLayer_plus_InterCid[row,iL,count] += Flow[iLayer] #grouping membrane flow rates in cell layers
                                else:
                                    UptakeLayer_minus_InterCid[row,iL,count] += Flow[iLayer]
        #error('Check UptakeLayer_plus & minus')
        if PileUp==1:
            #z_prox=0 #position of the proximal end of the root segment (initialization)
            if Ini_cond==2:
                Pile_cc_time_series = zeros((Nscenarios,Nxyl,int((float(Print_times[Nprint-1].get("value"))-t_resume)/dt_tracer)+Nprint+1),dtype=dp)
            else:
                Pile_cc_time_series = zeros((Nscenarios,Nxyl,int(float(Print_times[Nprint-1].get("value"))/dt_tracer)+Nprint+1),dtype=dp)
            for idt in range(size(Pile_cc_time_series,2)):
                for i in range(len(N_AxialNeumann_transi)):
                    Pile_cc_time_series[0][i][idt]=cc_AxialNeumann_transi[i][0] #The concentration at slice 0 is the one that is imposed at the distal end of the system
        
        
        for count in range(1,Nscenarios):
            
            print('Updating BC for scenario '+str(count))
            t00 = time.perf_counter()
            
            #Initializing the connectivity matrix including boundary conditions
            rhs = np.zeros((AxialLayers*Ntot,1))
            #rhs_x = np.zeros((AxialLayers*Ntot,1)) #Initializing the right-hand side matrix of xylem pressure potentials
            rhs_p = np.zeros((AxialLayers*Ntot,1)) #Initializing the right-hand side matrix of hydrostatic potentials for phloem BC
            rhs_e = np.zeros((AxialLayers*Ntot,1)) #Initializing the right-hand side matrix of cell elongation
            rhs_o = np.zeros((AxialLayers*Ntot,1)) #Initializing the right-hand side matrix of osmotic potentials
            if PileUp==2:
                s_cells = np.zeros((Ncells,Nmaturity))
                Os_cells = np.zeros((Ncells,Nmaturity)) #Initializing the cell osmotic potential vector
                s_membranes = np.zeros((Nmb,Nmaturity,1)) #Initializing the membrane reflection coefficient vector
                #Os_membranes = np.zeros((Nmb,Nmaturity,2)) #Initializing the osmotic potential storage side by side of membranes (0 for the wall, 1 for the protoplast)
                if Ana_Os:
                    Os_walls = np.zeros((Nwalls,AxialLayers)) #Initializing wall osmotic potentials
                else:
                    Os_walls = np.zeros((Nwalls,Nmaturity)) #Initializing wall osmotic potentials
            else:
                s_cells = np.zeros((Ncells,1))
                Os_cells = np.zeros((Ncells,1)) #Initializing the cell osmotic potential vector
                s_membranes = np.zeros((Nmb,1)) #Initializing the membrane reflection coefficient vector
                #Os_membranes = np.zeros((Nmb,2)) #Initializing the osmotic potential storage side by side of membranes (0 for the wall, 1 for the protoplast)
                if Ana_Os:
                    Os_walls = np.zeros((Nwalls,AxialLayers)) #Initializing wall osmotic potentials
                else:
                    Os_walls = np.zeros((Nwalls,1)) #Initializing the wall osmotic potential vector
            #rhs_s invariable between different scenarios but can vary for different hydraulic properties
            
            #Apoplastic & symplastic convective direction matrices initialization
            if Apo_Contagion==1 or Sym_Contagion==1:
                Cell_connec_flow=zeros((Ncells,14),dtype=int) #Flow direction across plasmodesmata, positive when entering the cell, negative otherwise
                Apo_connec_flow=zeros((NwallsJun,5),dtype=int) #Flow direction across cell walls, rows correspond to apoplastic nodes, and the listed nodes in each row receive convective flow from the row node
                nApo_connec_flow=zeros((NwallsJun,1),dtype=int)
            
            print('Scenario #'+str(count))
            
            #Reflection coefficients of membranes (undimensional)
            s_hetero[0][count]=int(Psi_cell_range[count].get("s_hetero")) #0:Uniform, 1: non-uniform, stele twice more permeable to solute, 2: non-uniform, cortex twice more permeable to solute
            s_factor[0][count]=float(Psi_cell_range[count].get("s_factor")) #(undimensional [0 -> 1]) multiplies all sigma values
            Elong_cell[0][count]=float(Elong_cell_range[count].get("midpoint_rate")) #Cell elongation rate (cm/d)
            Elong_cell_side_diff[0][count]=float(Elong_cell_range[count].get("side_rate_difference")) #Difference between cell elongation rates on the sides of the root in the EZ (cm/d)
            if s_hetero[0][count]==0:
                s_epi=s_factor[0][count]*1.0
                s_exo_epi=s_factor[0][count]*1.0
                s_exo_cortex=s_factor[0][count]*1.0
                s_cortex=s_factor[0][count]*1.0
                s_endo_cortex=s_factor[0][count]*1.0
                s_endo_peri=s_factor[0][count]*1.0
                s_peri=s_factor[0][count]*1.0
                s_stele=s_factor[0][count]*1.0
                s_comp=s_factor[0][count]*1.0
                s_sieve=s_factor[0][count]*1.0
            elif s_hetero[0][count]==1:
                s_epi=s_factor[0][count]*1.0
                s_exo_epi=s_factor[0][count]*1.0
                s_exo_cortex=s_factor[0][count]*1.0
                s_cortex=s_factor[0][count]*1.0
                s_endo_cortex=s_factor[0][count]*1.0
                s_endo_peri=s_factor[0][count]*0.5
                s_peri=s_factor[0][count]*0.5
                s_stele=s_factor[0][count]*0.5
                s_comp=s_factor[0][count]*0.5
                s_sieve=s_factor[0][count]*0.5
            elif s_hetero[0][count]==2:
                s_epi=s_factor[0][count]*0.5
                s_exo_epi=s_factor[0][count]*0.5
                s_exo_cortex=s_factor[0][count]*0.5
                s_cortex=s_factor[0][count]*0.5
                s_endo_cortex=s_factor[0][count]*0.5
                s_endo_peri=s_factor[0][count]*1.0
                s_peri=s_factor[0][count]*1.0
                s_stele=s_factor[0][count]*1.0
                s_comp=s_factor[0][count]*1.0
                s_sieve=s_factor[0][count]*1.0
            
            #Osmotic potentials (hPa)
            Os_hetero[0][count]=int(Psi_cell_range[count].get("Os_hetero")) #0:Uniform, 1: non-uniform no KNO3 treatment, 2: non-uniform with KNO3 treatment to help guttation
            Os_cortex[0][count]=float(Psi_cell_range[count].get("Os_cortex")) # Cortical cell osmotic potential (hPa)
            Os_sieve[0][count]=float(BC_sieve_range[count].get("osmotic"))
            if Os_hetero[0][count]==0:
                #Os_apo=-3000 #-0.3 MPa (Enns et al., 2000) applied stress
                #-0.80 MPa (Enns et al., 2000) concentration of cortical cells, no KNO3
                Os_epi=float(Os_cortex[0][count])
                Os_exo=float(Os_cortex[0][count])
                Os_c1=float(Os_cortex[0][count])
                Os_c2=float(Os_cortex[0][count])
                Os_c3=float(Os_cortex[0][count])
                Os_c4=float(Os_cortex[0][count])
                Os_c5=float(Os_cortex[0][count])
                Os_c6=float(Os_cortex[0][count])
                Os_c7=float(Os_cortex[0][count])
                Os_c8=float(Os_cortex[0][count])
                Os_endo=float(Os_cortex[0][count])
                Os_peri=float(Os_cortex[0][count])
                Os_stele=float(Os_cortex[0][count])
                Os_comp=(float(Os_sieve[0][count])+Os_cortex[0][count])/2 #Average phloem and parenchyma
                #Os_sieve=float(Os_cortex[0][count])
            elif Os_hetero[0][count]==1:
                Os_epi=-5000 #(Rygol et al. 1993) #float(Os_cortex[0][count]) #-0.80 MPa (Enns et al., 2000) concentration of cortical cells, no KNO3
                Os_exo=-5700 #(Rygol et al. 1993) #float(Os_cortex[0][count]) #-0.80 MPa (Enns et al., 2000) concentration of cortical cells, no KNO3
                Os_c1=-6400 #(Rygol et al. 1993)
                Os_c2=-7100 #(Rygol et al. 1993)
                Os_c3=-7800 #(Rygol et al. 1993)
                Os_c4=-8500 #(Rygol et al. 1993)
                Os_c5=-9000 #(Rygol et al. 1993)
                Os_c6=-9300 #(Rygol et al. 1993)
                Os_c7=-9000 #(Rygol et al. 1993)
                Os_c8=-8500 #(Rygol et al. 1993)
                Os_endo=-6200 #-0.62 MPa (Enns et al., 2000) concentration of endodermis cells, no KNO3
                Os_peri=-5000 #-0.50 MPa (Enns et al., 2000) concentration of pericycle cells, no KNO3
                Os_stele=-7400 #-0.74 MPa (Enns et al., 2000) concentration of xylem parenchyma cells, no KNO3
                Os_comp=(float(Os_sieve[0][count])-7400)/2 #Average phloem and parenchyma
                #Os_sieve=-14200 #-1.42 MPa (Pritchard, 1996) in barley phloem
            elif Os_hetero[0][count]==2:
                Os_epi=-11200 #(Rygol et al. 1993) #float(Os_cortex[0][count]) #-1.26 MPa (Enns et al., 2000) concentration of cortical cells, with KNO3
                Os_exo=-11500 #(Rygol et al. 1993) #float(Os_cortex[0][count]) #-1.26 MPa (Enns et al., 2000) concentration of cortical cells, with KNO3
                Os_c1=-11800 #(Rygol et al. 1993)
                Os_c2=-12100 #(Rygol et al. 1993)
                Os_c3=-12400 #(Rygol et al. 1993)
                Os_c4=-12700 #(Rygol et al. 1993)
                Os_c5=-12850 #(Rygol et al. 1993)
                Os_c6=-12950 #(Rygol et al. 1993)
                Os_c7=-12850 #(Rygol et al. 1993)
                Os_c8=-12700 #(Rygol et al. 1993)
                Os_endo=-10500 #-1.05 MPa (Enns et al., 2000) concentration of endodermis cells, with KNO3
                Os_peri=-9200 #-0.92 MPa (Enns et al., 2000) concentration of pericycle cells, with KNO3
                Os_stele=-12100 #-1.21 MPa (Enns et al., 2000) concentration of xylem parenchyma cells, with KNO3
                Os_comp=(float(Os_sieve[0][count])-12100)/2 #Average of phloem and parenchyma
                #Os_sieve=-14200 #-1.42 MPa (Pritchard, 1996) in barley phloem
            elif Os_hetero[0][count]==3:
                Os_epi=float(Os_cortex[0][count])
                Os_exo=float(Os_cortex[0][count])
                Os_c1=float(Os_cortex[0][count])
                Os_c2=float(Os_cortex[0][count])
                Os_c3=float(Os_cortex[0][count])
                Os_c4=float(Os_cortex[0][count])
                Os_c5=float(Os_cortex[0][count])
                Os_c6=float(Os_cortex[0][count])
                Os_c7=float(Os_cortex[0][count])
                Os_c8=float(Os_cortex[0][count])
                Os_endo=float((Os_cortex[0][count]-5000.0)/2.0)
                Os_peri=-5000.0 #Simple case with no stele pushing water out
                Os_stele=-5000.0
                Os_comp=(float(Os_sieve[0][count])-5000.0)/2 #Average phloem and parenchyma
                #Os_sieve=-5000.0
            
            jmb=0 #Index for membrane conductance vector
            for node, edges in G.adjacency_iter() : #adjacency_iter returns an iterator of (node, adjacency dict) tuples for all nodes. This is the fastest way to look at every edge. For directed graphs, only outgoing adjacencies are included.
                i=indice[node] #Node ID number
                #Here we count surrounding cell types in order to identify on which side of the endodermis or exodermis we are.
                count_endo=0 #total number of endodermis cells around the wall
                count_stele_overall=0 #total number of stelar cells around the wall
                count_exo=0 #total number of exodermis cells around the wall
                count_epi=0 #total number of epidermis cells around the wall
                count_cortex=0 #total number of cortical cells around the wall
                count_passage=0 #total number of passage cells around the wall
                if i<Nwalls: #wall ID
                    for neighbour, eattr in edges.items(): #Loop on connections (edges)
                        if eattr['path'] == 'membrane': #Wall connection
                            if any(passage_cell_ID==array((indice[neighbour])-NwallsJun)):
                                count_passage+=1
                            if G.node[neighbour]['cgroup']==3:#Endodermis
                                count_endo+=1
                            elif G.node[neighbour]['cgroup']>4:#Pericycle or stele
                                count_stele_overall+=1
                            elif G.node[neighbour]['cgroup']==4:#Cortex
                                count_cortex+=1
                            elif G.node[neighbour]['cgroup']==1:#Exodermis
                                count_exo+=1
                            elif G.node[neighbour]['cgroup']==2:#Epidermis
                                count_epi+=1
                for neighbour, eattr in edges.items(): #Loop on connections (edges)
                    j = (indice[neighbour]) #neighbouring node number
                    if j > i: #Only treating the information one way to save time
                        path = eattr['path'] #eattr is the edge attribute (i.e. connection type)
                        if path == "membrane": #Membrane connection
                            #Cell and wall osmotic potentials (cell types: 1=Exodermis;2=epidermis;3=endodermis;4=cortex;5=stele;16=pericycle)
                            rank=int(Cell_rank[int(j-NwallsJun)])
                            row=int(rank2row[rank])   
                            if rank==1:#Exodermis
                                #print('Caution, the analytical solution of radial solute distribution does not work with PileUp==2 yet')
                                Os_cells[j-NwallsJun,:]=Os_exo
                                if count_epi==1: #wall between exodermis and epidermis
                                    s_membranes[jmb,:]=s_exo_epi
                                elif count_epi==0: #wall between exodermis and cortex or between two exodermal cells
                                    s_membranes[jmb,:]=s_exo_cortex
                            elif rank==2:#Epidermis
                                Os_cells[j-NwallsJun,:]=Os_epi
                                s_membranes[jmb,:]=s_epi
                            elif rank==3:#Endodermis
                                Os_cells[j-NwallsJun,:]=Os_endo
                                if count_stele_overall==0 and count_cortex>0: #wall between endodermis and cortex or between two endodermal cells
                                    s_membranes[jmb,:]=s_endo_cortex
                                elif count_stele_overall>0 and count_endo>0: #wall between endodermis and pericycle
                                    s_membranes[jmb,:]=s_endo_peri
                                else: #Wall between endodermal cells
                                    s_membranes[jmb,:]=s_endo_peri
                            elif rank>=25 and rank<50:#Cortex
                                if r_discret[4]==1:
                                    Os_cells[j-NwallsJun,:]=Os_c8
                                elif row==row_outercortex-7:
                                    Os_cells[j-NwallsJun,:]=Os_c8
                                elif row==row_outercortex-6:
                                    Os_cells[j-NwallsJun,:]=Os_c7
                                elif row==row_outercortex-5:
                                    Os_cells[j-NwallsJun,:]=Os_c6
                                elif row==row_outercortex-4:
                                    Os_cells[j-NwallsJun,:]=Os_c5
                                elif row==row_outercortex-3:
                                    Os_cells[j-NwallsJun,:]=Os_c4
                                elif row==row_outercortex-2:
                                    Os_cells[j-NwallsJun,:]=Os_c3
                                elif row==row_outercortex-1:
                                    Os_cells[j-NwallsJun,:]=Os_c2
                                elif row==row_outercortex:
                                    Os_cells[j-NwallsJun,:]=Os_c1
                                else:
                                    Os_cells[j-NwallsJun,:]=Os_c8
                                s_membranes[jmb,:]=s_cortex
                                if j-NwallsJun in InterCid and XylOK>0: 
                                    if PileUp==2:
                                        Os_cells[j-NwallsJun,Barriers>0]=0 #Os_membranes[jmb,Barriers>0,1]=0 #The value here does not matter as s_membrane equals zero
                                        s_membranes[jmb,Barriers>0]=0
                                    else:
                                        Os_cells[j-NwallsJun,:]=0 #Os_membranes[jmb][1]=0
                                        s_membranes[jmb]=0
                            elif G.node[j]['cgroup']==5:#Stelar parenchyma
                                Os_cells[j-NwallsJun,:]=Os_stele
                                s_membranes[jmb,:]=s_stele
                            elif rank==16 or rank==17 or rank==18:#Pericycle
                                Os_cells[j-NwallsJun,:]=Os_peri
                                s_membranes[jmb,:]=s_peri
                            elif G.node[j]['cgroup']==11 or G.node[j]['cgroup']==23:#Phloem sieve tube cell
                                if not isnan(Os_sieve[0][count]):
                                    if PileUp==2:
                                        if j in listprotosieve:
                                            Os_cells[j-NwallsJun,:]=float(Os_sieve[0][count]) #Os_membranes[jmb,:,1]=float(Os_sieve[0][count])
                                        else:
                                            Os_cells[j-NwallsJun,Barriers>0]=float(Os_sieve[0][count]) #Os_membranes[jmb,Barriers>0,1]=float(Os_sieve[0][count])
                                            Os_cells[j-NwallsJun,Barriers==0]=Os_stele #Os_membranes[jmb,Barriers==0,1]=Os_stele
                                    else:
                                        if XylOK>0 or j in listprotosieve:
                                            Os_cells[j-NwallsJun,:]=float(Os_sieve[0][count]) #Os_membranes[jmb][1]=float(Os_sieve[0][count])
                                        else:
                                            Os_cells[j-NwallsJun,:]=Os_stele #Os_membranes[jmb][1]=Os_stele
                                else:
                                    Os_cells[j-NwallsJun,:]=Os_stele
                                s_membranes[jmb,:]=s_sieve
                            elif G.node[j]['cgroup']==12 or G.node[j]['cgroup']==26 or G.node[j]['cgroup']==24:#Companion cell
                                if not isnan(Os_sieve[0][count]):
                                    #if PileUp==2:
                                    Os_cells[j-NwallsJun,:]=Os_comp
                                else:
                                    Os_cells[j-NwallsJun,:]=Os_stele
                                s_membranes[jmb,:]=s_comp
                            elif G.node[j]['cgroup']==13 or G.node[j]['cgroup']==19 or G.node[j]['cgroup']==20:#Xylem cell or vessel
                                if PileUp==2:
                                    Os_cells[j-NwallsJun,Barriers==0]=Os_stele #Os_membranes[jmb,Barriers==0,1]=Os_stele
                                    Os_cells[j-NwallsJun,Barriers>0]=0.0 #Os_membranes[jmb,Barriers>0,1]=0 #The value here does not matter as s_membrane equals zero
                                    s_membranes[jmb,Barriers==0]=s_stele
                                    s_membranes[jmb,Barriers>0]=0.0
                                else:
                                    if XylOK==0:
                                        Os_cells[j-NwallsJun,:]=Os_stele #Os_membranes[jmb][1]=Os_stele
                                        s_membranes[jmb]=s_stele
                                    else:
                                        Os_cells[j-NwallsJun,:]=0.0 #Os_membranes[jmb][1]=0 #The value here does not matter as s_membrane equals zero
                                        s_membranes[jmb]=0.0
                            jmb+=1
            
            #Calculation of equivalent osmotic potential on each side of the endodermis, in the symplast and the apoplast
            #Estimate symplastic and apoplastic driving forces
            Os_sym_eq[iMaturity,:,count]=0.0
            Os_apo_eq[iMaturity,:,count]=0.0
            if PileUp==2 or AxialLayers>1:
                Os_sym_eq_lay[iMaturity*AxialLayers:(iMaturity+1)*AxialLayers,:,count]=0.0
                Os_apo_eq_lay[iMaturity*AxialLayers:(iMaturity+1)*AxialLayers,:,count]=0.0
            jmb=0 #Index for membrane vector
            for node, edges in G.adjacency_iter() : #adjacency_iter returns an iterator of (node, adjacency dict) tuples for all nodes. This is the fastest way to look at every edge. For directed graphs, only outgoing adjacencies are included.
                i = indice[node] #Node ID number
                if i<Nwalls: #wall ID
                    if r_rel[i]>=0: #cortical side
                        for neighbour, eattr in edges.items(): #Loop on connections (edges)
                            if eattr['path'] == 'membrane': #Wall connection
                                if PileUp==2 or AxialLayers>1:
                                    for iLayer in range(0,AxialLayers):
                                        temp1=STFmb[jmb,iLayer]*(Os_cells[neighbour-NwallsJun,int(MaturL[iLayer])]*s_membranes[jmb,int(MaturL[iLayer])]) #(Os_membranes[jmb,int(MaturL[iLayer]),1]*s_membranes[jmb,int(MaturL[iLayer])])
                                        temp2=STFmb[jmb,iLayer]*(Os_soil[0][count]*s_membranes[jmb,int(MaturL[iLayer])])
                                        Os_sym_eq_lay[iLayer+iMaturity*AxialLayers,0,count]+=temp1
                                        Os_apo_eq_lay[iLayer+iMaturity*AxialLayers,0,count]+=temp2
                                        Os_sym_eq[iMaturity,0,count]+=temp1
                                        Os_apo_eq[iMaturity,0,count]+=temp2
                                else:
                                    Os_sym_eq[iMaturity,0,count]+=STFmb[jmb,0]*(Os_cells[neighbour-NwallsJun,0]*s_membranes[jmb,0]) #(Os_membranes[jmb,int(MaturL[iLayer]),1]*s_membranes[jmb,int(MaturL[iLayer])])
                                    Os_apo_eq[iMaturity,0,count]+=STFmb[jmb,0]*(Os_soil[0][count]*s_membranes[jmb,0])
                                jmb+=1
                        Os_walls[i,:]=Os_soil[0][count]
                    else: #Stelar side
                        for neighbour, eattr in edges.items(): #Loop on connections (edges)
                            if eattr['path'] == 'membrane': #Wall connection
                                if PileUp==2 or AxialLayers>1:
                                    for iLayer in range(0,AxialLayers):
                                        temp1=STFmb[jmb,iLayer]*(Os_cells[neighbour-NwallsJun,int(MaturL[iLayer])]*s_membranes[jmb,int(MaturL[iLayer])]) #(Os_membranes[jmb,int(MaturL[iLayer]),1]*s_membranes[jmb,int(MaturL[iLayer])])
                                        if BarriersL[iLayer]==0:
                                            temp2=STFmb[jmb,iLayer]*(Os_soil[0][count]*s_membranes[jmb,int(MaturL[iLayer])])
                                        else:
                                            temp2=STFmb[jmb,iLayer]*(Os_xyl[0][count]*s_membranes[jmb,int(MaturL[iLayer])])
                                        Os_sym_eq_lay[iLayer+iMaturity*AxialLayers,1,count]-=temp1
                                        Os_apo_eq_lay[iLayer+iMaturity*AxialLayers,1,count]-=temp2
                                        Os_sym_eq[iMaturity,1,count]-=temp1
                                        Os_apo_eq[iMaturity,1,count]-=temp2
                                else:
                                    Os_sym_eq[iMaturity,1,count]-=STFmb[jmb,0]*(Os_cells[neighbour-NwallsJun,0]*s_membranes[jmb,0]) #(Os_membranes[jmb,int(MaturL[iLayer]),1]*s_membranes[jmb,int(MaturL[iLayer])])
                                    if Barrier==0:
                                        Os_apo_eq[iMaturity,1,count]-=STFmb[jmb,0]*(Os_soil[0][count]*s_membranes[jmb,0])
                                    else:
                                        Os_apo_eq[iMaturity,1,count]-=STFmb[jmb,0]*(Os_xyl[0][count]*s_membranes[jmb,0])
                                jmb+=1
                        if PileUp==2:
                            if Ana_Os:
                                Os_walls[i,np.where(BarriersL==0)]=Os_soil[0][count]
                                Os_walls[i,np.where(BarriersL>0)]=Os_xyl[0][count]
                            else:
                                Os_walls[i,np.where(Barriers==0)]=Os_soil[0][count]
                                Os_walls[i,np.where(Barriers>0)]=Os_xyl[0][count]
                        else:
                            Os_walls[i,np.where(Barrier==0)]=Os_soil[0][count]
                            Os_walls[i,np.where(Barrier>0)]=Os_xyl[0][count]
            if PileUp==2 or AxialLayers>1:
                Os_sym_eq_lay[iMaturity*AxialLayers:(iMaturity+1)*AxialLayers,0,count]/=SUF_Layer[:,0]
                Os_apo_eq_lay[iMaturity*AxialLayers:(iMaturity+1)*AxialLayers,0,count]/=SUF_Layer[:,0]
                Os_sym_eq_lay[iMaturity*AxialLayers:(iMaturity+1)*AxialLayers,1,count]/=SUF_Layer[:,0]
                Os_apo_eq_lay[iMaturity*AxialLayers:(iMaturity+1)*AxialLayers,1,count]/=SUF_Layer[:,0]
            
            #Soil and xylem water potentials
            #Psi_soil[0][count]=float(Psi_soil_range[count].get("pressure_left")) #Soil pressure potential (hPa)
            Psi_xyl[1][iMaturity][count]=float(BC_xyl_range[count].get("pressure_prox")) #Xylem pressure potential (hPa) at proximal end
            Psi_xyl[0][iMaturity][count]=float(BC_xyl_range[count].get("pressure_dist")) #Xylem pressure potential (hPa) at distal end
            dPsi_xyl[iMaturity][count]=float(BC_xyl_range[count].get("deltaP_prox")) #Xylem pressure potential change as compared to equilibrium pressure (hPa)
            Flow_xyl[1][0][count]=float(BC_xyl_range[count].get("flowrate_prox")) #Xylem flow rate on proximal end of the segment (cm^3/d)
            Flow_xyl[0][0][count]=float(BC_xyl_range[count].get("flowrate_dist")) #Xylem flow rate on distal end of the segment (cm^3/d)
            i_dist=0
            no_dist_BCx=False
            if not isnan(Psi_xyl[0][iMaturity][count]):
                i_dist+=1
            if not isnan(Flow_xyl[0][0][count]):
                i_dist+=1
            if i_dist==0:
                no_dist_BCx=True
            elif i_dist>1:
                error('May not have more than 1 water BC type at the distal end of xylem vessels')
            i_prox=0
            no_prox_BCx=False
            if not isnan(Psi_xyl[1][iMaturity][count]):
                i_prox+=1
            if not isnan(dPsi_xyl[iMaturity][count]):
                i_prox+=1
            if not isnan(Flow_xyl[1][0][count]):
                i_prox+=1
            if i_prox==0:
                no_prox_BCx=True
            elif i_prox>1:
                error('May not have more than 1 water BC type at the proximal end of xylem vessels')
            if not isnan(Flow_xyl[1][0][count]):
                tot_flow_prox=Flow_xyl[1][0][count]
                sum_area=0
                i=0
                for cid in listxyl:
                    area=cellarea[cid-NwallsJun]
                    Flow_xyl[1][i+1][count]=tot_flow_prox*area
                    sum_area+=area
                    i+=1
                i=0
                for cid in listxyl:
                    Flow_xyl[1][i+1][count]/=sum_area #Total xylem flow rate partitioned proportionnally to xylem cross-section area
                    i+=1
                if Flow_xyl[1][0][count]==0.0 and (Flow_xyl[0][0][count]==0.0 or isnan(Flow_xyl[0][0][count])):
                    iEquil_xyl=count
                if Ana_Os:
                    u=zeros((AxialLayers,2))
                    tot_flow_dist=0.0
                    if not isnan(Flow_xyl[0][0][count]):
                        tot_flow_dist=Flow_xyl[0][0][count]
                    if AxialLayers==1:
                        #Estimate the radial distribution of solutes later on from "u"
                        #First estimate water radial velocity in the apoplast
                        u[0][0]=(tot_flow_prox-tot_flow_dist)/(Totheight*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[0][0] #Cortex (cm/d)
                        u[0][1]=(tot_flow_prox-tot_flow_dist)/(Totheight*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[1][0] #Stele (cm/d)
                    else:
                        if PileUp==2:
                            temp=(tot_flow_prox-tot_flow_dist)*np.divide(SUF_Layer,heights*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[0][0]
                            u[:,0]=temp[:,0] #Cortex (cm/d)
                            temp=(tot_flow_prox-tot_flow_dist)*np.divide(SUF_Layer,heights*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[1][0]
                            u[:,1]=temp[:,0] #Stele (cm/d)
                        else:
                            temp=(tot_flow_prox-tot_flow_dist)*SUF_Layer/(heights*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[0][0]
                            u[:][0]=temp[:,0] #Cortex (cm/d)
                            temp=(tot_flow_prox-tot_flow_dist)*SUF_Layer/(heights*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[1][0]
                            u[:][1]=temp[:,0] #Stele (cm/d)
            elif not isnan(dPsi_xyl[iMaturity][count]):
                if not isnan(iEquil_xyl):
                    Psi_xyl[1][iMaturity][count]=Psi_xyl[1][iMaturity][iEquil_xyl]+dPsi_xyl[iMaturity][count]
                else:
                    print('Error: Cannot have xylem pressure change relative to equilibrium without having a prior scenario with equilibrium xylem boundary condition')
            if not isnan(Flow_xyl[0][0][count]):
                tot_flow_dist=Flow_xyl[0][0][count]
                sum_area=0
                i=0
                for cid in listxyl:
                    area=cellarea[cid-NwallsJun]
                    Flow_xyl[0][i+1][count]=tot_flow_dist*area
                    sum_area+=area
                    i+=1
                i=0
                for cid in listxyl: 
                    Flow_xyl[0][i+1][count]/=sum_area #Total xylem flow rate partitioned proportionnally to xylem cross-section area
                    i+=1
            
            if not isnan(Psi_xyl[1][iMaturity][count]):
                if Ana_Os:
                    if PileUp==2 or AxialLayers>1:
                        #Estimate the radial distribution of solutes
                        #First estimate total flow rate (cm^3/d) from BC & kr
                        u=zeros((AxialLayers,2))
                        tot_flow_dist=0.0
                        tot_flow_prox=0.0
                        #Express flow_ini for each Layer, assuming perfect apoplastic diffusion as a starting point
                        lay_flow_ini=kr_tot[iMaturity][0]*perimeter*Totheight*1.0E-04*SUF_Layer[:,0]*(Psi_soil[0][count]-Psi_xyl[1][iMaturity][count]+Os_apo_eq_lay[iMaturity*AxialLayers:(iMaturity+1)*AxialLayers,0,count]-Os_sym_eq_lay[iMaturity*AxialLayers:(iMaturity+1)*AxialLayers,0,count]-Os_apo_eq_lay[iMaturity*AxialLayers:(iMaturity+1)*AxialLayers,1,count]+Os_sym_eq_lay[iMaturity*AxialLayers:(iMaturity+1)*AxialLayers,1,count])
                        #Convergence loop of water radial velocity and solute apoplastic convection-diffusion
                        for iLayer in range(AxialLayers):
                            iter_n=0
                            lay_flow1=1E-18 # Adjusting tot_flow1 to avoid division by zero in first iteration
                            lay_flow2=lay_flow_ini[iLayer]
                            print('Layer',iLayer,'flow_rates pre & post =',lay_flow1,lay_flow2,'relative error',abs(lay_flow1-lay_flow2)/abs(lay_flow1),' iter_n =',iter_n)
                            lay_flow_best1=0.0
                            dlay_flow_best1=lay_flow2-lay_flow1
                            lay_flow_best2=NAN
                            dlay_flow_best2=NAN
                            stuck=False
                            while (abs(lay_flow1-lay_flow2)/max(abs(lay_flow1),abs(lay_flow2))>0.0002 and iter_n<50) or iter_n<3:
                                iter_n+=1
                                if stuck:
                                    lay_flow1=lay_flow_best1*0.9
                                    stuck=False
                                elif iter_n>1:
                                    lay_flow1=lay_flow_best1-dlay_flow_best1*(lay_flow_best1-lay_flow_best2)/(dlay_flow_best1-dlay_flow_best2)
                                else:
                                    lay_flow1=lay_flow2*0.5
                                u[iLayer,0]=lay_flow1/(Totheight*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[0][0] #Cortex apoplastic water velocity (cm/d) positive inwards
                                u[iLayer,1]=lay_flow1/(Totheight*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[1][0] #Stele apoplastic water velocity (cm/d) positive inwards
                                #Then estimate the radial solute distribution from an analytical solution (C(x)=C0+C0*(exp(u*x/D)-1)/(u/D*exp(u*x/D)-exp(u*L/D)+1)
                                Os_apo_eq_lay[iLayer+iMaturity*AxialLayers,0,count]=0.0
                                Os_apo_eq_lay[iLayer+iMaturity*AxialLayers,1,count]=0.0
                                jmb=0 #Index for membrane vector
                                for node, edges in G.adjacency_iter() : #adjacency_iter returns an iterator of (node, adjacency dict) tuples for all nodes. This is the fastest way to look at every edge. For directed graphs, only outgoing adjacencies are included.
                                    i = indice[node] #Node ID number
                                    if i<Nwalls: #wall ID
                                        for neighbour, eattr in edges.items(): #Loop on connections (edges)
                                            if eattr['path'] == 'membrane': #Wall connection
                                                if r_rel[i]>=0: #cortical side
                                                    Os_apo=float(Os_soil[0][count]*exp(u[iLayer,0]*abs(r_rel[i])*L_diff[0]/Os_soil[4][count]))
                                                    Os_apo_eq_lay[iLayer+iMaturity*AxialLayers,0,count]+=STFmb[jmb,iLayer]*(Os_apo*s_membranes[jmb,int(MaturL[iLayer])])
                                                    Os_apo_eq[iMaturity,0,count]+=STFmb[jmb,iLayer]*(Os_apo*s_membranes[jmb,int(MaturL[iLayer])])
                                                    Os_walls[i,iLayer]=Os_apo #Os_membranes[jmb][0]=Os_apo
                                                else: #Stelar side
                                                    Os_apo=float(Os_xyl[0][count]*exp(-u[iLayer,1]*abs(r_rel[i])*L_diff[1]/Os_xyl[4][count]))
                                                    Os_apo_eq_lay[iLayer+iMaturity*AxialLayers,1,count]-=STFmb[jmb,iLayer]*(Os_apo*s_membranes[jmb,int(MaturL[iLayer])])
                                                    Os_apo_eq[iMaturity,1,count]-=STFmb[jmb,iLayer]*(Os_apo*s_membranes[jmb,int(MaturL[iLayer])])
                                                    Os_walls[i,iLayer]=Os_apo #Os_membranes[jmb][0]=Os_apo
                                                jmb+=1
                                #Os_sym_eq_lay[iLayer,0,count]/=SUF_Layer[iLayer,0]
                                Os_apo_eq_lay[iLayer+iMaturity*AxialLayers,0,count]/=SUF_Layer[iLayer,0]
                                #Os_sym_eq_lay[iLayer,1,count]/=SUF_Layer[iLayer,0]
                                Os_apo_eq_lay[iLayer+iMaturity*AxialLayers,1,count]/=SUF_Layer[iLayer,0]
                                lay_flow2=float(kr_tot[iMaturity][0]*perimeter*Totheight*1.0E-04*SUF_Layer[iLayer,0]*(Psi_soil[0][count]-Psi_xyl[1][iMaturity][count]+Os_apo_eq_lay[iLayer+iMaturity*AxialLayers,0,count]-Os_sym_eq_lay[iLayer+iMaturity*AxialLayers,0,count]-Os_apo_eq_lay[iLayer+iMaturity*AxialLayers,1,count]+Os_sym_eq_lay[iLayer+iMaturity*AxialLayers,1,count]))
                                print('lay_flow1',lay_flow1,'lay_flow2',lay_flow2,'layer',iLayer)
                                dlay_flow=lay_flow2-lay_flow1
                                if iter_n==1:
                                    if abs(dlay_flow)>abs(dlay_flow_best1):
                                        lay_flow_best2=lay_flow1
                                        dlay_flow_best2=dlay_flow
                                    else:
                                        dlay_flow_best2=dlay_flow_best1
                                        lay_flow_best2=lay_flow_best1
                                        #lay_flow_criteria_best[iLayer]=float(abs(dlay_flow)/abs(lay_flow1[iLayer]))
                                        lay_flow_best1=lay_flow1
                                        dlay_flow_best1=dlay_flow
                                else:
                                    if abs(dlay_flow)<abs(dlay_flow_best1):
                                        dlay_flow_best2=dlay_flow_best1
                                        lay_flow_best2=lay_flow_best1
                                        dlay_flow_best1=dlay_flow
                                        lay_flow_best1=lay_flow1
                                        #lay_flow_criteria_best[iLayer]=float(abs(dlay_flow)/abs(lay_flow1[iLayer]))
                                    elif abs(dlay_flow)<abs(dlay_flow_best2):
                                        dlay_flow_best2=dlay_flow
                                        lay_flow_best2=lay_flow1
                                    else:
                                        stuck=True
                                print('Layer',iLayer,'flow_rates pre & post =',lay_flow1,lay_flow2,'relative error',abs(lay_flow1-lay_flow2)/max(abs(lay_flow1),abs(lay_flow2)),' iter_n =',iter_n)
                                if isnan(lay_flow2):
                                    error('NaN flow rate')
                            lay_flow_criteria_best=float(abs(dlay_flow_best1)/abs(lay_flow_best1))
                            if lay_flow_criteria_best<0.0002:
                                print('Success with convergence of solute radial distribution: relative error =',lay_flow_criteria_best,' flow rate =',lay_flow_best1)
                            else:
                                print('Closest to convergence of solute radial distribution: relative error =',lay_flow_criteria_best,' flow rate =',lay_flow_best1)
                            u[iLayer,0]=lay_flow_best1/(Totheight*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[0][0] #Cortex (cm/d)
                            u[iLayer,1]=lay_flow_best1/(Totheight*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[1][0] #Stele (cm/d)
                            tot_flow_prox+=lay_flow_best1
                            ##Then estimate osmotic potentials in radial walls later on: C(x)=C0+C0*(exp(u*x/D)-1)/(u/D*exp(u*x/D)-exp(u*L/D)+1)
                        Os_apo_eq[iMaturity,0,count]=0.0
                        Os_apo_eq[iMaturity,1,count]=0.0
                        for iLayer in range(AxialLayers):
                            jmb=0 #Index for membrane vector
                            for node, edges in G.adjacency_iter() : #adjacency_iter returns an iterator of (node, adjacency dict) tuples for all nodes. This is the fastest way to look at every edge. For directed graphs, only outgoing adjacencies are included.
                                i = indice[node] #Node ID number
                                if i<Nwalls: #wall ID
                                    for neighbour, eattr in edges.items(): #Loop on connections (edges)
                                        if eattr['path'] == 'membrane': #Wall connection
                                            if r_rel[i]>=0: #cortical side
                                                Os_apo_eq[iMaturity,0,count]+=STFmb[jmb,iLayer]*(Os_walls[i,iLayer]*s_membranes[jmb,int(MaturL[iLayer])])
                                            else: #Stelar side
                                                Os_apo_eq[iMaturity,1,count]-=STFmb[jmb,iLayer]*(Os_walls[i,iLayer]*s_membranes[jmb,int(MaturL[iLayer])])
                                            jmb+=1
                    else: #Single axial layer
                        #Estimate the radial distribution of solutes
                        #First estimate total flow rate (cm^3/d) from BC & kr
                        tot_flow1=0.0
                        u=zeros((2,1))
                        iter_n=0
                        tot_flow_ini=kr_tot[iMaturity][0]*perimeter*Totheight*1.0E-04*(Psi_soil[0][count]+Os_soil[0][count]-Os_sym_eq[iMaturity,0,count]-Psi_xyl[1][iMaturity][count]-Os_xyl[0][count]+Os_sym_eq[iMaturity,1,count]) 
                        tot_flow2=tot_flow_ini
                        print('flow_rate =',tot_flow_ini,' iter_n =',iter_n)
                        if tot_flow1==0.0:
                            tot_flow1=1E-15
                            print('Adjusting tot_flow1 to avoid division by zero in first iteration, tot_flow1 =',tot_flow1)
                        tot_flow_log=[]
                        tot_flow_best1=0
                        dtot_flow_best1=tot_flow_ini
                        tot_flow_best2=tot_flow_ini
                        dtot_flow_best2=nan
                        tot_flow_criteria_best=abs((tot_flow1-tot_flow2)/tot_flow1)
                        stuck=False
                        #Convergence loop of water radial velocity and solute apoplastic convection-diffusion
                        while (abs((tot_flow1-tot_flow2)/tot_flow1)>0.0002 and iter_n<50) or iter_n<3:
                            iter_n+=1
                            if stuck:
                                tot_flow1=tot_flow_best1*0.9
                                stuck=False
                            elif iter_n>1:
                                tot_flow1=tot_flow_best1-dtot_flow_best1*(tot_flow_best1-tot_flow_best2)/(dtot_flow_best1-dtot_flow_best2)
                            else:
                                tot_flow1=tot_flow2*0.1
                            #elif iter_n>1 and sign(tot_flow1/tot_flow2)==1:
                            #    tot_flow1=tot_flow1+(tot_flow2-tot_flow1)*min(0.5,max(0.01,10**(-abs((tot_flow1-tot_flow2)/tot_flow1))))#(tot_flow1+tot_flow2)/2
                            #else: #iter_n==2 and tot_flows have opposite signs
                            #    tot_flow1=tot_flow1*0.9
                            #Then estimate water radial velocity in the apoplast
                            u[0][0]=tot_flow1/(Totheight*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[0][0] #Cortex apoplastic water velocity (cm/d) positive inwards
                            u[1][0]=tot_flow1/(Totheight*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[1][0] #Stele apoplastic water velocity (cm/d) positive inwards
                            #Then estimate the radial solute distribution from an analytical solution (C(x)=C0+C0*(exp(u*x/D)-1)/(u/D*exp(u*x/D)-exp(u*L/D)+1)
                            Os_apo_eq[iMaturity,0,count]=0.0 #cortex side
                            Os_apo_eq[iMaturity,1,count]=0.0 #stele side
                            #temp1=0.0
                            #temp2=0.0
                            jmb=0 #Index for membrane vector
                            for node, edges in G.adjacency_iter() : #adjacency_iter returns an iterator of (node, adjacency dict) tuples for all nodes. This is the fastest way to look at every edge. For directed graphs, only outgoing adjacencies are included.
                                i = indice[node] #Node ID number
                                if i<Nwalls: #wall ID
                                    for neighbour, eattr in edges.items(): #Loop on connections (edges)
                                        if eattr['path'] == 'membrane': #Wall connection
                                            if r_rel[i]>=0: #cortical side
                                                Os_apo=Os_soil[0][count]*exp(u[0][0]*abs(r_rel[i])*L_diff[0]/Os_soil[4][count])
                                                Os_apo_eq[iMaturity,0,count]+=STFmb[jmb,iMaturity]*(Os_apo*s_membranes[jmb])
                                            else: #Stelar side
                                                Os_apo=Os_xyl[0][count]*exp(-u[1][0]*abs(r_rel[i])*L_diff[1]/Os_xyl[4][count])
                                                Os_apo_eq[iMaturity,1,count]-=STFmb[jmb,iMaturity]*(Os_apo*s_membranes[jmb])
                                            Os_walls[i,0]=Os_apo #Os_membranes[jmb][0]=Os_apo
                                            jmb+=1
                            tot_flow2=kr_tot[iMaturity][0]*perimeter*Totheight*1.0E-04*(Psi_soil[0][count]+Os_apo_eq[iMaturity,0,count]-Os_sym_eq[iMaturity,0,count]-Psi_xyl[1][iMaturity][count]-Os_apo_eq[iMaturity,1,count]+Os_sym_eq[iMaturity,1,count])
                            dtot_flow=tot_flow2-tot_flow1
                            if iter_n==1:
                                if abs(dtot_flow)>abs(dtot_flow_best1):
                                    dtot_flow_best2=dtot_flow
                                else:
                                    dtot_flow_best2=dtot_flow_best1
                                    dtot_flow_best1=dtot_flow
                                    tot_flow_best2=tot_flow_best1
                                    tot_flow_best1=tot_flow1
                                    tot_flow_criteria_best=abs(dtot_flow/tot_flow1)
                            else:
                                if abs(dtot_flow)<abs(dtot_flow_best1):
                                    dtot_flow_best2=dtot_flow_best1
                                    tot_flow_best2=tot_flow_best1
                                    dtot_flow_best1=dtot_flow
                                    tot_flow_best1=tot_flow1
                                    tot_flow_criteria_best=abs(dtot_flow/tot_flow1)
                                elif abs(dtot_flow)<abs(dtot_flow_best2):
                                    dtot_flow_best2=dtot_flow
                                    tot_flow_best2=tot_flow1
                                else:
                                    stuck=True
                            print('flow_rates pre & post =',tot_flow1,tot_flow2,'relative dif',abs(tot_flow1-tot_flow2)/abs(tot_flow1),' iter_n =',iter_n)
                            if isnan(tot_flow2):
                                error('NaN flow rate')
                            tot_flow_log.append((tot_flow1,tot_flow2,abs((tot_flow1-tot_flow2)/tot_flow1)))
                        if tot_flow_criteria_best<0.0002:
                            print('Success with convergence of solute radial distribution: relative error =',tot_flow_criteria_best,' flow rate =',tot_flow_best1)
                        else:
                            print('Closest to convergence of solute radial distribution: relative error =',tot_flow_criteria_best,' flow rate =',tot_flow_best1)
                        u[0][0]=tot_flow_best1/(Totheight*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[0][0] #Cortex (cm/d)
                        u[1][0]=tot_flow_best1/(Totheight*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[1][0] #Stele (cm/d)
                        ##Then estimate osmotic potentials in radial walls later on: C(x)=C0+C0*(exp(u*x/D)-1)/(u/D*exp(u*x/D)-exp(u*L/D)+1)
            
            #Elongation BC
            if PileUp==2:
                iLayer=0
                for iM in range(Nmaturity):
                    if Barriers[iM]==0: #No elongation from the Casparian strip on
                        for wid in range(Nwalls):
                            temp=lengths[wid]*thickness/2*1.0E-08*(Elong_cell[0][count]+(x_rel[wid]-0.5)*Elong_cell_side_diff[0][count])*Water_fraction_apo #cm^3/d Cell wall horizontal surface assumed to be rectangular (junctions are pointwise elements)
                            rhs_e[iLayer*Ntot+wid:int(iLayer+Nlayers[iM])*Ntot:Ntot][0]=temp
                        for cid in range(Ncells):
                            if cellarea[cid]>cellperimeter[cid]*thickness/2:
                                temp=(cellarea[cid]-cellperimeter[cid]*thickness/2)*1.0E-8*(Elong_cell[0][count]+(x_rel[NwallsJun+cid]-0.5)*Elong_cell_side_diff[0][count])*Water_fraction_sym #cm^3/d Wall thickness removed from cell horizontal area to obtain protoplast horizontal area
                                rhs_e[iLayer*Ntot+NwallsJun+cid:int(iLayer+Nlayers[iM])*Ntot:Ntot][0]=temp
                            else:
                                rhs_e[iLayer*Ntot+NwallsJun+cid:int(iLayer+Nlayers[iM])*Ntot:Ntot][0]=0 #The cell elongation virtually does not imply water influx, though its walls do (typically intercellular spaces)
                    iLayer+=int(Nlayers[iM])
            else:
                if Barrier==0: #No elongation from the Casparian strip on
                    for wid in range(Nwalls):
                        temp=lengths[wid]*thickness/2*1.0E-08*(Elong_cell[0][count]+(x_rel[wid]-0.5)*Elong_cell_side_diff[0][count])*Water_fraction_apo #cm^3/d Cell wall horizontal surface assumed to be rectangular (junctions are pointwise elements)
                        for iLayer in range(AxialLayers):
                            rhs_e[iLayer*Ntot+wid][0]=temp
                    for cid in range(Ncells):
                        if cellarea[cid]>cellperimeter[cid]*thickness/2:
                            temp=(cellarea[cid]-cellperimeter[cid]*thickness/2)*1.0E-8*(Elong_cell[0][count]+(x_rel[NwallsJun+cid]-0.5)*Elong_cell_side_diff[0][count])*Water_fraction_sym #cm^3/d Wall thickness removed from cell horizontal area to obtain protoplast horizontal area
                            for iLayer in range(AxialLayers):
                                rhs_e[iLayer*Ntot+NwallsJun+cid][0]=temp
                        else:
                            for iLayer in range(AxialLayers):
                                rhs_e[iLayer*Ntot+NwallsJun+cid][0]=0 #The cell elongation virtually does not imply water influx, though its walls do (typically intercellular spaces)
            
            #Phloem BC
            Psi_sieve[1][iMaturity][count]=float(BC_sieve_range[count].get("pressure_prox")) #Phloem sieve element pressure potential (hPa) at proximal end
            Psi_sieve[0][iMaturity][count]=float(BC_sieve_range[count].get("pressure_dist")) #Phloem sieve element pressure potential (hPa) at distal end
            dPsi_sieve[iMaturity][count]=float(BC_sieve_range[count].get("deltaP_prox")) #Phloem pressure potential change as compared to equilibrium pressure (hPa)
            Flow_sieve[1][0][count]=float(BC_sieve_range[count].get("flowrate_prox")) #Phloem flow rate (cm^3/d) at proximal end
            Flow_sieve[0][0][count]=float(BC_sieve_range[count].get("flowrate_dist")) #Phloem flow rate (cm^3/d) at distal end
            i_dist=0
            no_dist_BCp=False
            if not isnan(Psi_sieve[0][iMaturity][count]):
                i_dist+=1
            if not isnan(Flow_sieve[0][0][count]):
                i_dist+=1
            if i_dist==0:
                no_dist_BCp=True
            elif i_dist>1:
                error('May not have more than 1 water BC type at the distal end of xylem vessels')
            i_prox=0
            no_prox_BCp=False
            if not isnan(Psi_sieve[1][iMaturity][count]):
                i_prox+=1
            if not isnan(dPsi_sieve[iMaturity][count]):
                i_prox+=1
            if not isnan(Flow_sieve[1][0][count]):
                i_prox+=1
            if i_prox==0:
                no_prox_BCp=True
            elif i_prox>1:
                error('May not have more than 1 water BC type at the proximal end of xylem vessels')
            if not isnan(Flow_sieve[1][0][count]):
                if XylOK==0: #Essential difference is that there is only protophloem in EZ
                    if Flow_sieve[1][0][count]==0 and (Flow_sieve[0][0][count]==0 or isnan(Flow_sieve[0][0][count])): #This flag creates conditions with no soil water contributing to cell elongation (all the water then comes from the phloem)
                        tot_flow_prox=-float(sum(rhs_e)) #"Equilibrium condition" with phloem water fully used for elongation (no xylem in the EZ)
                    else:
                        tot_flow_prox=Flow_sieve[1][0][count]
                    sum_area=0
                    i=1
                    for cid in listsieve:
                        if cid in listprotosieve:
                            area=cellarea[cid-NwallsJun]
                            Flow_sieve[1][i][count]=tot_flow_prox*area
                            sum_area+=area
                        else:
                            Flow_sieve[1][i][count]=0
                        i+=1
                    i=1
                    for cid in listsieve:
                        if cid in listprotosieve:
                            Flow_sieve[1][i][count]/=sum_area #Total phloem flow rate partitioned proportionnally to phloem cross-section area
                        i+=1
                elif XylOK>0:
                    tot_flow_prox=Flow_sieve[1][0][count]
                    sum_area=0
                    i=1
                    for cid in listsieve:
                        area=cellarea[cid-NwallsJun]
                        Flow_sieve[1][i][count]=tot_flow_prox*area
                        sum_area+=area
                        i+=1
                    i=1
                    for cid in listsieve:
                        Flow_sieve[1][i][count]/=sum_area #Total phloem flow rate partitioned proportionnally to phloem cross-section area
                        i+=1
                if Flow_sieve[1][0][count]==0.0 and Flow_sieve[0][0][count]==0.0:
                    iEquil_sieve=count
            elif not isnan(dPsi_sieve[iMaturity][count]):
                if not isnan(iEquil_sieve):
                    Psi_sieve[1][iMaturity][count]=Psi_sieve[1][iMaturity][iEquil_sieve]+dPsi_sieve[iMaturity][count]
                else:
                    print('Error: Cannot have phloem pressure change relative to equilibrium without having a prior scenario with equilibrium phloem boundary condition')
            if not isnan(Flow_sieve[0][0][count]):
                if PileUp==2:
                    if Barriers[0]==0:
                        tot_flow_dist=Flow_sieve[0][0][count]
                        sum_area=0
                        i=1
                        for cid in listsieve:
                            if cid in listprotosieve:
                                area=cellarea[cid-NwallsJun]
                                Flow_sieve[0][i][count]=tot_flow_dist*area
                                sum_area+=area
                            else:
                                Flow_sieve[0][i][count]=0
                            i+=1
                        i=1
                        for cid in listsieve:
                            if cid in listprotosieve:
                                Flow_sieve[0][i][count]/=sum_area #Total phloem flow rate partitioned proportionnally to phloem cross-section area
                            i+=1
                    else:
                        tot_flow_dist=Flow_sieve[0][0][count]
                        sum_area=0
                        i=1
                        for cid in listsieve:
                            area=cellarea[cid-NwallsJun]
                            Flow_sieve[0][i][count]=tot_flow_dist*area
                            sum_area+=area
                            i+=1
                        i=1
                        for cid in listsieve:
                            Flow_sieve[0][i][count]/=sum_area #Total phloem flow rate partitioned proportionnally to phloem cross-section area
                            i+=1
                else:
                    if XylOK==0:
                        tot_flow_dist=Flow_sieve[0][0][count]
                        sum_area=0
                        i=1
                        for cid in listsieve:
                            if cid in listprotosieve:
                                area=cellarea[cid-NwallsJun]
                                Flow_sieve[0][i][count]=tot_flow_dist*area
                                sum_area+=area
                            else:
                                Flow_sieve[0][i][count]=0
                            i+=1
                        i=1
                        for cid in listsieve:
                            if cid in listprotosieve:
                                Flow_sieve[0][i][count]/=sum_area #Total phloem flow rate partitioned proportionnally to phloem cross-section area
                            i+=1
                    elif XylOK>0:
                        tot_flow_dist=Flow_sieve[0][0][count]
                        sum_area=0
                        i=1
                        for cid in listsieve:
                            area=cellarea[cid-NwallsJun]
                            Flow_sieve[0][i][count]=tot_flow_dist*area
                            sum_area+=area
                            i+=1
                        i=1
                        for cid in listsieve:
                            Flow_sieve[0][i][count]/=sum_area #Total phloem flow rate partitioned proportionnally to phloem cross-section area
                            i+=1
            
            jmb=0 #Index for membrane conductance vector
            for node, edges in G.adjacency_iter() : #adjacency_iter returns an iterator of (node, adjacency dict) tuples for all nodes. This is the fastest way to look at every edge. For directed graphs, only outgoing adjacencies are included.
                i=indice[node] #Node ID number
                #Here we count surrounding cell types in order to identify on which side of the endodermis or exodermis we are.
                count_endo=0 #total number of endodermis cells around the wall
                count_stele_overall=0 #total number of stelar cells around the wall
                count_exo=0 #total number of exodermis cells around the wall
                count_epi=0 #total number of epidermis cells around the wall
                count_cortex=0 #total number of cortical cells around the wall
                count_passage=0 #total number of passage cells around the wall
                if i<Nwalls: #wall ID
                    if Os_soil[2][count] == 1 or (Nmaturity==1 and int(Maturityrange[0].get("Barrier"))==0): #Left-right gradient for apoplastic osmotic potential
                        Os_soil_local=float(Os_soil[0][count]*(1-x_rel[i])+Os_soil[1][count]*x_rel[i])
                    elif Os_soil[2][count] == 2: #Central symmetrical gradient for apoplastic osmotic potential
                        if Os_soil[4][count] == 0: #Not the analytical solution
                            Os_soil_local=float(Os_soil[0][count]+(Os_soil[1][count]-Os_soil[0][count])*abs(r_rel[i])**Os_soil[3][count])
                        else:
                            if r_rel[i]>=0: #cortical side
                                if AxialLayers>1:
                                    Os_soil_local=Os_soil[0][count]*exp(u[:,0]*abs(r_rel[i])*L_diff[0]/Os_soil[4][count])
                                    Os_soil_local[np.where(BarriersL==0)]=ones(size(np.where(BarriersL==0)))*float(Os_soil[0][count]*(1-x_rel[i])+Os_soil[1][count]*x_rel[i])
                                else:
                                    Os_soil_local=Os_soil[0][count]*exp(u[0][0]*abs(r_rel[i])*L_diff[0]/Os_soil[4][count])
                                    Os_soil_local[np.where(Barrier==0)]=float(Os_soil[0][count]*(1-x_rel[i])+Os_soil[1][count]*x_rel[i])
                    if Os_xyl[2][count] == 2:
                        if Os_xyl[4][count] == 0: #Not the analytical solution
                            Os_xyl_local=float(Os_xyl[1][count]+(Os_xyl[0][count]-Os_xyl[1][count])*(1-abs(r_rel[i]))**Os_xyl[3][count])
                        else:
                            if r_rel[i]<0: #cortical side
                                if AxialLayers>1:
                                    Os_xyl_local=Os_xyl[0][count]*exp(-u[:,1]*abs(r_rel[i])*L_diff[1]/Os_xyl[4][count])
                                    Os_xyl_local[np.where(BarriersL==0)]=ones(size(np.where(BarriersL==0)))*float(Os_soil[0][count]*(1-x_rel[i])+Os_soil[1][count]*x_rel[i])
                                else:
                                    Os_xyl_local=Os_xyl[0][count]*exp(-u[1][0]*abs(r_rel[i])*L_diff[1]/Os_xyl[4][count])
                                    Os_xyl_local[np.where(Barrier==0)]=float(Os_soil[0][count]*(1-x_rel[i])+Os_soil[1][count]*x_rel[i])
                    elif Os_xyl[2][count] == 1:
                        Os_xyl_local=float((Os_xyl[0][count]+Os_xyl[1][count])/2)
                    for neighbour, eattr in edges.items(): #Loop on connections (edges)
                        if eattr['path'] == 'membrane': #Wall connection
                            if any(passage_cell_ID==array((indice[neighbour])-NwallsJun)):
                                count_passage+=1
                            if G.node[neighbour]['cgroup']==3:#Endodermis
                                count_endo+=1
                            elif G.node[neighbour]['cgroup']>4:#Pericycle or stele
                                count_stele_overall+=1
                            elif G.node[neighbour]['cgroup']==4:#Cortex
                                count_cortex+=1
                            elif G.node[neighbour]['cgroup']==1:#Exodermis
                                count_exo+=1
                            elif G.node[neighbour]['cgroup']==2:#Epidermis
                                count_epi+=1
                for neighbour, eattr in edges.items(): #Loop on connections (edges)
                    j = (indice[neighbour]) #neighbouring node number
                    if j > i: #Only treating the information one way to save time
                        path = eattr['path'] #eattr is the edge attribute (i.e. connection type)
                        if path == "membrane": #Membrane connection
                            #Cell and wall osmotic potentials (cell types: 1=Exodermis;2=epidermis;3=endodermis;4=cortex;5=stele;16=pericycle)
                            rank=int(Cell_rank[int(j-NwallsJun)])
                            row=int(rank2row[rank])
                            if rank==1:#Exodermis
                                if PileUp==2:
                                    #Os_membranes[jmb,:,1]=Os_exo
                                    OsCellLayer[row,:,count]+=Os_exo
                                    nOsCellLayer[row,:,count]+=1
                                    OsCellLayer[row+1,:,count]+=Os_exo
                                    nOsCellLayer[row+1,:,count]+=1
                                else:
                                    #Os_membranes[jmb][1]=Os_exo
                                    s_cells[j-NwallsJun] = s_epi
                                    OsCellLayer[row,iMaturity,count]+=Os_exo
                                    nOsCellLayer[row,iMaturity,count]+=1
                                    OsCellLayer[row+1,iMaturity,count]+=Os_exo
                                    nOsCellLayer[row+1,iMaturity,count]+=1
                                if Ana_Os and AxialLayers>1:
                                    Os_walls[i,:]=Os_soil_local[:] #may differ across layers
                                else:
                                    Os_walls[i,:]=Os_soil_local #may not differ across layers
                                Os_cells[j-NwallsJun,:]=Os_exo
                                if count_epi==1: #wall between exodermis and epidermis
                                    s_membranes[jmb,:]=s_exo_epi
                                    OsWallLayer[row+1][iMaturity][count]+=mean(Os_soil_local)
                                    nOsWallLayer[row+1][iMaturity][count]+=1
                                elif count_epi==0: #wall between exodermis and cortex or between two exodermal cells
                                    s_membranes[jmb,:]=s_exo_cortex
                                    OsWallLayer[row][iMaturity][count]+=mean(Os_soil_local)
                                    nOsWallLayer[row][iMaturity][count]+=1
                                s_cells[j-NwallsJun,:] = (s_exo_epi+s_exo_cortex)/2
                            elif rank==2:#Epidermis
                                Os_cells[j-NwallsJun,:]=Os_epi
                                if PileUp==2:
                                    #Os_membranes[jmb,:,1]=Os_epi
                                    OsCellLayer[row,:,count]+=Os_epi
                                    nOsCellLayer[row,:,count]+=1
                                else:
                                    #Os_membranes[jmb][1]=Os_epi
                                    OsCellLayer[row,iMaturity,count]+=Os_epi
                                    nOsCellLayer[row,iMaturity,count]+=1
                                OsWallLayer[row][iMaturity][count]+=mean(Os_soil_local)
                                nOsWallLayer[row][iMaturity][count]+=1
                                if Ana_Os and AxialLayers>1:
                                    Os_walls[i,:]=Os_soil_local[:] #may differ across layers
                                else:
                                    Os_walls[i,:]=Os_soil_local
                                s_membranes[jmb,:]=s_epi
                                s_cells[j-NwallsJun,:] = s_epi
                            elif rank==3:#Endodermis
                                Os_cells[j-NwallsJun,:]=Os_endo
                                if PileUp==2:
                                    #Os_membranes[jmb,:,1]=Os_endo
                                    OsCellLayer[row,:,count]+=Os_endo
                                    nOsCellLayer[row,:,count]+=1
                                    OsCellLayer[row+3,:,count]+=Os_endo
                                    nOsCellLayer[row+3,:,count]+=1
                                else:
                                    #Os_membranes[jmb][1]=Os_endo
                                    OsCellLayer[row,iMaturity,count]+=Os_endo
                                    nOsCellLayer[row,iMaturity,count]+=1
                                    OsCellLayer[row+3,iMaturity,count]+=Os_endo
                                    nOsCellLayer[row+3,iMaturity,count]+=1
                                if count_stele_overall==0 and count_cortex>0: #wall between endodermis and cortex or between two endodermal cells
                                    if Ana_Os and AxialLayers>1:
                                        Os_walls[i,:]=Os_soil_local[:] #may differ across layers
                                    else:
                                        Os_walls[i,:]=Os_soil_local
                                    s_membranes[jmb,:]=s_endo_cortex
                                    #Not including the osmotic potential of walls that are located at the same place as the casparian strip
                                    OsWallLayer[row+3][iMaturity][count]+=mean(Os_soil_local)
                                    nOsWallLayer[row+3][iMaturity][count]+=1
                                elif count_stele_overall>0 and count_endo>0: #wall between endodermis and pericycle
                                    if XylOK==0: #No apoplastic barrier
                                        if Ana_Os and AxialLayers>1:
                                            Os_walls[i,:]=Os_soil_local[:] #may differ across layers
                                        else:
                                            Os_walls[i,:]=Os_soil_local
                                        OsWallLayer[row][iMaturity][count]+=mean(Os_soil_local)
                                        nOsWallLayer[row][iMaturity][count]+=1
                                    else:
                                        if Ana_Os and AxialLayers>1:
                                            Os_walls[i,:]=Os_xyl_local[:] #may differ across layers
                                        else:
                                            Os_walls[i,0]=Os_xyl_local #float(Os_xyl[0][count])
                                        OsWallLayer[row][iMaturity][count]+=mean(Os_xyl_local) #float(Os_xyl[0][count])
                                        nOsWallLayer[row][iMaturity][count]+=1
                                    s_membranes[jmb,:]=s_endo_peri
                                else: #Wall between endodermal cells
                                    if XylOK==0: #No apoplastic barrier
                                        if Ana_Os and AxialLayers>1:
                                            Os_walls[i,:]=Os_soil_local[:] #may differ across layers
                                        else:
                                            Os_walls[i,:]=Os_soil_local
                                        OsWallLayer[row][iMaturity][count]+=mean(Os_soil_local)
                                        nOsWallLayer[row][iMaturity][count]+=1
                                    else:
                                        if Ana_Os and AxialLayers>1:
                                            Os_walls[i,:]=Os_xyl_local[:] #may differ across layers
                                        else:
                                            Os_walls[i,0]=Os_xyl_local #float(Os_xyl[0][count])
                                        OsWallLayer[row][iMaturity][count]+=mean(Os_xyl_local) #float(Os_xyl[0][count])
                                        nOsWallLayer[row][iMaturity][count]+=1
                                    s_membranes[jmb,:]=s_endo_peri
                                s_cells[j-NwallsJun,:] = (s_endo_peri+s_endo_cortex)/2
                            elif rank>=25 and rank<50:#Cortex
                                if r_discret[4]==1:
                                    Os_cells[j-NwallsJun,:]=Os_c8
                                    Ostemp=Os_c8
                                elif row==row_outercortex-7:
                                    Os_cells[j-NwallsJun,:]=Os_c8
                                    Ostemp=Os_c8
                                elif row==row_outercortex-6:
                                    Os_cells[j-NwallsJun,:]=Os_c7
                                    Ostemp=Os_c7
                                elif row==row_outercortex-5:
                                    Os_cells[j-NwallsJun,:]=Os_c6
                                    Ostemp=Os_c6
                                elif row==row_outercortex-4:
                                    Os_cells[j-NwallsJun,:]=Os_c5
                                    Ostemp=Os_c5
                                elif row==row_outercortex-3:
                                    Os_cells[j-NwallsJun,:]=Os_c4
                                    Ostemp=Os_c4
                                elif row==row_outercortex-2:
                                    Os_cells[j-NwallsJun,:]=Os_c3
                                    Ostemp=Os_c3
                                elif row==row_outercortex-1:
                                    Os_cells[j-NwallsJun,:]=Os_c2
                                    Ostemp=Os_c2
                                elif row==row_outercortex:
                                    Os_cells[j-NwallsJun,:]=Os_c1
                                    Ostemp=Os_c1
                                else:
                                    Os_cells[j-NwallsJun,:]=Os_c8
                                    Ostemp=Os_c8
                                if Ana_Os and AxialLayers>1:
                                    Os_walls[i,:]=Os_soil_local[:] #may differ across layers
                                else:
                                    Os_walls[i,:]=Os_soil_local
                                s_membranes[jmb,:]=s_cortex
                                s_cells[j-NwallsJun,:] = s_cortex
                                OsWallLayer[row][iMaturity][count]+=mean(Os_soil_local)
                                nOsWallLayer[row][iMaturity][count]+=1
                                if j-NwallsJun in InterCid and XylOK>0: 
                                    if PileUp==2:
                                        Os_cells[j-NwallsJun,Barriers>0]=Os_soil[0][count] #Os_soil_local[Barriers>0]
                                        if Ana_Os:
                                            Os_walls[i,BarriersL>0]=Os_soil_local[BarriersL>0] #may differ across layers
                                        else:
                                            Os_walls[i,:]=Os_soil_local
                                        #Os_membranes[jmb,Barriers>0,1]=Os_soil_local
                                        #Os_membranes[jmb,Barriers>0,0]=Os_soil_local
                                        s_membranes[jmb,Barriers>0]=0
                                        s_cells[j-NwallsJun,Barriers>0]=0
                                        OsCellLayer[row,Barriers==0,count]+=Ostemp
                                        nOsCellLayer[row,Barriers==0,count]+=1
                                        OsWallLayer[row,iMaturity,count]+=mean(Os_soil_local)
                                        nOsWallLayer[row,iMaturity,count]+=1
                                    else:
                                        Os_cells[j-NwallsJun,0]=Os_soil_local
                                        if Ana_Os and AxialLayers>1:
                                            Os_walls[i,:]=Os_soil_local[:] #may differ across layers
                                        else:
                                            Os_walls[i,:]=Os_soil_local
                                        #Os_membranes[jmb][1]=Os_soil_local
                                        #Os_membranes[jmb][0]=Os_soil_local
                                        s_membranes[jmb]=0
                                        s_cells[j-NwallsJun]=0
                                        OsCellLayer[row,iMaturity,count]+=Ostemp
                                        nOsCellLayer[row,iMaturity,count]+=1
                                        OsWallLayer[row][iMaturity][count]+=mean(Os_soil_local)
                                        nOsWallLayer[row][iMaturity][count]+=1
                                else:
                                    if PileUp==2:
                                        OsCellLayer[row,:,count]+=Ostemp
                                        nOsCellLayer[row,:,count]+=1
                                    else:
                                        OsCellLayer[row,iMaturity,count]+=Ostemp
                                        nOsCellLayer[row,iMaturity,count]+=1
                            elif G.node[j]['cgroup']==5:#Stelar parenchyma
                                Os_cells[j-NwallsJun,:]=Os_stele
                                if PileUp==2:
                                    #Os_membranes[jmb,:,1]=Os_stele
                                    OsCellLayer[row,:,count]+=Os_stele
                                    nOsCellLayer[row,:,count]+=1
                                else:
                                    #Os_membranes[jmb][1]=Os_stele
                                    OsCellLayer[row][iMaturity][count]+=Os_stele
                                    nOsCellLayer[row][iMaturity][count]+=1
                                if XylOK==0: #No apoplastic barrier
                                    if Ana_Os and AxialLayers>1:
                                        Os_walls[i,:]=Os_soil_local[:] #may differ across layers
                                    else:
                                        Os_walls[i,:]=Os_soil_local
                                    OsWallLayer[row][iMaturity][count]+=mean(Os_soil_local)
                                    nOsWallLayer[row][iMaturity][count]+=1
                                else:
                                    if Ana_Os and AxialLayers>1:
                                        Os_walls[i,:]=Os_xyl_local[:] #may differ across layers
                                    else:
                                        Os_walls[i,0]=Os_xyl_local #float(Os_xyl[0][count])
                                    OsWallLayer[row][iMaturity][count]+=mean(Os_xyl_local) #float(Os_xyl[0][count])
                                    nOsWallLayer[row][iMaturity][count]+=1
                                s_membranes[jmb,:]=s_stele
                                s_cells[j-NwallsJun,:] = s_stele
                            elif rank==16 or rank==17 or rank==18:#Pericycle
                                Os_cells[j-NwallsJun,:]=Os_peri
                                if PileUp==2:
                                    #Os_membranes[jmb,:,1]=Os_peri
                                    OsCellLayer[row,:,count]+=Os_peri
                                    nOsCellLayer[row,:,count]+=1
                                else:
                                    #Os_membranes[jmb][1]=Os_peri
                                    OsCellLayer[row][iMaturity][count]+=Os_peri
                                    nOsCellLayer[row][iMaturity][count]+=1
                                if XylOK==0: #No apoplastic barrier
                                    if Ana_Os and AxialLayers>1:
                                        Os_walls[i,:]=Os_soil_local[:] #may differ across layers
                                    else:
                                        Os_walls[i,:]=Os_soil_local
                                    OsWallLayer[row][iMaturity][count]+=Os_soil_local
                                    nOsWallLayer[row][iMaturity][count]+=1
                                else:
                                    if Ana_Os and AxialLayers>1:
                                        Os_walls[i,:]=Os_xyl_local[:] #may differ across layers
                                    else:
                                        Os_walls[i,0]=Os_xyl_local #float(Os_xyl[0][count])
                                    OsWallLayer[row][iMaturity][count]+=mean(Os_xyl_local) #float(Os_xyl[0][count])
                                    nOsWallLayer[row][iMaturity][count]+=1
                                s_membranes[jmb,:]=s_peri
                                s_cells[j-NwallsJun,:] = s_peri
                            elif G.node[j]['cgroup']==11 or G.node[j]['cgroup']==23:#Phloem sieve tube cell
                                if not isnan(Os_sieve[0][count]):
                                    if PileUp==2:
                                        if j in listprotosieve:
                                            Os_cells[j-NwallsJun,:]=float(Os_sieve[0][count])
                                            #Os_membranes[jmb,:,1]=float(Os_sieve[0][count])
                                            OsCellLayer[row,:,count]+=float(Os_sieve[0][count])
                                            nOsCellLayer[row,:,count]+=1
                                        else:
                                            Os_cells[j-NwallsJun,Barriers>0]=float(Os_sieve[0][count])
                                            Os_cells[j-NwallsJun,Barriers==0]=Os_stele
                                            #Os_membranes[jmb,Barriers>0,1]=float(Os_sieve[0][count])
                                            #Os_membranes[jmb,Barriers==0,1]=Os_stele
                                            OsCellLayer[row,Barriers>0,count]+=float(Os_sieve[0][count])
                                            OsCellLayer[row,Barriers==0,count]+=Os_stele
                                            nOsCellLayer[row,:,count]+=1
                                    else:
                                        if XylOK>0 or j in listprotosieve:
                                            Os_cells[j-NwallsJun,0]=float(Os_sieve[0][count])
                                            #Os_membranes[jmb][1]=float(Os_sieve[0][count])
                                            OsCellLayer[row][iMaturity][count]+=float(Os_sieve[0][count])
                                            nOsCellLayer[row][iMaturity][count]+=1
                                        else:
                                            Os_cells[j-NwallsJun,0]=Os_stele
                                            #Os_membranes[jmb][1]=Os_stele
                                            OsCellLayer[row][iMaturity][count]+=Os_stele
                                            nOsCellLayer[row][iMaturity][count]+=1
                                else:
                                    Os_cells[j-NwallsJun,:]=Os_stele
                                    if PileUp==2:
                                        #Os_membranes[jmb,:,1]=Os_stele
                                        OsCellLayer[row,:,count]+=Os_stele
                                        nOsCellLayer[row,:,count]+=1
                                    else:
                                        #Os_membranes[jmb][1]=Os_stele
                                        OsCellLayer[row][iMaturity][count]+=Os_stele
                                        nOsCellLayer[row][iMaturity][count]+=1
                                if XylOK==0: #No apoplastic barrier
                                    if Ana_Os and AxialLayers>1:
                                        Os_walls[i,:]=Os_soil_local[:] #may differ across layers
                                    else:
                                        Os_walls[i,:]=Os_soil_local
                                    OsWallLayer[row][iMaturity][count]+=mean(Os_soil_local)
                                    nOsWallLayer[row][iMaturity][count]+=1
                                else:
                                    if Ana_Os and AxialLayers>1:
                                        Os_walls[i,:]=Os_xyl_local[:] #may differ across layers
                                    else:
                                        Os_walls[i,0]=Os_xyl_local #float(Os_xyl[0][count])
                                    OsWallLayer[row][iMaturity][count]+=mean(Os_xyl_local) #float(Os_xyl[0][count])
                                    nOsWallLayer[row][iMaturity][count]+=1
                                s_membranes[jmb,:]=s_sieve
                                s_cells[j-NwallsJun,:] = s_sieve
                            elif G.node[j]['cgroup']==12 or G.node[j]['cgroup']==26 or G.node[j]['cgroup']==24:#Companion cell
                                if not isnan(Os_sieve[0][count]):
                                    Os_cells[j-NwallsJun,:]=Os_comp
                                    if PileUp==2:
                                        #Os_membranes[jmb,:,1]=Os_comp
                                        OsCellLayer[row,:,count]+=Os_comp
                                        nOsCellLayer[row,:,count]+=1
                                    else:
                                        #Os_membranes[jmb][1]=Os_comp
                                        OsCellLayer[row][iMaturity][count]+=Os_comp
                                        nOsCellLayer[row][iMaturity][count]+=1
                                else:
                                    Os_cells[j-NwallsJun,:]=Os_stele
                                    if PileUp==2:
                                        #Os_membranes[jmb,:,1]=Os_stele
                                        OsCellLayer[row,:,count]+=Os_stele
                                        nOsCellLayer[row,:,count]+=1
                                    else:
                                        #Os_membranes[jmb][1]=Os_stele
                                        OsCellLayer[row][iMaturity][count]+=Os_stele
                                        nOsCellLayer[row][iMaturity][count]+=1
                                if XylOK==0: #No apoplastic barrier
                                    if Ana_Os and AxialLayers>1:
                                        Os_walls[i,:]=Os_soil_local[:] #may differ across layers
                                    else:
                                        Os_walls[i,:]=Os_soil_local
                                    OsWallLayer[row][iMaturity][count]+=mean(Os_soil_local)
                                    nOsWallLayer[row][iMaturity][count]+=1
                                else:
                                    if Ana_Os and AxialLayers>1:
                                        Os_walls[i,:]=Os_xyl_local[:] #may differ across layers
                                    else:
                                        Os_walls[i,0]=Os_xyl_local #float(Os_xyl[0][count])
                                    OsWallLayer[row][iMaturity][count]+=mean(Os_xyl_local)
                                    nOsWallLayer[row][iMaturity][count]+=1
                                s_membranes[jmb,:]=s_comp
                                s_cells[j-NwallsJun,:] = s_comp
                            elif G.node[j]['cgroup']==13 or G.node[j]['cgroup']==19 or G.node[j]['cgroup']==20:#Xylem cell or vessel
                                if PileUp==2:
                                    Os_cells[j-NwallsJun,Barriers==0]=Os_stele
                                    Os_cells[j-NwallsJun,Barriers>0]=float(Os_xyl[0][count])#Os_xyl_local
                                    #Os_membranes[jmb,Barriers==0,1]=Os_stele
                                    #Os_membranes[jmb,Barriers>0,0]=Os_xyl_local
                                    #Os_membranes[jmb,Barriers>0,1]=Os_xyl_local
                                    OsCellLayer[row,Barriers==0,count]+=Os_stele
                                    nOsCellLayer[row,Barriers==0,count]+=1
                                    s_membranes[jmb,Barriers==0]=s_stele
                                    s_membranes[jmb,Barriers>0]=0.0
                                    s_cells[j-NwallsJun,Barriers==0] = s_stele
                                    s_cells[j-NwallsJun,Barriers>0]=0.0
                                    if XylOK==0:
                                        if Ana_Os and AxialLayers>1:
                                            Os_walls[i,:]=Os_soil_local[:] #may differ across layers
                                        else:
                                            Os_walls[i,:]=Os_soil_local
                                        OsWallLayer[row][iMaturity][count]+=mean(Os_soil_local)
                                        nOsWallLayer[row][iMaturity][count]+=1
                                    else:
                                        if Ana_Os and AxialLayers>1:
                                            Os_walls[i,BarriersL>0]=Os_xyl_local[BarriersL>0] #may differ across layers
                                            Os_walls[i,BarriersL==0]=Os_soil_local[BarriersL==0] #may differ across layers
                                        else:
                                            Os_walls[i,0]=Os_xyl_local #float(Os_xyl[0][count])
                                        OsWallLayer[row][iMaturity][count]+=mean(Os_xyl_local) #float(Os_xyl[0][count])
                                        nOsWallLayer[row][iMaturity][count]+=1
                                else:
                                    if XylOK==0:
                                        Os_cells[j-NwallsJun,0]=Os_stele
                                        #Os_membranes[jmb][1]=Os_stele
                                        OsCellLayer[row][iMaturity][count]+=Os_stele
                                        nOsCellLayer[row][iMaturity][count]+=1
                                        if Ana_Os and AxialLayers>1:
                                            Os_walls[i,:]=Os_soil_local[:] #may differ across layers
                                        else:
                                            Os_walls[i,:]=Os_soil_local
                                        s_membranes[jmb]=s_stele
                                        s_cells[j-NwallsJun] = s_stele
                                        OsWallLayer[row][iMaturity][count]+=mean(Os_soil_local)
                                        nOsWallLayer[row][iMaturity][count]+=1
                                    else:
                                        Os_cells[j-NwallsJun,0]=Os_xyl_local
                                        #Os_membranes[jmb][0]=Os_xyl_local
                                        #Os_membranes[jmb][1]=Os_xyl_local
                                        #Os_membranes[jmb][1]=Os_xyl_local
                                        if Ana_Os and AxialLayers>1:
                                            Os_walls[i,:]=Os_xyl_local[:] #may differ across layers
                                        else:
                                            Os_walls[i,0]=Os_xyl_local #float(Os_xyl[0][count])
                                        s_membranes[jmb]=0.0
                                        s_cells[j-NwallsJun]=0.0
                                        OsWallLayer[row][iMaturity][count]+=mean(Os_xyl_local) #float(Os_xyl[0][count])
                                        nOsWallLayer[row][iMaturity][count]+=1
                            if PileUp==2:
                                if Ana_Os and AxialLayers>1:
                                    iLayer=0
                                    #for iM in range(Nmaturity):
                                    for iLayer in range(AxialLayers):
                                        iM=int(MaturL[iLayer])
                                        rhs_o[iLayer*Ntot+i]+= Kmb[jmb,iM]*s_membranes[jmb,iM]*(Os_walls[i,iLayer] - Os_cells[j-NwallsJun,iM]) #Wall node
                                        rhs_o[iLayer*Ntot+j]+= Kmb[jmb,iM]*s_membranes[jmb,iM]*(Os_cells[j-NwallsJun,iM] - Os_walls[i,iLayer]) #Cell node 
                                        #rhs_o[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot]+= Kmb[jmb,iM]*s_membranes[jmb,iM]*(Os_walls[i,iLayer:int(iLayer+Nlayers[iM])] - Os_cells[j-NwallsJun,iM]) #Wall node
                                        #rhs_o[iLayer*Ntot+j:int(iLayer+Nlayers[iM])*Ntot:Ntot]+= Kmb[jmb,iM]*s_membranes[jmb,iM]*(Os_cells[j-NwallsJun,iM] - Os_walls[i,iLayer:int(iLayer+Nlayers[iM])]) #Cell node 
                                        #iLayer+=int(Nlayers[iM])
                                else:
                                    iLayer=0
                                    for iM in range(Nmaturity):
                                        rhs_o[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot]+= Kmb[jmb,iM]*s_membranes[jmb,iM]*(Os_walls[i,iM] - Os_cells[j-NwallsJun,iM]) #Wall node
                                        rhs_o[iLayer*Ntot+j:int(iLayer+Nlayers[iM])*Ntot:Ntot]+= Kmb[jmb,iM]*s_membranes[jmb,iM]*(Os_cells[j-NwallsJun,iM] - Os_walls[i,iM]) #Cell node 
                                        iLayer+=int(Nlayers[iM])
                            else:
                                K=Kmb[jmb]
                                if Ana_Os and AxialLayers>1:
                                    for iLayer in range(AxialLayers):
                                        rhs_o[iLayer*Ntot+i]+= K*s_membranes[jmb]*(Os_walls[i,iLayer] - Os_cells[j-NwallsJun,0]) #Wall node
                                        rhs_o[iLayer*Ntot+j]+= K*s_membranes[jmb]*(Os_cells[j-NwallsJun,0] - Os_walls[i,iLayer]) #Cell node 
                                else:
                                    temp1=K*s_membranes[jmb]*(Os_walls[i,0] - Os_cells[j-NwallsJun,0]) #Wall node
                                    temp2=K*s_membranes[jmb]*(Os_cells[j-NwallsJun,0] - Os_walls[i,0]) #Cell node 
                                    for iLayer in range(AxialLayers):
                                        rhs_o[iLayer*Ntot+i]+= temp1
                                        rhs_o[iLayer*Ntot+j]+= temp2
                            jmb+=1
            if PileUp==2:
                iLayer=int(Nlayers[0])
                for iM in range(1,Nmaturity):
                    if (Barriers[iM]>0 and Barriers[iM-1]==0):
                        for i in listxyl:
                            rhs_o[iLayer*Ntot+i]+=K_axial[(iLayer-1)*Ntot+i]*s_stele*(Os_cells[i-NwallsJun,iM] - Os_cells[i-NwallsJun,iM-1])
                            rhs_o[(iLayer-1)*Ntot+i]+=K_axial[(iLayer-1)*Ntot+i]*s_stele*(Os_cells[i-NwallsJun,iM-1] - Os_cells[i-NwallsJun,iM])
                        for i in listsieve:
                            if not i in listprotosieve:
                                rhs_o[iLayer*Ntot+i]+=K_axial[(iLayer-1)*Ntot+i]*s_sieve*(Os_cells[i-NwallsJun,iM] - Os_cells[i-NwallsJun,iM-1])
                                rhs_o[(iLayer-1)*Ntot+i]+=K_axial[(iLayer-1)*Ntot+i]*s_sieve*(Os_cells[i-NwallsJun,iM-1] - Os_cells[i-NwallsJun,iM])
                    iLayer+=int(Nlayers[iM])
            for row in range(int(r_discret[0])):
                if nOsWallLayer[row][iMaturity][count]>0:
                    OsWallLayer[row][iMaturity][count]=OsWallLayer[row][iMaturity][count]/nOsWallLayer[row][iMaturity][count]
                if PileUp==2:
                    for iM in range(Nmaturity):
                        if nOsCellLayer[row,iM,count]>0:
                            OsCellLayer[row,iM,count]=OsCellLayer[row,iM,count]/nOsCellLayer[row,iM,count]
                else:
                    if nOsCellLayer[row][iMaturity][count]>0:
                        OsCellLayer[row][iMaturity][count]=OsCellLayer[row][iMaturity][count]/nOsCellLayer[row][iMaturity][count]
            
            #Adding up all BC
            #Elongation BC
            rhs += rhs_e
            
            #Osmotic BC
            rhs += rhs_o
            
            #Soil BC
            for iLayer in range(AxialLayers):
                rhs[iLayer*Ntot:(iLayer+1)*Ntot] += np.multiply(rhs_s[iLayer*Ntot:(iLayer+1)*Ntot],Psi_soil[0][count]*(1-x_rel)+Psi_soil[1][count]*x_rel)
            
            #Xylem BC
            if XylOK>0: #No mature xylem before the Casparian strip stage
                if not isnan(Psi_xyl[1][iMaturity][count]): #Pressure xylem BC at proximal end of segment
                    for cid in listxyl:
                        rhs[-Ntot+cid][0] -= K_axial[-Ntot+cid]*Psi_xyl[1][iMaturity][count]  #Axial conductance of xylem vessels
                        if sparseM==1:
                            matrix_W2[-Ntot+cid,-Ntot+cid] -= K_axial[-Ntot+cid]
                        else:
                            matrix_W[cid][cid] -= K_axial[cid]
                elif not isnan(Flow_xyl[1][0][count]): #Flow xylem BC proximal end of segment
                    i=0
                    for cid in listxyl:
                        rhs[-Ntot+cid][0] += Flow_xyl[1][i+1][count] #(cm^3/d) proximal end of segment
                        i+=1
                if not isnan(Psi_xyl[0][iMaturity][count]): #Pressure xylem BC at distal end of segment
                    for cid in listxyl:
                        rhs[cid][0] -= K_axial[cid]*Psi_xyl[0][iMaturity][count]  #Axial conductance of xylem vessels
                        if sparseM==1:
                            matrix_W2[cid,cid] -= K_axial[cid]
                        else:
                            matrix_W[cid][cid] -= K_axial[cid]
                elif not isnan(Flow_xyl[0][0][count]): #Flow xylem BC at distal end of segment
                    i=0
                    for cid in listxyl:
                        rhs[cid][0] -= Flow_xyl[0][i+1][count] #(cm^3/d) distal end of segment
                        i+=1
            
            #Phloem BC
            if PileUp==2:
                #Proximal phloem BC
                if XylOK==0: #Protophloem only
                    if not isnan(Psi_sieve[1][iMaturity][count]):
                        for cid in listprotosieve:
                            rhs[-Ntot+cid][0] -= K_axial[-Ntot+cid]*Psi_sieve[1][iMaturity][count]  #Axial conductance of phloem sieve tube
                            if sparseM==1:
                                matrix_W2[-Ntot+cid,-Ntot+cid] -= K_axial[-Ntot+cid]
                            else:
                                matrix_W[cid][cid] -= K_axial[cid]
                    elif not isnan(Flow_sieve[1][0][count]):
                        i=1
                        for cid in listsieve: #only Flow_sieve[i] corresponding to listprotosieve should be non-null
                            rhs[-Ntot+cid][0] += Flow_sieve[1][i][count] #(cm^3/d)
                            i+=1
                else:#Includes mature phloem
                    if not isnan(Psi_sieve[1][iMaturity][count]): #Then there is a phloem BC in scenarios (assuming that we did not pick scenarios with and others without)
                        for cid in listsieve: #both proto and metaphloem
                            rhs[-Ntot+cid][0] -= K_axial[-Ntot+cid]*Psi_sieve[1][iMaturity][count]  #Axial conductance of xylem vessels
                            if sparseM==1:
                                matrix_W2[-Ntot+cid,-Ntot+cid] -= K_axial[-Ntot+cid]
                            else:
                                matrix_W[cid][cid] -= K_axial[-Ntot+cid]
                    elif not isnan(Flow_sieve[1][0][count]):
                        i=1
                        for cid in listsieve:
                            rhs[-Ntot+cid][0] += Flow_sieve[1][i][count] #(cm^3/d)
                            i+=1
                #Distal phloem BC
                if Barriers[0]==0: #Only protophloem at distal end
                    if not isnan(Psi_sieve[0][iMaturity][count]):
                        for cid in listprotosieve:
                            rhs[cid][0] -= K_axial[cid]*Psi_sieve[0][iMaturity][count]  #Axial conductance of phloem sieve tube
                            if sparseM==1:
                                matrix_W2[cid,cid] -= K_axial[cid]
                            else:
                                matrix_W[cid][cid] -= K_axial[cid]
                    elif not isnan(Flow_sieve[0][0][count]):
                        i=1
                        for cid in listsieve:  #only Flow_sieve[i] corresponding to listprotosieve should be non-null
                            rhs[cid][0] -= Flow_sieve[0][i][count] #(cm^3/d)
                            i+=1
                else: #Includes mature phloem at distal end
                    if not isnan(Psi_sieve[0][iMaturity][count]): #Then there is a phloem BC in scenarios (assuming that we did not pick scenarios with and others without)
                        for cid in listsieve: #both proto and metaphloem
                            rhs[cid][0] -= K_axial[cid]*Psi_sieve[0][iMaturity][count]  #Axial conductance of xylem vessels
                            if sparseM==1:
                                matrix_W2[cid,cid] -= K_axial[cid]
                            else:
                                matrix_W[cid][cid] -= K_axial[cid]
                    elif not isnan(Flow_sieve[0][0][count]):
                        i=1
                        for cid in listsieve:
                            rhs[cid][0] -= Flow_sieve[0][i][count] #(cm^3/d)
                            i+=1
            else: #PileUp not == 2
                if XylOK==0: #Protophloem only
                    #Proximal phloem BC
                    if not isnan(Psi_sieve[1][iMaturity][count]):
                        for cid in listprotosieve:
                            rhs[-Ntot+cid][0] -= K_axial[-Ntot+cid]*Psi_sieve[1][iMaturity][count]  #Axial conductance of phloem sieve tube
                            if sparseM==1:
                                matrix_W2[-Ntot+cid,-Ntot+cid] -= K_axial[-Ntot+cid]
                            else:
                                matrix_W[cid][cid] -= K_axial[cid]
                    elif not isnan(Flow_sieve[1][0][count]):
                        i=1
                        for cid in listsieve: #only Flow_sieve[i] corresponding to listprotosieve should be non-null
                            rhs[-Ntot+cid][0] += Flow_sieve[1][i][count] #(cm^3/d)
                            i+=1
                    #Distal phloem BC
                    if not isnan(Psi_sieve[0][iMaturity][count]):
                        for cid in listprotosieve:
                            rhs[cid][0] -= K_axial[cid]*Psi_sieve[0][iMaturity][count]  #Axial conductance of phloem sieve tube
                            if sparseM==1:
                                matrix_W2[cid,cid] -= K_axial[cid]
                            else:
                                matrix_W[cid][cid] -= K_axial[cid]
                    elif not isnan(Flow_sieve[0][0][count]):
                        i=1
                        for cid in listsieve: #only Flow_sieve[i] corresponding to listprotosieve should be non-null
                            rhs[cid][0] -= Flow_sieve[0][i][count] #(cm^3/d)
                            i+=1
                elif XylOK>0: #Includes mature phloem
                    #Proximal phloem BC
                    if not isnan(Psi_sieve[1][iMaturity][count]): #Then there is a phloem BC in scenarios (assuming that we did not pick scenarios with and others without)
                        for cid in listsieve: #both proto and metaphloem
                            rhs[-Ntot+cid][0] -= K_axial[-Ntot+cid]*Psi_sieve[1][iMaturity][count]  #Axial conductance of xylem vessels
                            if sparseM==1:
                                matrix_W2[-Ntot+cid,-Ntot+cid] -= K_axial[-Ntot+cid]
                            else:
                                matrix_W[cid][cid] -= K_axial[-Ntot+cid]
                    elif not isnan(Flow_sieve[1][0][count]):
                        i=1
                        for cid in listsieve:
                            rhs[-Ntot+cid][0] += Flow_sieve[1][i][count] #(cm^3/d)
                            i+=1
                    #Distal phloem BC
                    if not isnan(Psi_sieve[0][iMaturity][count]): #Then there is a phloem BC in scenarios (assuming that we did not pick scenarios with and others without)
                        for cid in listsieve: #both proto and metaphloem
                            rhs[cid][0] -= K_axial[cid]*Psi_sieve[0][iMaturity][count]  #Axial conductance of xylem vessels
                            if sparseM==1:
                                matrix_W2[cid,cid] -= K_axial[cid]
                            else:
                                matrix_W[cid][cid] -= K_axial[cid]
                    elif not isnan(Flow_sieve[0][0][count]):
                        i=1
                        for cid in listsieve: #Because PileUp is not equal to 2, the maturity level is uniform
                            rhs[cid][0] -= Flow_sieve[0][i][count] #(cm^3/d)
                            i+=1
            
            t01 = time.perf_counter()
            print(t01-t00, 'seconds process time to adjusting BC for scenario '+str(count))
            
            ##################################################
            ##Solve Doussan equation, results in soln matrix##
            ##################################################
            
            print("Solving water flow")
            t00 = time.perf_counter() #Lasts about 45 sec in the 24mm root
            if sparseM==1:
                soln = spsolve(matrix_W2.tocsr(),rhs)
            else:
                soln = np.linalg.solve(matrix_W,rhs) #Solving the equation to get potentials inside the network
                #Verification that computation was correct
                verif1=np.allclose(np.dot(matrix_W,soln),rhs)
            t01 = time.perf_counter()
            print(t01-t00, 'seconds process time to solve water flow for scenario '+str(count))
            
            if not matrix_analysis==1:
                #Removing Xylem and phloem BC terms in "matrix" in case they would change in the next scenario
                print("Removing axial convective components from matrix W")
                t00 = time.perf_counter()
                if XylOK>0:
                    if not isnan(Psi_xyl[1][iMaturity][count]): #Pressure xylem BC
                        if sparseM==1:
                            for cid in listxyl:
                                matrix_W2[-Ntot+cid,-Ntot+cid] += K_axial[-Ntot+cid]
                        else:
                            for cid in listxyl:
                                matrix_W[cid][cid] += K_axial[-Ntot+cid]
                    if not isnan(Psi_xyl[0][iMaturity][count]): #Pressure xylem BC
                        if sparseM==1:
                            for cid in listxyl:
                                matrix_W2[cid,cid] += K_axial[cid]
                        else:
                            for cid in listxyl:
                                matrix_W[cid][cid] += K_axial[cid]
                if XylOK==0: #Protophloem only
                    if not isnan(Psi_sieve[1][iMaturity][count]):
                        if sparseM==1:
                            for cid in listprotosieve:
                                matrix_W2[-Ntot+cid,-Ntot+cid] += K_axial[-Ntot+cid]
                        else:
                            for cid in listprotosieve:
                                matrix_W[cid][cid] += K_axial[-Ntot+cid]
                    if not isnan(Psi_sieve[0][iMaturity][count]):
                        if sparseM==1:
                            for cid in listprotosieve:
                                matrix_W2[cid,cid] += K_axial[cid]
                        else:
                            for cid in listprotosieve:
                                matrix_W[cid][cid] += K_axial[cid]
                elif XylOK>0: #Includes mature phloem
                    if not isnan(Psi_sieve[1][iMaturity][count]): #Then there is a phloem BC in scenarios (assuming that we did not pick scenarios with and others without)
                        if sparseM==1:
                            for cid in listsieve: #both proto and metaphloem
                                matrix_W2[-Ntot+cid,-Ntot+cid] += K_axial[-Ntot+cid]
                        else:
                            for cid in listsieve: #both proto and metaphloem
                                matrix_W[cid][cid] += K_axial[-Ntot+cid]
                    if not isnan(Psi_sieve[0][iMaturity][count]): #Then there is a phloem BC in scenarios (assuming that we did not pick scenarios with and others without)
                        if sparseM==1:
                            for cid in listsieve: #both proto and metaphloem
                                matrix_W2[cid,cid] += K_axial[cid]
                        else:
                            for cid in listsieve: #both proto and metaphloem
                                matrix_W[cid][cid] += K_axial[cid]
            
            #Flow rates at interfaces
            Q_soil=[]
            Q_soil_neg=[]
            Q_soil_pos=[]
            for ind in Borderwall:
                for iLayer in range(AxialLayers):
                    Q=rhs_s[iLayer*Ntot+ind]*(soln[Ntot*iLayer+ind]-(Psi_soil[0][count]*(1-x_rel[ind])+Psi_soil[1][count]*x_rel[ind]))
                    Q_soil.append(Q) #(cm^3/d) Positive for water flowing into the root, rhs_s is minus the conductance at the soil root interface
                    if Q>0:
                        Q_soil_pos.append(Q)
                        Q_soil_neg.append(0.0)
                    else:
                        Q_soil_neg.append(Q)
                        Q_soil_pos.append(0.0)
                    #if Apo_Contagion==2:
                    #    if Sym_Contagion==2:
                    #            if Q<0.0: #Solute flowing out of the root
                    #                if sparseM==1:
                    #                    for iLayer in range(AxialLayers):
                    #                        matrix_C2[iLayer*Ntot+ind,iLayer*Ntot+ind]+=Q
                    #                else:
                    #                    matrix_C[ind][ind]+=Q
                    #    else:
                    #        if ind not in N_Dirichlet_ini:
                    #            if Q<0.0:
                    #                matrix_ApoC[ind][ind]+=Q
            
            t01 = time.perf_counter()
            print(t01-t00, 'seconds process time to remove axial convective components from matrices W and C')
            
            for ind in Borderjunction:
                for iLayer in range(AxialLayers):
                    Q=rhs_s[iLayer*Ntot+ind]*(soln[Ntot*iLayer+ind]-(Psi_soil[0][count]*(1-x_rel[ind])+Psi_soil[1][count]*x_rel[ind]))
                    Q_soil.append(Q) #(cm^3/d) Positive for water flowing into the root
                    if Q>0:
                        Q_soil_pos.append(Q)
                        Q_soil_neg.append(0.0)
                    else:
                        Q_soil_neg.append(Q)
                        Q_soil_pos.append(0.0)
            
            Q_xyl=zeros((2,Nxyl)) #Flow rate in xylem tubes (cm^3/d) positive towards the shoot
            if XylOK>0:
                if not isnan(Psi_xyl[1][iMaturity][count]): #Xylem pressure BC
                    i=0
                    for cid in listxyl:
                        Q_xyl[1][i]=K_axial[-Ntot+cid]*(soln[-Ntot+cid]-Psi_xyl[1][iMaturity][count])
                        if not isnan(Psi_xyl[0][iMaturity][count]):
                            Q_xyl[0][i]=K_axial[cid]*(Psi_xyl[0][iMaturity][count]-soln[cid])
                        rank=int(Cell_rank[cid-NwallsJun])
                        row=int(rank2row[rank])
                        Q_xyl_layer[row][iMaturity][count] += -Q_xyl[1][i]+Q_xyl[0][i] #Segment xylem water balance (negative if water from other BC got into the xylem)
                        i+=1
                elif not isnan(Flow_xyl[1][0][count]): #Xylem flow BC
                    i=0
                    for cid in listxyl:
                        if test_mass_balance==1:
                            if not isnan(Flow_xyl[0][0][count]):
                                Q_xyl[0][i]=Flow_xyl[0][i+1][count]
                            Q_xyl[1][i]=Q_xyl[0][i]+sum(Q_soil)*Flow_xyl[1][i+1][count]/Flow_xyl[1][0][count]  #(cm^3/d)
                        else:
                            if not isnan(Flow_xyl[0][0][count]):
                                Q_xyl[0][i]=Flow_xyl[0][i+1][count]
                            Q_xyl[1][i]=Flow_xyl[1][i+1][count]  #(cm^3/d)
                        rank=int(Cell_rank[cid-NwallsJun])
                        row=int(rank2row[rank])
                        Q_xyl_layer[row][iMaturity][count] += -Q_xyl[1][i]+Q_xyl[0][i]
                        i+=1
            
            Q_sieve=zeros((2,Nsieve)) #(cm^3/d) Flow positive upwards
            #Proximal phloem BC
            if not isnan(Psi_sieve[1][iMaturity][count]): #Phloem proximal pressure BC
                i=0
                for cid in listsieve: #Q will be 0 for metaphloem if Barrier==0 because of flag
                    if (cid in listprotosieve) or XylOK:
                        Q_sieve[1][i]=K_axial[-Ntot+cid]*(soln[-Ntot+cid]-Psi_sieve[1][iMaturity][count]) #(cm^3/d) Negative for water flowing into phloem tubes at the proximal end
                    i+=1
            elif not isnan(Flow_sieve[1][0][count]): #Phloem proximal flow BC
                i=0
                for cid in listsieve:
                    if cid in listprotosieve or XylOK:
                        Q_sieve[1][i]=Flow_sieve[1][i+1][count]
                    i+=1
            #Distal phloem BC
            if not isnan(Psi_sieve[0][iMaturity][count]):
                i=0
                for cid in listsieve: #Q_sieve will be 0 for metaphloem if Barrier==0 (using the flag)
                    if PileUp==2:
                        if (cid in listprotosieve) or Barriers[0]>0:
                            Q_sieve[0][i]=K_axial[cid]*(Psi_sieve[0][iMaturity][count]-soln[cid]) #(cm^3/d) Positive for water flowing shootwards
                    else:
                        if (cid in listprotosieve) or Barrier>0:
                            Q_sieve[0][i]=K_axial[cid]*(Psi_sieve[0][iMaturity][count]-soln[cid]) #(cm^3/d) Positive for water flowing shootwards
                    i+=1
            elif not isnan(Flow_sieve[0][0][count]):
                i=0
                for cid in listsieve: #Q will be 0 for metaphloem if Barrier==0 (using the flag)
                    if PileUp==2:
                        if (cid in listprotosieve) or Barriers[0]>0:
                            Q_sieve[0][i]=Flow_sieve[0][i+1][count]
                    else:
                        if (cid in listprotosieve) or Barrier>0:
                            Q_sieve[0][i]=Flow_sieve[0][i+1][count]
                    i+=1           
            i=0
            for cid in listsieve:
                rank=int(Cell_rank[cid-NwallsJun])
                row=int(rank2row[rank])
                Q_sieve_layer[row][iMaturity][count] += -Q_sieve[1][i]+Q_sieve[0][i]
                i+=1
            #Elongation flux
            Q_elong=-rhs_e #(cm^3/d) The elongation flux virtually disappears from the cross-section => negative
            for cid in range(Ncells):
                rank=int(Cell_rank[cid])
                row=int(rank2row[rank])
                for iLayer in range(AxialLayers):
                    Q_elong_layer[row][iMaturity][count] += Q_elong[iLayer*Ntot+NwallsJun+cid]
            Q_tot[iMaturity][count]=sum(Q_soil) #(cm^3/d) Total flow rate at root surface
            if abs(float(Q_tot[iMaturity][0]))>1.0E-17: #Changed int() to float()
                for ind in range(NwallsJun,Ntot): #NwallsJun is the index of the first cell
                    cellnumber1=ind-NwallsJun
                    rank = int(Cell_rank[cellnumber1])
                    row = int(rank2row[rank])
                    if rank == 1: #Exodermis
                        for iLayer in range(AxialLayers):
                            iL=TopLayer-AxialLayers+iLayer
                            PsiCellLayer[row,iL,count] += soln[iLayer*Ntot+ind]/AxialLayers*(STFcell_plus[cellnumber1,iL]+abs(STFcell_minus[cellnumber1,iL]))/(STFlayer_plus[row,iL]+abs(STFlayer_minus[row,iL])+STFlayer_plus[row+1,iL]+abs(STFlayer_minus[row+1,iL])) #(hPa)
                            PsiCellLayer[row+1,iL,count] += soln[iLayer*Ntot+ind]/AxialLayers*(STFcell_plus[cellnumber1,iL]+abs(STFcell_minus[cellnumber1,iL]))/(STFlayer_plus[row,iL]+abs(STFlayer_minus[row,iL])+STFlayer_plus[row+1,iL]+abs(STFlayer_minus[row+1,iL])) #(hPa)
                    elif rank == 3: #Endodermis
                        if any(passage_cell_ID==array(cellnumber1)):# and Barrier==2: #Passage cell
                            for iLayer in range(AxialLayers):
                                iL=TopLayer-AxialLayers+iLayer
                                PsiCellLayer[row+1,iL,count] += soln[iLayer*Ntot+ind]/AxialLayers*(STFcell_plus[cellnumber1,iL]+abs(STFcell_minus[cellnumber1,iL]))/(STFlayer_plus[row+1,iL]+abs(STFlayer_minus[row+1,iL])+STFlayer_plus[row+2,iL]+abs(STFlayer_minus[row+2,iL])) #(hPa)
                                PsiCellLayer[row+2,iL,count] += soln[iLayer*Ntot+ind]/AxialLayers*(STFcell_plus[cellnumber1,iL]+abs(STFcell_minus[cellnumber1,iL]))/(STFlayer_plus[row+1,iL]+abs(STFlayer_minus[row+1,iL])+STFlayer_plus[row+2,iL]+abs(STFlayer_minus[row+2,iL])) #(hPa)
                        else:
                            for iLayer in range(AxialLayers):
                                iL=TopLayer-AxialLayers+iLayer
                                PsiCellLayer[row,iL,count] += soln[iLayer*Ntot+ind]/AxialLayers*(STFcell_plus[cellnumber1,iL]+abs(STFcell_minus[cellnumber1,iL]))/(STFlayer_plus[row,iL]+abs(STFlayer_minus[row,iL])+STFlayer_plus[row+3,iL]+abs(STFlayer_minus[row+3,iL])) #(hPa)
                                PsiCellLayer[row+3,iL,count] += soln[iLayer*Ntot+ind]/AxialLayers*(STFcell_plus[cellnumber1,iL]+abs(STFcell_minus[cellnumber1,iL]))/(STFlayer_plus[row,iL]+abs(STFlayer_minus[row,iL])+STFlayer_plus[row+3,iL]+abs(STFlayer_minus[row+3,iL])) #(hPa)
                                #if not Barrier==2:
                                #    PsiCellLayer[row+1][iL][count] = nan
                                #    PsiCellLayer[row+2][iL][count] = nan
                    elif (ind not in listxyl) or XylOK==0: #Not for mature xylem
                        for iLayer in range(AxialLayers):
                            iL=TopLayer-AxialLayers+iLayer
                            PsiCellLayer[row,iL,count] += soln[iLayer*Ntot+ind]/AxialLayers*(STFcell_plus[cellnumber1,iL]+abs(STFcell_minus[cellnumber1,iL]))/(STFlayer_plus[row,iL]+abs(STFlayer_minus[row,iL])) #(hPa)
            
            if count==iEquil_xyl:
                if XylOK>0 and isnan(Psi_xyl[1][iMaturity][count]):
                    Psi_xyl[1][iMaturity][count]=0.0
                    for cid in listxyl:
                        Psi_xyl[1][iMaturity][count]+=soln[-Ntot+cid]/Nxyl
            if count==iEquil_sieve:
                if XylOK>0:
                    if isnan(Psi_sieve[1][iMaturity][count]):
                        Psi_sieve[1][iMaturity][count]=0.0
                        for cid in listsieve:
                            Psi_sieve[1][iMaturity][count]+=soln[-Ntot+cid]/Nsieve #Average of phloem water pressures
                elif XylOK==0:
                    if isnan(Psi_sieve[1][iMaturity][count]):
                        Psi_sieve[1][iMaturity][count]=0.0
                        for cid in listprotosieve:
                            Psi_sieve[1][iMaturity][count]+=soln[-Ntot+cid]/Nprotosieve #Average of protophloem water pressures
            
            print("Uptake rate per unit root length: soil ",(sum(Q_soil)/(Totheight)/1.0E-04),"cm^2/d, xylem ",((sum(Q_xyl[0])-sum(Q_xyl[1]))/(Totheight)/1.0E-04),"cm^2/d, phloem ",((sum(Q_sieve[0])-sum(Q_sieve[1]))/(Totheight)/1.0E-04),"cm^2/d, elongation ",(sum(Q_elong)/(Totheight)/1.0E-04),"cm^2/d")
            print("Uptake rate: soil ",(sum(Q_soil)),"cm^3/d, xylem ",((sum(Q_xyl[0])-sum(Q_xyl[1]))),"cm^3/d, phloem ",((sum(Q_sieve[0])-sum(Q_sieve[1]))),"cm^3/d, elongation ",(sum(Q_elong)),"cm^3/d")
            if not isnan(sum(Q_sieve[0])-sum(Q_sieve[1])):
                print("Mass balance error:",(sum(Q_soil)+sum(Q_xyl[0])-sum(Q_xyl[1])+sum(Q_sieve[0])-sum(Q_sieve[1])+sum(Q_elong)),"cm^3/d")
            else:
                print("Mass balance error:",(sum(Q_soil)+sum(Q_xyl[0])-sum(Q_xyl[1])+sum(Q_elong)),"cm^3/d")
            
            #################################################################
            ##Calul of Fluxes between nodes and Creating the edge_flux_list##
            #################################################################
            
            print("Creating a list of water fluxes through walls, membranes and PD, & adjust matrix_C")
            t00 = time.perf_counter() #Lasts about 973 sec in the 24mm root
            
            #Filling the fluxes list
            MembraneFlowDensity=[]
            WallFlowDensity=[]
            WallFlowDensity_cos=[]
            PlasmodesmFlowDensity=[]
            Fjw_list=[]
            Fcw_list=[]
            Fcc_list=[]
            Faxial_sym_apo_list=[] #Axial symplastic and apoplastic flow
            Faxial_TM_list=[] #Axial transmembrane flow
            psi_o_cell=empty((AxialLayers),dtype=dp)
            jmb=0 #Index for membrane conductance vector
            for node, edges in G.adjacency_iter() : #adjacency_iter returns an iterator of (node, adjacency dict) tuples for all nodes. This is the fastest way to look at every edge. For directed graphs, only outgoing adjacencies are included.
                i = indice[node] #Node ID number
                psi=empty((AxialLayers),dtype=dp)
                for iLayer in range(AxialLayers): #The command psi=sol[...] changes soln when psi is updated
                    psi[iLayer] = float(soln[i+iLayer*Ntot]) #Node water potential
                #ind_o_cell = nan #Opposite cell index
                #Here we count surrounding cell types in order to know if the wall is part of an apoplastic barrier, as well as to know on which side of the exodermis or endodermis the membrane is located
                count_endo=0 #total number of endodermis cells around the wall
                count_peri=0 #total number of pericycle cells around the wall
                count_PPP=0 #total number of PPP cells arount the wall
                count_exo=0 #total number of exodermis cells around the wall
                count_epi=0 #total number of epidermis cells around the wall
                count_stele=0 #total number of stelar parenchyma cells around the wall
                count_stele_overall=0 #total number of stele cells (of any type) around the wall
                count_comp=0 #total number of companion cells around the wall
                count_sieve=0 #total number of stelar parenchyma cells around the wall
                count_xyl=0 #total number of xylem cells around the wall
                count_cortex=0 #total number of phloem sieve cells around the wall
                count_passage=0 #total number of passage cells around the wall
                count_interC=0 #total number of intercellular spaces around the wall
                if PileUp==2:
                    PDs=ones((AxialLayers,1))
                noPD=False #Initializes the flag for wall connected to an intercellular space -> does not have plasmodesmata
                if i<Nwalls: #wall ID
                    for neighbour, eattr in edges.items(): #Loop on connections (edges)
                        if eattr['path'] == 'membrane': #Wall connection
                            if any(passage_cell_ID==array((indice[neighbour])-NwallsJun)):
                                count_passage+=1
                            if any(InterCid==array((indice[neighbour])-NwallsJun)):
                                count_interC+=1
                            if G.node[neighbour]['cgroup']==3:#Endodermis
                                count_endo+=1
                            elif G.node[neighbour]['cgroup']==13 or G.node[neighbour]['cgroup']==19 or G.node[neighbour]['cgroup']==20:#Xylem cell or vessel
                                count_xyl+=1
                            elif G.node[neighbour]['cgroup']==16 or G.node[neighbour]['cgroup']==21:#Pericycle or stele
                                count_peri+=1
                                if neighbour in PPP:
                                    count_PPP+=1
                            elif G.node[neighbour]['cgroup']==1:#Exodermis
                                count_exo+=1
                            elif G.node[neighbour]['cgroup']==2:#Epidermis
                                count_epi+=1
                            elif G.node[neighbour]['cgroup']==4:#Cortex
                                count_cortex+=1
                            elif G.node[neighbour]['cgroup']==5:#Stelar parenchyma
                                count_stele+=1
                            elif G.node[neighbour]['cgroup']==11 or G.node[neighbour]['cgroup']==23:#Phloem sieve tube
                                count_sieve+=1
                            elif G.node[neighbour]['cgroup']==12 or G.node[neighbour]['cgroup']==26 or G.node[neighbour]['cgroup']==24:#Companion cell
                                count_comp+=1
                            if G.node[neighbour]['cgroup']>4:#Stele overall
                                count_stele_overall+=1
                psi_o_cell[:] = nan #Opposite cell water potential
                for neighbour, eattr in edges.items(): #Loop on connections (edges)
                    j = indice[neighbour] #Neighbouring node ID number
                    #if j > i: #Only treating the information one way to save time
                    psin = empty((AxialLayers),dtype=dp)
                    for iLayer in range(AxialLayers):
                        psin[iLayer] = float(soln[j+iLayer*Ntot]) #Neighbouring node water potential
                    path = eattr['path'] #eattr is the edge attribute (i.e. connection type)
                    if i<Nwalls:
                        if Paraview==1 or ParTrack==1 or Apo_Contagion>0 or Sym_Contagion>0:
                            if path == "wall":
                                if PileUp==2:
                                    K=empty((AxialLayers))
                                    iLayer=0
                                    for iM in range(Nmaturity):
                                        temp=1.0E-04*((eattr['lat_dist']*Xwalls+heights[iLayer]-thickness)*thickness)/eattr['length'] #Wall section to length ratio (cm)
                                        K[iLayer:iLayer+int(Nlayers[iM])] = kwi[i,iM]*temp #Junction-Wall conductance (cm^3/hPa/d)
                                        iLayer+=int(Nlayers[iM])
                                    ##K was vector in array, so K[0] is the vector we want, just a question of format
                                    #Fjw = K[0] * (psin - psi) * sign(j-i) #(cm^3/d) Water flow rate positive from junction to wall
                                    Fjw = K * (psin - psi) * sign(j-i) #(cm^3/d) Water flow rate positive from junction to wall
                                else:
                                    temp=1.0E-04*((eattr['lat_dist']*Xwalls+height-thickness)*thickness)/eattr['length']
                                    #K = eattr['kw']*1.0E-04*((eattr['lat_dist']+height)*eattr['thickness']-square(eattr['thickness']))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)
                                    if Xylem_pieces and ((count_interC>=2 and Barrier>0) or count_xyl==2): #"Fake wall" splitting an intercellular space or a xylem cell in two
                                        K = 1.0E-16 #Non conductive
                                    elif count_cortex>=2: #wall between two cortical cells
                                        K = kw_cortex_cortex*temp #Junction-Wall conductance (cm^3/hPa/d)
                                    elif count_endo>=2: #wall between two endodermis cells
                                        K = kw_endo_endo*temp #Junction-Wall conductance (cm^3/hPa/d)
                                    #elif count_stele_overall>0 and count_endo>0: #wall between endodermis and pericycle
                                    #    if count_passage>0:
                                    #        K = kw_passage*temp
                                    #    else:
                                    #        K = kw_endo_peri*temp #Junction-Wall conductance (cm^3/hPa/d)
                                    #elif count_stele_overall==0 and count_endo==1: #wall between endodermis and cortex
                                    #    if count_passage>0:
                                    #        K = kw_passage*temp
                                    #    else:
                                    #        K = kw_endo_cortex*temp #Junction-Wall conductance (cm^3/hPa/d)
                                    elif count_exo>=2: #wall between two exodermis cells
                                        K = kw_exo_exo*temp #Junction-Wall conductance (cm^3/hPa/d)
                                    else: #other walls
                                        K = kw*temp #Junction-Wall conductance (cm^3/hPa/d)
                                    Fjw = K * (psin - psi) * sign(j-i) #(cm^3/d) Water flow rate positive from junction to wall
                                Fjw_list.append((i,j,Fjw))
                                if Paraview==1 or ParTrack==1:
                                    #The ordering in WallFlowDensity will correspond to the one of ThickWallsX, saved for display only
                                    cos_angle=(position[i][0]-position[j][0])/(hypot(position[j][0]-position[i][0],position[j][1]-position[i][1])) #Vectors junction1-wall
                                    if PileUp==2:
                                        WallFlowDensity.append((i,j, Fjw / (((eattr['lat_dist']*Xwalls+heights[:,0]-thickness)*thickness)*1.0E-08))) # (cm/d) Positive towards lower node ID 
                                        WallFlowDensity_cos.append((i,j, cos_angle * Fjw / (((eattr['lat_dist']*Xwalls+heights[:,0]-thickness)*thickness)*1.0E-08))) # (cm/d) Positive towards lower node ID 
                                    else:
                                        WallFlowDensity.append((i,j, Fjw / (((eattr['lat_dist']*Xwalls+height-thickness)*thickness)*1.0E-08))) # (cm/d) Positive towards lower node ID 
                                        WallFlowDensity_cos.append((i,j, cos_angle * Fjw / (((eattr['lat_dist']*Xwalls+height-thickness)*thickness)*1.0E-08))) # (cm/d) Positive towards lower node ID 
                                #if Ana_Os and Os_soil[5][count]*Os_xyl[5][count]==1:
                                if Apo_Contagion==2:
                                    if Sym_Contagion==2: # Apo & Sym contagion
                                        for iLayer in range(AxialLayers):
                                            if Fjw[iLayer]>0: #Flow from junction to wall
                                                    if sparseM==1:
                                                        matrix_C2[iLayer*Ntot+i,iLayer*Ntot+j] += Fjw[iLayer]
                                                    else:
                                                        matrix_C[i][j] += Fjw[iLayer]
                                                    if sparseM==1:
                                                        matrix_C2[iLayer*Ntot+j,iLayer*Ntot+j] -= Fjw[iLayer]
                                                    else:
                                                        matrix_C[j][j] -= Fjw[iLayer]
                                            else: #Flow from wall to junction
                                                    if sparseM==1:
                                                        matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] += Fjw[iLayer]
                                                    else:
                                                        matrix_C[i][i] += Fjw[iLayer]
                                                    if sparseM==1:
                                                        matrix_C2[iLayer*Ntot+j,iLayer*Ntot+i] -= Fjw[iLayer]
                                                    else:
                                                        matrix_C[j][i] -= Fjw[iLayer]
                                    else: #Only Apo contagion
                                        if Fjw>0: #Flow from junction to wall
                                            if i not in N_Dirichlet_ini:
                                                matrix_ApoC[i][j] += Fjw
                                            if j not in N_Dirichlet_ini:
                                                matrix_ApoC[j][j] -= Fjw
                                        else: #Flow from wall to junction
                                            if i not in N_Dirichlet_ini:
                                                matrix_ApoC[i][i] += Fjw
                                            if j not in N_Dirichlet_ini:
                                                matrix_ApoC[j][i] -= Fjw
                                
                                if Apo_Contagion==1:
                                    if Fjw>0:
                                        Apo_connec_flow[j][nApo_connec_flow[j]]=i
                                        nApo_connec_flow[j]+=1
                                    elif Fjw<0:
                                        Apo_connec_flow[i][nApo_connec_flow[i]]=j
                                        nApo_connec_flow[i]+=1
                            elif path == "membrane": #Membrane connection
                                if PileUp==2:
                                    Fcw=empty((AxialLayers))
                                    iLayer=0
                                    if Ana_Os:
                                        for iM in range(Nmaturity):
                                            Fcw[iLayer:iLayer+int(Nlayers[iM])] = Kmb[jmb,iM] * (psi[iLayer:iLayer+int(Nlayers[iM])] - psin[iLayer:iLayer+int(Nlayers[iM])] + s_membranes[jmb,iM]*(Os_walls[i,iLayer:iLayer+int(Nlayers[iM])] - Os_cells[j-NwallsJun,iM])) #(cm^3/d) Water flow rate positive from wall to protoplast
                                            if G.node[j]['cgroup']==13 or G.node[j]['cgroup']==19 or G.node[j]['cgroup']==20: #xylem cell or vessel
                                                if Barriers[iM]>0: #Xylem vessel
                                                    PDs[iLayer:int(iLayer+Nlayers[iM])]=0
                                                elif Barriers[iM]==0: #Xylem cell
                                                    if (count_xyl==2 and Xylem_pieces):
                                                        PDs[iLayer:int(iLayer+Nlayers[iM])]=0
                                            elif (j-NwallsJun in InterCid) and Barriers[iM]>0: #the neighbour is an intercellular space "cell"
                                                PDs[iLayer:int(iLayer+Nlayers[iM])]=0
                                            iLayer+=int(Nlayers[iM])
                                        Fcw_list.append((i,j,-Fcw,s_membranes[jmb,:])) #Water flow rate positive from protoplast to wall
                                        if Paraview==1 or ParTrack==1:
                                            #Flow densities calculation
                                            #The ordering in MembraneFlowDensity will correspond to the one of ThickWalls, saved for display only 
                                            MembraneFlowDensity.append(Fcw / (1.0E-08*(heights[:,0]+eattr['dist'])*eattr['length']))
                                    else:
                                        for iM in range(Nmaturity):
                                            Fcw[iLayer:iLayer+int(Nlayers[iM])] = Kmb[jmb,iM] * (psi[iLayer:iLayer+int(Nlayers[iM])] - psin[iLayer:iLayer+int(Nlayers[iM])] + s_membranes[jmb,iM]*(Os_walls[i,iM] - Os_cells[j-NwallsJun,iM])) #(cm^3/d) Water flow rate positive from wall to protoplast
                                            if G.node[j]['cgroup']==13 or G.node[j]['cgroup']==19 or G.node[j]['cgroup']==20: #xylem cell or vessel
                                                if Barriers[iM]>0: #Xylem vessel
                                                    PDs[iLayer:int(iLayer+Nlayers[iM])]=0
                                                elif Barriers[iM]==0: #Xylem cell
                                                    if (count_xyl==2 and Xylem_pieces):
                                                        PDs[iLayer:int(iLayer+Nlayers[iM])]=0
                                            elif (j-NwallsJun in InterCid) and Barriers[iM]>0: #the neighbour is an intercellular space "cell"
                                                PDs[iLayer:int(iLayer+Nlayers[iM])]=0
                                            iLayer+=int(Nlayers[iM])
                                        Fcw_list.append((i,j,-Fcw,s_membranes[jmb,:])) #Water flow rate positive from protoplast to wall
                                        if Paraview==1 or ParTrack==1:
                                            #Flow densities calculation
                                            #The ordering in MembraneFlowDensity will correspond to the one of ThickWalls, saved for display only 
                                            MembraneFlowDensity.append(Fcw / (1.0E-08*(heights[:,0]+eattr['dist'])*eattr['length']))
                                else:
                                    Fcw=empty((AxialLayers))
                                    if Ana_Os and AxialLayers>1:
                                        for iLayer in range(AxialLayers):
                                            Fcw[iLayer] = Kmb[jmb] * (psi[iLayer] - psin[iLayer] + s_membranes[jmb]*(Os_walls[i,iLayer] - Os_cells[j-NwallsJun,0])) #(cm^3/d) Water flow rate positive from wall to protoplast
                                    else:
                                        Fcw = Kmb[jmb] * (psi - psin + s_membranes[jmb]*(Os_walls[i,0] - Os_cells[j-NwallsJun,0])) #(cm^3/d) Water flow rate positive from wall to protoplast
                                    if Barrier>0: #Xylem vessel
                                        noPD=True
                                    elif Barrier==0: #Xylem cell
                                        if (count_xyl==2 and Xylem_pieces):
                                            noPD=True
                                    Fcw_list.append((i,j,-Fcw,s_membranes[jmb])) #Water flow rate positive from protoplast to wall
                                    if Paraview==1 or ParTrack==1:
                                        #Flow densities calculation
                                        #The ordering in MembraneFlowDensity will correspond to the one of ThickWalls, saved for display only 
                                        MembraneFlowDensity.append(Fcw / (1.0E-08*(height+eattr['dist'])*eattr['length']))
                                ####Solute convection across membranes####
                                if Apo_Contagion==2 and Sym_Contagion==2:
                                    for iLayer in range(AxialLayers):
                                        if Fcw[iLayer]>0: #Flow from wall to protoplast
                                                if D2O1==1:#Solute that moves across membranes like water 
                                                    if sparseM==1:
                                                        matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] -= Fcw[iLayer]
                                                    else:
                                                        matrix_C[i][i] -= Fcw[iLayer]
                                                else: #Solute that moves across membranes independently of water (the membrane is possibly not one) 
                                                    if sparseM==1:
                                                        if PileUp==2:
                                                            matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] -= Fcw[iLayer]*(1-s_membranes[jmb,MaturL[iLayer]])
                                                        else:
                                                            matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] -= Fcw[iLayer]*(1-s_membranes[jmb])
                                                    else:
                                                        matrix_C[i][i] -= Fcw[iLayer]*(1-s_membranes[jmb])
                                                if D2O1==1:#Solute that moves across membranes like water 
                                                    if sparseM==1:
                                                        matrix_C2[iLayer*Ntot+j,iLayer*Ntot+i] += Fcw[iLayer]
                                                    else:
                                                        matrix_C[j][i] += Fcw[iLayer]
                                                else: #Solute that moves across membranes independently of water (the membrane is possibly not one) 
                                                    if sparseM==1:
                                                        if PileUp==2:
                                                            matrix_C2[iLayer*Ntot+j,iLayer*Ntot+i] += Fcw[iLayer]*(1-s_membranes[jmb,MaturL[iLayer]])
                                                        else:
                                                            matrix_C2[iLayer*Ntot+j,iLayer*Ntot+i] += Fcw[iLayer]*(1-s_membranes[jmb])
                                                    else:
                                                        matrix_C[j][i] += Fcw[iLayer]*(1-s_membranes[jmb])
                                        else: #Flow from protoplast to wall
                                                if D2O1==1:#Solute that moves across membranes like water 
                                                    if sparseM==1:
                                                        matrix_C2[iLayer*Ntot+j,iLayer*Ntot+j] += Fcw[iLayer]
                                                    else:
                                                        matrix_C[j][j] += Fcw[iLayer]
                                                else: #Solute that moves across membranes independently of water (the membrane is possibly not one) 
                                                    if sparseM==1:
                                                        if PileUp==2:
                                                            matrix_C2[iLayer*Ntot+j,iLayer*Ntot+j] += Fcw[iLayer]*(1-s_membranes[jmb,MaturL[iLayer]])
                                                        else:
                                                            matrix_C2[iLayer*Ntot+j,iLayer*Ntot+j] += Fcw[iLayer]*(1-s_membranes[jmb])
                                                    else:
                                                        matrix_C[j][j] += Fcw[iLayer]*(1-s_membranes[jmb])
                                                if D2O1==1:#Solute that moves across membranes like water 
                                                    if sparseM==1:
                                                        matrix_C2[iLayer*Ntot+i,iLayer*Ntot+j] -= Fcw[iLayer]
                                                    else:
                                                        matrix_C[i][j] -= Fcw[iLayer]
                                                else: #Solute that moves across membranes independently of water (the membrane is possibly not one) 
                                                    if sparseM==1:
                                                        if PileUp==2:
                                                            matrix_C2[iLayer*Ntot+i,iLayer*Ntot+j] -= Fcw[iLayer]*(1-s_membranes[jmb,MaturL[iLayer]])
                                                        else:
                                                            matrix_C2[iLayer*Ntot+i,iLayer*Ntot+j] -= Fcw[iLayer]*(1-s_membranes[jmb])
                                                    else:
                                                        matrix_C[i][j] -= Fcw[iLayer]*(1-s_membranes[jmb])
                                    
                                #Macroscopic distributed parameter for transmembrane flow
                                #Discretization based on cell layers and apoplasmic barriers
                                rank = int(Cell_rank[j-NwallsJun])
                                row = int(rank2row[rank])
                                if rank == 1 and count_epi > 0: #Outer exodermis
                                    row += 1
                                if rank == 3 and count_cortex > 0: #Outer endodermis
                                    if any(passage_cell_ID==array(j-NwallsJun)):# and Barrier==2:
                                        row += 2
                                    else:
                                        row += 3
                                elif rank == 3 and count_stele_overall > 0: #Inner endodermis
                                    if any(passage_cell_ID==array(j-NwallsJun)):# and Barrier==2:
                                        row += 1
                                #Flow = K * (psi - psin + s_membranes[jmb]*(Os_walls[i] - Os_cells[j-NwallsJun]))
                                if PileUp==2:
                                    iLayer=0
                                    for iM in range(Nmaturity):
                                        if ((j-NwallsJun not in InterCid) and (j not in listxyl)) or Barriers[iM]==0: #No aerenchyma in the elongation zone
                                            for iL in range(iLayer,int(iLayer+Nlayers[iM])):
                                                if Fcw[iL] > 0 :
                                                    UptakeLayer_plus[row][iL][count] += Fcw[iL] #grouping membrane flow rates in cell layers
                                                else:
                                                    UptakeLayer_minus[row][iL][count] += Fcw[iL]
                                        if Kmb[jmb,iM]>1.0e-18: #Not an impermeable wall
                                            PsiWallLayer[row,iLayer:int(iLayer+Nlayers[iM]),count] += psi[iLayer:int(iLayer+Nlayers[iM])]
                                            NWallLayer[row,iMaturity,count] += 1
                                        iLayer+=int(Nlayers[iM])
                                else:
                                    if ((j-NwallsJun not in InterCid) and (j not in listxyl)) or Barrier==0: #No aerenchyma in the elongation zone
                                        for iLayer in range(AxialLayers):
                                            iL=TopLayer-AxialLayers+iLayer
                                            if Fcw[iLayer] > 0 : #Previously "Flow"
                                                UptakeLayer_plus[row][iL][count] += Fcw[iLayer] #grouping membrane flow rates in cell layers #Previously "Flow"
                                            else:
                                                UptakeLayer_minus[row][iL][count] += Fcw[iLayer] #Previously "Flow"
                                    
                                    if K>1.0e-18: #Not an impermeable wall
                                        for iLayer in range(AxialLayers):
                                            iL=TopLayer-AxialLayers+iLayer
                                            PsiWallLayer[row,iL,count] += psi[iLayer]
                                        NWallLayer[row,iMaturity,count] += 1
                                jmb+=1
                                #Plasmodesmata part
                                if isnan(psi_o_cell[0]):
                                    for iLayer in range(AxialLayers):
                                        psi_o_cell[iLayer]=float(psin[iLayer])
                                    ind_o_cell=j
                                else:
                                    if PD_freq_source==2:
                                        #We have to find back the wall crossed by the PD from the 2 cells ID
                                        wids=[]
                                        for r in Cell2Wall_loop[j-NwallsJun]:
                                            wids.append(r.get("id"))
                                        for r in Cell2Wall_loop[ind_o_cell-NwallsJun]:
                                            wid1=r.get("id")
                                            if wid1 in wids:
                                                break
                                        if PileUp==2:
                                            temp=Kpl*(psin - psi_o_cell)*float(Wall_PD[int(wid1)])*heights[:,0]/1.0E-4 # Units: cm^2/d, based on total number of plasmodesmata on this wall
                                        else:
                                            temp=Kpl*(psin - psi_o_cell)*float(Wall_PD[int(wid1)])*height/1.0E-4 # Units: cm^2/d, based on total number of plasmodesmata on this wall
                                    else:
                                        if noPD: #No plasmodesmata because the wall i is connected to an intercellular space or xylem vessel
                                            temp=zeros((AxialLayers)) #The ordering in PlasmodesmFlowDensity will correspond to the one of ThickWalls except for boderline walls, saved for display only                        
                                        elif count_epi==1 and count_exo==1: #wall between epidermis and exodermis
                                            temp=Kpl*Fplxheight_epi_exo * (psin - psi_o_cell)
                                        elif (count_exo==1 or count_epi==1) and count_cortex==1: #wall between exodermis and cortex
                                            temp1=float(Kplrange[iPD].get("cortex_factor"))
                                            temp=Kpl*2*temp1/(temp1+1)*Fplxheight_outer_cortex * Length_outer_cortex_tot / Length_outer_cortex_nospace * (psin - psi_o_cell)
                                        elif count_cortex==2: #wall between cortical cells
                                            temp1=float(Kplrange[iPD].get("cortex_factor"))
                                            temp=Kpl*temp1*Fplxheight_cortex_cortex * Length_cortex_cortex_tot / Length_cortex_cortex_nospace * (psin - psi_o_cell)
                                        elif count_cortex==1 and count_endo==1: #wall between cortex and endodermis
                                            temp1=float(Kplrange[iPD].get("cortex_factor"))
                                            temp=Kpl*2*temp1/(temp1+1)*Fplxheight_cortex_endo * Length_cortex_endo_tot / Length_cortex_endo_nospace * (psin - psi_o_cell)
                                        elif count_endo==2: #wall between endodermal cells
                                            temp=Kpl*Fplxheight_endo_endo * (psin - psi_o_cell)
                                        elif count_stele_overall>0 and count_endo>0: #wall between endodermis and pericycle
                                            if count_PPP>0:
                                                temp1=float(Kplrange[iPD].get("PPP_factor"))
                                            else:
                                                temp1=1
                                            temp=Kpl*2*temp1/(temp1+1)*Fplxheight_endo_peri * (psin - psi_o_cell)
                                        elif count_stele==2: #wall between stelar parenchyma cells
                                            temp=Kpl*Fplxheight_stele_stele * (psin - psi_o_cell)
                                        elif count_peri>0 and count_stele==1: #wall between stele and pericycle
                                            if count_PPP>0:
                                                temp1=float(Kplrange[iPD].get("PPP_factor"))
                                            else:
                                                temp1=1
                                            temp=Kpl*2*temp1/(temp1+1)*Fplxheight_peri_stele * (psin - psi_o_cell)
                                        elif count_comp==1 and count_stele==1: #wall between stele and companion cell
                                            temp1=float(Kplrange[iPD].get("PCC_factor"))
                                            temp=Kpl*2*temp1/(temp1+1)*Fplxheight_stele_comp * (psin - psi_o_cell)
                                        elif count_peri==1 and count_comp==1: #wall between pericycle and companion cell
                                            temp1=float(Kplrange[iPD].get("PCC_factor"))
                                            if count_PPP>0:
                                                temp2=float(Kplrange[iPD].get("PPP_factor"))
                                            else:
                                                temp2=1
                                            temp=Kpl*2*temp1*temp2/(temp1+temp2)*Fplxheight_peri_comp * (psin - psi_o_cell)
                                        elif count_comp==2: #wall between companion cells 
                                            temp1=float(Kplrange[iPD].get("PCC_factor"))
                                            temp=Kpl*temp1*Fplxheight_comp_comp * (psin - psi_o_cell)
                                        elif count_comp==1 and count_sieve==1: #wall between companion cell and sieve tube
                                            temp1=float(Kplrange[iPD].get("PCC_factor"))
                                            temp=Kpl*2*temp1/(temp1+1)*Fplxheight_comp_sieve * (psin - psi_o_cell)
                                        elif count_peri==1 and count_sieve==1: #wall between stele and sieve tube
                                            temp=Kpl*Fplxheight_peri_sieve * (psin - psi_o_cell)
                                        elif count_stele==1 and count_sieve==1: #wall between stele and pericycle
                                            if count_PPP>0:
                                                temp1=float(Kplrange[iPD].get("PPP_factor"))
                                            else:
                                                temp1=1
                                            temp=Kpl*2*temp1/(temp1+1)*Fplxheight_stele_sieve * (psin - psi_o_cell)
                                        else: #Default plasmodesmatal frequency
                                            temp=Kpl*Fplxheight * (psin - psi_o_cell)  #The ordering in PlasmodesmFlowDensity will correspond to the one of ThickWalls except for boderline walls, saved for display only 
                                    if PileUp==2:
                                        PlasmodesmFlowDensity.append(PDs[:,0]*temp/(1.0E-04*heights[:,0]))
                                        PlasmodesmFlowDensity.append(-PDs[:,0]*temp/(1.0E-04*heights[:,0]))
                                        Fcc=PDs[:,0]*temp*1.0E-04*eattr['length']*sign(j-ind_o_cell)
                                    else:
                                        PlasmodesmFlowDensity.append(temp/(1.0E-04*height))
                                        PlasmodesmFlowDensity.append(-temp/(1.0E-04*height))
                                        Fcc=temp*1.0E-04*eattr['length']*sign(j-ind_o_cell)
                                    #if (ind_o_cell==454 or ind_o_cell==448) and (j==454 or j==448):
                                    #    print('Fcc PD xyl =',Fcc)
                                    if ind_o_cell<j:
                                        Fcc_list.append((ind_o_cell,j,Fcc)) #(cm^3/d) Water flow rate positive from high index to low index cell
                                    else:
                                        Fcc_list.append((j,ind_o_cell,Fcc))
                                    #if Ana_Os:
                                    if Sym_Contagion==2: #Convection across plasmodesmata
                                        if Apo_Contagion==2: #Apo & Sym Contagion
                                            for iLayer in range(AxialLayers):
                                                if Fcc[iLayer]>0: #Flow from high index to low index cell
                                                    if ind_o_cell<j: #From j to ind_o_cell
                                                            if sparseM==1:
                                                                matrix_C2[iLayer*Ntot+j,iLayer*Ntot+j] -= Fcc[iLayer]
                                                            else:
                                                                matrix_C[j][j] -= Fcc[iLayer]
                                                            if sparseM==1:
                                                                matrix_C2[iLayer*Ntot+ind_o_cell,iLayer*Ntot+j] += Fcc[iLayer]
                                                            else:
                                                                matrix_C[ind_o_cell][j] += Fcc[iLayer]
                                                    else: #From ind_o_cell to j
                                                            if sparseM==1:
                                                                matrix_C2[iLayer*Ntot+ind_o_cell,iLayer*Ntot+ind_o_cell] -= Fcc[iLayer]
                                                            else:
                                                                matrix_C[ind_o_cell][ind_o_cell] -= Fcc[iLayer]
                                                            if sparseM==1:
                                                                matrix_C2[iLayer*Ntot+j,iLayer*Ntot+ind_o_cell] += Fcc[iLayer]
                                                            else:
                                                                matrix_C[j][ind_o_cell] += Fcc[iLayer]
                                                else: #Flow from low index to high index cell
                                                    if ind_o_cell<j: #From ind_o_cell to j
                                                            if sparseM==1:
                                                                matrix_C2[iLayer*Ntot+ind_o_cell,iLayer*Ntot+ind_o_cell] += Fcc[iLayer]
                                                            else:
                                                                matrix_C[ind_o_cell][ind_o_cell] += Fcc[iLayer]
                                                            if sparseM==1:
                                                                matrix_C2[iLayer*Ntot+j,iLayer*Ntot+ind_o_cell] -= Fcc[iLayer]
                                                            else:
                                                                matrix_C[j][ind_o_cell] -= Fcc[iLayer]
                                                    else: #From j to ind_o_cell
                                                            if sparseM==1:
                                                                matrix_C2[iLayer*Ntot+j,iLayer*Ntot+j] += Fcc[iLayer]
                                                            else:
                                                                matrix_C[j][j] += Fcc[iLayer]
                                                            if sparseM==1:
                                                                matrix_C2[iLayer*Ntot+ind_o_cell,iLayer*Ntot+j] -= Fcc[iLayer]
                                                            else:
                                                                matrix_C[ind_o_cell][j] -= Fcc[iLayer]
                                        else: #Only Sym contagion
                                            if Fcc>0: #Flow from high index to low index cell
                                                if ind_o_cell<j: #From j to ind_o_cell
                                                    if j not in N_Dirichlet_ini:
                                                        matrix_SymC[j-NwallsJun][j-NwallsJun] -= Fcc
                                                    if ind_o_cell not in N_Dirichlet_ini:
                                                        matrix_SymC[ind_o_cell-NwallsJun][j-NwallsJun] += Fcc
                                                else: #From ind_o_cell to j
                                                    if ind_o_cell not in N_Dirichlet_ini:
                                                        matrix_SymC[ind_o_cell-NwallsJun][ind_o_cell-NwallsJun] -= Fcc
                                                    if j not in N_Dirichlet_ini:
                                                        matrix_SymC[j-NwallsJun][ind_o_cell-NwallsJun] += Fcc
                                            else: #Flow from low index to high index cell
                                                if ind_o_cell<j: #From ind_o_cell to j
                                                    if ind_o_cell not in N_Dirichlet_ini:
                                                        matrix_SymC[ind_o_cell-NwallsJun][ind_o_cell-NwallsJun] += Fcc
                                                    if j not in N_Dirichlet_ini:
                                                        matrix_SymC[j-NwallsJun][ind_o_cell-NwallsJun] -= Fcc
                                                else: #From j to ind_o_cell
                                                    if j not in N_Dirichlet_ini:
                                                        matrix_SymC[j-NwallsJun][j-NwallsJun] += Fcc
                                                    if ind_o_cell not in N_Dirichlet_ini:
                                                        matrix_SymC[ind_o_cell-NwallsJun][j-NwallsJun] -= Fcc
                                    
                                    if Sym_Contagion==1:
                                        itemp=0
                                        while not Cell_connec[ind_o_cell-NwallsJun][itemp] == j-NwallsJun:
                                            itemp+=1
                                        Cell_connec_flow[ind_o_cell-NwallsJun][itemp]=sign(temp)
                                        itemp=0
                                        while not Cell_connec[j-NwallsJun][itemp] == ind_o_cell-NwallsJun:
                                            itemp+=1
                                        Cell_connec_flow[j-NwallsJun][itemp]=-sign(temp)
                        elif Paraview==0 and ParTrack==0:
                            if path == "membrane": #Membrane connection
                                #Discretization based on cell layers and apoplasmic barriers
                                rank = int(Cell_rank[j-NwallsJun])
                                row = int(rank2row[rank])
                                if rank == 1 and count_epi > 0: #Outer exodermis
                                    row += 1
                                if rank == 3 and count_cortex > 0: #Outer endodermis
                                    if any(passage_cell_ID==array(j-NwallsJun)):# and Barrier==2:
                                        row += 2
                                    else:
                                        row += 3
                                elif rank == 3 and count_stele_overall > 0: #Inner endodermis
                                    if any(passage_cell_ID==array(j-NwallsJun)):# and Barrier==2:
                                        row += 1
                                if PileUp==2:
                                    Fcw=empty((AxialLayers))
                                    iLayer=0
                                    if Ana_Os:
                                        for iM in range(Nmaturity):
                                            Fcw[iLayer:iLayer+int(Nlayers[iM])] = Kmb[jmb,iM] * (psi[iLayer:iLayer+int(Nlayers[iM])] - psin[iLayer:iLayer+int(Nlayers[iM])] + s_membranes[jmb,iM]*(Os_walls[i,iLayer:iLayer+int(Nlayers[iM])] - Os_cells[j-NwallsJun,iM])) #(cm^3/d) Water flow rate positive from wall to protoplast
                                            if ((j-NwallsJun not in InterCid) and (j not in listxyl)) or Barriers[iM]==0: #No aerenchyma in the elongation zone
                                                for iL in range(iLayer,int(iLayer+Nlayers[iM])):
                                                    if Fcw[iL] > 0 :
                                                        UptakeLayer_plus[row][iL][count] += Fcw[iL] #grouping membrane flow rates in cell layers
                                                    else:
                                                        UptakeLayer_minus[row][iL][count] += Fcw[iL]
                                            if Kmb[jmb,iM]>1.0e-18: #Not an impermeable wall
                                                PsiWallLayer[row,iLayer:int(iLayer+Nlayers[iM]),count] += psi[iLayer:int(iLayer+Nlayers[iM])]
                                                NWallLayer[row,iMaturity,count] += 1
                                            iLayer+=int(Nlayers[iM])
                                    else:
                                        for iM in range(Nmaturity):
                                            Fcw[iLayer:iLayer+int(Nlayers[iM])] = Kmb[jmb,iM] * (psi[iLayer:iLayer+int(Nlayers[iM])] - psin[iLayer:iLayer+int(Nlayers[iM])] + s_membranes[jmb,iM]*(Os_walls[i,iM] - Os_cells[j-NwallsJun,iM])) #(cm^3/d) Water flow rate positive from wall to protoplast
                                            if ((j-NwallsJun not in InterCid) and (j not in listxyl)) or Barriers[iM]==0: #No aerenchyma in the elongation zone
                                                for iL in range(iLayer,int(iLayer+Nlayers[iM])):
                                                    if Fcw[iL] > 0 :
                                                        UptakeLayer_plus[row][iL][count] += Fcw[iL] #grouping membrane flow rates in cell layers
                                                    else:
                                                        UptakeLayer_minus[row][iL][count] += Fcw[iL]
                                            if Kmb[jmb,iM]>1.0e-18: #Not an impermeable wall
                                                PsiWallLayer[row,iLayer:int(iLayer+Nlayers[iM]),count] += psi[iLayer:int(iLayer+Nlayers[iM])]
                                                NWallLayer[row,iMaturity,count] += 1
                                            iLayer+=int(Nlayers[iM])
                                else:
                                    K=Kmb[jmb]
                                    #Flow densities calculation
                                    #Macroscopic distributed parameter for transmembrane flow
                                    if AxialLayers>1 and Ana_Os:
                                        Flow = K * (psi - psin + s_membranes[jmb]*(Os_walls[i,:] - Os_cells[j-NwallsJun,0]))
                                    else:
                                        Flow = K * (psi - psin + s_membranes[jmb]*(Os_walls[i,0] - Os_cells[j-NwallsJun,0]))
                                    
                                    if ((j-NwallsJun not in InterCid) and (j not in listxyl)) or Barrier==0:
                                        for iLayer in range(AxialLayers):
                                            iL=TopLayer-AxialLayers+iLayer
                                            if Flow[iLayer] > 0 :
                                                UptakeLayer_plus[row][iL][count] += Flow[iLayer] #grouping membrane flow rates in cell layers
                                            else:
                                                UptakeLayer_minus[row][iL][count] += Flow[iLayer]
                                    
                                    if K>1.0e-18: #Not an impermeable wall
                                        for iLayer in range(AxialLayers):
                                            iL=TopLayer-AxialLayers+iLayer
                                            PsiWallLayer[row,iL,count] += psi[iLayer]
                                        NWallLayer[row,iMaturity,count] += 1
                                jmb+=1
            
            #Axial solute convection
            if Sym_Contagion==2 and Apo_Contagion==2:
                if PileUp==2:
                    for i in range(NwallsJun):
                        Faxial_sym_apo_list.append((K_axial[i:-Ntot:Ntot].T*(soln[i:(AxialLayers-1)*Ntot:Ntot]-soln[i+Ntot:AxialLayers*Ntot:Ntot])).T) #Axial flow in walls and xylem vessels (cm3/d, positive upwards)
                    for i in range(NwallsJun,Ntot):
                        if i in listprotosieve or (not ((i in listxyl) or (i in listsieve))):
                            #The osmotic potential is currently uniform axially except from immature to mature xylem and metaphloem, so its effect cancels out in other cell types
                            #print(soln[i:(AxialLayers-1)*Ntot:Ntot],soln[i+Ntot:AxialLayers*Ntot:Ntot])
                            temp=empty((AxialLayers-1,1))
                            temp2=empty((AxialLayers-1,1))
                            iLayer=0
                            if Nmaturity>1:
                                for iM in range(Nmaturity-1):
                                    temp[iLayer:int(iLayer+Nlayers[iM])]=(K_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot].T*Fraction_notTM[i-NwallsJun,iM]*(soln[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot]-soln[(iLayer+1)*Ntot+i:int(iLayer+1+Nlayers[iM])*Ntot:Ntot])).T
                                    temp2[iLayer:int(iLayer+Nlayers[iM])]=(K_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot].T*(1-Fraction_notTM[i-NwallsJun,iM])*(soln[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot]-soln[(iLayer+1)*Ntot+i:int(iLayer+1+Nlayers[iM])*Ntot:Ntot])).T
                                    iLayer+=int(Nlayers[iM])
                                iM+=1
                            else:
                                iM=0
                            #Last maturity zone
                            temp[iLayer:int(iLayer+Nlayers[iM]-1)]=(K_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM]-1)*Ntot:Ntot].T*Fraction_notTM[i-NwallsJun,iM]*(soln[iLayer*Ntot+i:int(iLayer+Nlayers[iM]-1)*Ntot:Ntot]-soln[(iLayer+1)*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot])).T
                            temp2[iLayer:int(iLayer+Nlayers[iM]-1)]=(K_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM]-1)*Ntot:Ntot].T*(1-Fraction_notTM[i-NwallsJun,iM])*(soln[iLayer*Ntot+i:int(iLayer+Nlayers[iM]-1)*Ntot:Ntot]-soln[(iLayer+1)*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot])).T
                            Faxial_sym_apo_list.append(temp)
                            Faxial_TM_list.append(temp2)
                            #Faxial_sym_apo_list.append((K_axial[i:-Ntot:Ntot].T*Fraction_notTM[i-NwallsJun,MaturL]*(soln[i:(AxialLayers-1)*Ntot:Ntot]-soln[i+Ntot:AxialLayers*Ntot:Ntot])).T)
                            #Faxial_TM_list.append((K_axial[i:-Ntot:Ntot].T*(1-Fraction_notTM[i-NwallsJun,iM])*(soln[i:(AxialLayers-1)*Ntot:Ntot]-soln[i+Ntot:AxialLayers*Ntot:Ntot])).T)
                        else: #Metaphloem and xylem, mature or not (both have a transition between levels of maturity 0 and >0)
                            temp=empty((AxialLayers-1,1))
                            temp2=empty((AxialLayers-1,1))
                            iLayer=0
                            if Nmaturity>1:
                                for iM in range(Nmaturity-1):
                                    if Barriers[iM]==0 and Barriers[iM+1]>0: #transition between immature and mature vascular elements
                                        temp[iLayer:int(iLayer+Nlayers[iM]-1)]=(K_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM]-1)*Ntot:Ntot].T*Fraction_notTM[i-NwallsJun,iM]*(soln[iLayer*Ntot+i:int(iLayer+Nlayers[iM]-1)*Ntot:Ntot]-soln[(iLayer+1)*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot])).T
                                        temp2[iLayer:int(iLayer+Nlayers[iM]-1)]=(K_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM]-1)*Ntot:Ntot].T*(1-Fraction_notTM[i-NwallsJun,iM])*(soln[iLayer*Ntot+i:int(iLayer+Nlayers[iM]-1)*Ntot:Ntot]-soln[(iLayer+1)*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot])).T
                                        if i in listsieve:
                                            #At this transition, symplastic and transmembrane fluxes might have opposite directions. Here we save the net axial flow rate
                                            temp[int(iLayer+Nlayers[iM]-1)]=K_axial[int(iLayer+Nlayers[iM]-1)*Ntot+i]*Fraction_notTM[i-NwallsJun,iM]*(soln[int(iLayer+Nlayers[iM]-1)*Ntot+i]-soln[int(iLayer+Nlayers[iM])*Ntot+i])
                                            temp2[int(iLayer+Nlayers[iM]-1)]=K_axial[int(iLayer+Nlayers[iM]-1)*Ntot+i]*(1-Fraction_notTM[i-NwallsJun,iM])*(soln[int(iLayer+Nlayers[iM]-1)*Ntot+i]-soln[int(iLayer+Nlayers[iM])*Ntot+i]+s_sieve*(Os_cells[i-NwallsJun,iM] - Os_cells[i-NwallsJun,iM+1]))
                                        else: #i in listxyl
                                            #At this transition, only the transmembrane pathway exists
                                            temp[int(iLayer+Nlayers[iM]-1)]=0.0
                                            temp2[int(iLayer+Nlayers[iM]-1)]=K_axial[int(iLayer+Nlayers[iM]-1)*Ntot+i]*(soln[int(iLayer+Nlayers[iM]-1)*Ntot+i]-soln[int(iLayer+Nlayers[iM])*Ntot+i]+s_stele*(Os_cells[i-NwallsJun,iM] - Os_cells[i-NwallsJun,iM+1]))
                                    else:#The osmotic potential is currently uniform axially except at transitions between immature and mature xylem or metaphloem, so here its effect cancels out
                                        temp[iLayer:int(iLayer+Nlayers[iM])]=(K_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot].T*Fraction_notTM[i-NwallsJun,iM]*(soln[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot]-soln[(iLayer+1)*Ntot+i:int(iLayer+1+Nlayers[iM])*Ntot:Ntot])).T
                                        temp2[iLayer:int(iLayer+Nlayers[iM])]=(K_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot].T*(1-Fraction_notTM[i-NwallsJun,iM])*(soln[iLayer*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot]-soln[(iLayer+1)*Ntot+i:int(iLayer+1+Nlayers[iM])*Ntot:Ntot])).T
                                    iLayer+=int(Nlayers[iM])
                                #iM = Nmaturity, shorter component (does not include the last iLayer)
                                iM+=1
                            else:
                                iM=0
                            #Last maturity zone
                            if Barriers[iM]>0: #Mature xylem vessel or metaphloem
                                temp[iLayer:int(iLayer+Nlayers[iM]-1)]=(K_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM]-1)*Ntot:Ntot].T*Fraction_notTM[i-NwallsJun,iM]*(soln[iLayer*Ntot+i:int(iLayer+Nlayers[iM]-1)*Ntot:Ntot]-soln[(iLayer+1)*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot])).T
                                temp2[iLayer:int(iLayer+Nlayers[iM]-1)]=(K_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM]-1)*Ntot:Ntot].T*(1-Fraction_notTM[i-NwallsJun,iM])*(soln[iLayer*Ntot+i:int(iLayer+Nlayers[iM]-1)*Ntot:Ntot]-soln[(iLayer+1)*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot])).T
                            else:
                                #The osmotic potential is currently uniform axially, so its effect cancels out
                                #print(soln[i:(AxialLayers-1)*Ntot:Ntot],soln[i+Ntot:AxialLayers*Ntot:Ntot])
                                temp[iLayer:int(iLayer+Nlayers[iM]-1)]=(K_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM]-1)*Ntot:Ntot].T*Fraction_notTM[i-NwallsJun,iM]*(soln[iLayer*Ntot+i:int(iLayer+Nlayers[iM]-1)*Ntot:Ntot]-soln[(iLayer+1)*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot])).T
                                temp2[iLayer:int(iLayer+Nlayers[iM]-1)]=(K_axial[iLayer*Ntot+i:int(iLayer+Nlayers[iM]-1)*Ntot:Ntot].T*(1-Fraction_notTM[i-NwallsJun,iM])*(soln[iLayer*Ntot+i:int(iLayer+Nlayers[iM]-1)*Ntot:Ntot]-soln[(iLayer+1)*Ntot+i:int(iLayer+Nlayers[iM])*Ntot:Ntot])).T
                            Faxial_sym_apo_list.append(temp)
                            Faxial_TM_list.append(temp2)
                else:
                    for i in range(Ntot):
                        if i<NwallsJun:
                            Faxial_sym_apo_list.append(K_axial[i]*(soln[i:(AxialLayers-1)*Ntot:Ntot]-soln[i+Ntot:AxialLayers*Ntot:Ntot])) #Axial flow in walls and xylem vessels (cm3/d, positive upwards)
                        elif Barrier>0 and (i in listxyl): #Mature xylem vessel
                            Faxial_sym_apo_list.append(K_axial[i]*Fraction_notTM[i-NwallsJun]*(soln[i:(AxialLayers-1)*Ntot:Ntot]-soln[i+Ntot:AxialLayers*Ntot:Ntot]))
                            Faxial_TM_list.append(K_axial[i]*(1-Fraction_notTM[i-NwallsJun])*(soln[i:(AxialLayers-1)*Ntot:Ntot]-soln[i+Ntot:AxialLayers*Ntot:Ntot]))
                        else:
                            #The osmotic potential is currently uniform axially, so its effect cancels out
                            #print(soln[i:(AxialLayers-1)*Ntot:Ntot],soln[i+Ntot:AxialLayers*Ntot:Ntot])
                            Faxial_sym_apo_list.append(K_axial[i]*Fraction_notTM[i-NwallsJun]*(soln[i:(AxialLayers-1)*Ntot:Ntot]-soln[i+Ntot:AxialLayers*Ntot:Ntot]))
                            Faxial_TM_list.append(K_axial[i]*(1-Fraction_notTM[i-NwallsJun])*(soln[i:(AxialLayers-1)*Ntot:Ntot]-soln[i+Ntot:AxialLayers*Ntot:Ntot]))
                #for cid in listxyl:
                #    Faxial_list_BC.append(K_axial[i]*(soln[cid+(AxialLayers-1)*Ntot]-Psi_xyl[1][iMaturity][count]))
                
                i=0 #node index (0 to Ntot-1)
                for Faxial in Faxial_sym_apo_list:
                    for iLayer in range(AxialLayers-1):
                        Q=Faxial[iLayer] #Positive when pressure i (distal) > pressure i+Ntot (proximal)
                        if Q>0: #Flow from i upwards
                            if sparseM==1:
                                    matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] -= Q
                                    matrix_C2[(iLayer+1)*Ntot+i,iLayer*Ntot+i] += Q
                            else:
                                    matrix_C[iLayer*Ntot+i][iLayer*Ntot+i] -= Q
                                    matrix_C[(iLayer+1)*Ntot+i][iLayer*Ntot+i] += Q
                        else: #Flow from i+Ntot downwards
                            if sparseM==1:
                                    matrix_C2[(iLayer+1)*Ntot+i,(iLayer+1)*Ntot+i] += Q  #changed sign 12/2/2019
                                    matrix_C2[iLayer*Ntot+i,(iLayer+1)*Ntot+i] -= Q #changed sign 12/2/2019
                            else:
                                    matrix_C[(iLayer+1)*Ntot+i][(iLayer+1)*Ntot+i] += Q  #changed sign 12/2/2019
                                    matrix_C[iLayer*Ntot+i][(iLayer+1)*Ntot+i] -= Q #changed sign 12/2/2019
                    i+=1
                
                i=NwallsJun #cell node index (NwallsJun to Ntot-1)
                for Faxial in Faxial_TM_list:
                    for iLayer in range(AxialLayers-1):
                        Q=Faxial[iLayer] #Positive when pressure i (distal) > pressure i+Ntot (proximal)
                        if Q>0: #Flow from i upwards
                            if sparseM==1:
                                if D2O1==1:#Solute that moves across membranes like water 
                                    matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] -= Q
                                    matrix_C2[(iLayer+1)*Ntot+i,iLayer*Ntot+i] += Q
                                else:#Solute that moves across membranes independently of water (the membrane is possibly not one) 
                                    if PileUp==2:
                                        matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] -= Q*(1-s_cells[i-NwallsJun,MaturL[iLayer]])
                                        matrix_C2[(iLayer+1)*Ntot+i,iLayer*Ntot+i] += Q*(1-s_cells[i-NwallsJun,MaturL[iLayer]])
                                    else:
                                        matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] -= Q*(1-s_cells[i-NwallsJun])
                                        matrix_C2[(iLayer+1)*Ntot+i,iLayer*Ntot+i] += Q*(1-s_cells[i-NwallsJun])
                            else:
                                    matrix_C[iLayer*Ntot+i][iLayer*Ntot+i] -= Q
                                    matrix_C[(iLayer+1)*Ntot+i][iLayer*Ntot+i] += Q
                        else: #Flow from i+Ntot downwards
                            if sparseM==1:
                                if D2O1==1:#Solute that moves across membranes like water 
                                    matrix_C2[(iLayer+1)*Ntot+i,(iLayer+1)*Ntot+i] += Q  #changed sign 12/2/2019
                                    matrix_C2[iLayer*Ntot+i,(iLayer+1)*Ntot+i] -= Q #changed sign 12/2/2019
                                else:
                                    if PileUp==2:
                                        matrix_C2[(iLayer+1)*Ntot+i,(iLayer+1)*Ntot+i] += Q*(1-s_cells[i-NwallsJun,MaturL[iLayer]])  #changed sign 12/2/2019
                                        matrix_C2[iLayer*Ntot+i,(iLayer+1)*Ntot+i] -= Q*(1-s_cells[i-NwallsJun,MaturL[iLayer]]) #changed sign 12/2/2019
                                    else:
                                        matrix_C2[(iLayer+1)*Ntot+i,(iLayer+1)*Ntot+i] += Q*(1-s_cells[i-NwallsJun])  #changed sign 12/2/2019
                                        matrix_C2[iLayer*Ntot+i,(iLayer+1)*Ntot+i] -= Q*(1-s_cells[i-NwallsJun]) #changed sign 12/2/2019
                            else:
                                    matrix_C[(iLayer+1)*Ntot+i][(iLayer+1)*Ntot+i] += Q  #changed sign 12/2/2019
                                    matrix_C[iLayer*Ntot+i][(iLayer+1)*Ntot+i] -= Q #changed sign 12/2/2019
                    i+=1
            
            ###############################################
            #Setting solute BC (AxialNeumann, SoilNeumann)#
            ###############################################
            if Sym_Contagion==2 and Apo_Contagion==2:
                if test_mass_balance==1:
                    for i in range(AxialLayers*Ntot):
                        temp=sum(matrix_C2[i])
                        print('mass_balance pre-BC',i,temp)
                
                rhs_C = np.zeros((AxialLayers*Ntot),dtype=dp) #Initializing the right-hand side matrix of solute apoplastic concentrations
                
                #Solute production BC type
                if PileUp==2:
                    for cid in N_Prod_ini:
                        rhs_C[cid] -= cc_Prod_ini[N_Prod_ini.index(cid)]*Volume_W2[cid] #units: mol/day
                else:
                    for cid in N_Prod_ini:
                        rhs_C[cid] -= cc_Prod_ini[N_Prod_ini.index(cid)]*Volume_W[cid] #units: mol/day
                
                #AxialNeumann BC type
                if XylOK>0:
                    i=0
                    for cid in listxyl:
                        Q_dist=float(Q_xyl[0][i]) #Positive upwards
                        Q_prox=float(Q_xyl[1][i])
                        
                        if Q_prox>0: #Flow from proximal xylem vessel toward shoot
                            if sparseM==1:
                                    matrix_C2[-Ntot+cid,-Ntot+cid] -= Q_prox
                            else:
                                    matrix_C[-Ntot+cid][-Ntot+cid] -= Q_prox
                        else: #Flow from shoot towards proximal xylem vessel
                            if cid-NwallsJun in N_AxialNeumann_ini: #if cid not in Nodes_C_Dirichlet_ini:
                                if PileUp==1:
                                    print('No proximal solute flux BC when using the PileUp mode')
                                    #rhs_C[-Ntot+cid][0] += Q_prox*Sym_cc_flows_ini[h][iMaturity][count-1][i][1]
                                else:
                                    rhs_C[-Ntot+cid] += Q_prox*cc_AxialNeumann_ini[N_AxialNeumann_ini.index(cid-NwallsJun)][1]
                        if Q_dist<0: #Flow from distal xylem vessel toward root tip
                            if sparseM==1:
                                    matrix_C2[cid,cid] += Q_dist
                            else:
                                    matrix_C[cid][cid] += Q_dist
                        else: #Flow from root tip towards distal xylem vessel
                            if cid-NwallsJun in N_AxialNeumann_ini: #if cid not in Nodes_C_Dirichlet_ini:
                                if PileUp==1:
                                    rhs_C[cid] -= Q_dist*cc_AxialNeumann_ini[h][iMaturity][count-1][N_AxialNeumann_ini.index(cid-NwallsJun)][0]
                                else:
                                    rhs_C[cid] -= Q_dist*cc_AxialNeumann_ini[N_AxialNeumann_ini.index(cid-NwallsJun)][0] #possibly use of "" instead of i in case more than xylem vessels taken into account
                        i+=1
                
                #Phloem BC for solute
                i=0
                for cid in listsieve:
                    if Q_sieve[1][i]>0: #Proximal phloem BC shootward
                        if sparseM==1:
                            matrix_C2[-Ntot+cid,-Ntot+cid] -= Q_sieve[1][i]
                        else:
                            matrix_C[-Ntot+cid][-Ntot+cid] -= Q_sieve[1][i]
                    else: #Flow from shoot towards proximal phloem sieve tube
                        if cid-NwallsJun in N_AxialNeumann_ini: #if cid not in Nodes_C_Dirichlet_ini:
                            if PileUp==1:
                                print('No proximal solute flux BC when using the PileUp==1 mode')
                                #rhs_C[-Ntot+cid][0] += Q_prox*Sym_cc_flows_ini[h][iMaturity][count-1][i][1]
                            else:
                                rhs_C[-Ntot+cid] += Q_sieve[1][i]*cc_AxialNeumann_ini[N_AxialNeumann_ini.index(cid-NwallsJun)][1]
                    if Q_sieve[0][i]<0: #Distal phloem BC tipward
                        if sparseM==1:
                            matrix_C2[cid,cid] += Q_sieve[0][i]
                        else:
                            matrix_C[cid][cid] += Q_sieve[0][i]
                    else: #Flow from tip towards distal phloem sieve tube
                        if cid-NwallsJun in N_AxialNeumann_ini: #if cid not in Nodes_C_Dirichlet_ini:
                            if PileUp==1:
                                rhs_C[cid] -= Q_sieve[0][i]*cc_AxialNeumann_ini[h][iMaturity][count-1][N_AxialNeumann_ini.index(cid-NwallsJun)][0]
                            else:
                                rhs_C[cid] -= Q_sieve[0][i]*cc_AxialNeumann_ini[N_AxialNeumann_ini.index(cid-NwallsJun)][0]
                    i+=1
                
                #ElongNeumann BC type, accounting for "elongation sink"
                if Apo_Contagion==2 and Sym_Contagion==2:
                    if sparseM==1:
                        for iLayer in range(AxialLayers):
                            matrix_C2[iLayer*Ntot:(iLayer+1)*Ntot,iLayer*Ntot:(iLayer+1)*Ntot] -= spdiags(rhs_e[iLayer*Ntot:(iLayer+1)*Ntot].T, 0, Ntot, Ntot).toarray()
                    else:
                        matrix_C -= diag(rhs_e[:,0])
                
                #SoilNeumann BC type
                if Apo_Contagion==2:
                    if Sym_Contagion==2:
                        if PileUp==2:
                            i=0
                            for ind in Borderwall:
                                for iLayer in range(AxialLayers):
                                    if ind+Ntot*iLayer in N_SoilNeumann_ini:
                                        #Convection
                                        Q=float(Q_soil[i])
                                        if Q<0.0: #Solute flowing out of the root
                                            matrix_C2[iLayer*Ntot+ind,iLayer*Ntot+ind]+=Q
                                        else: #Solute flowing into the root
                                            rhs_C[iLayer*Ntot+ind]-=Q*cc_SoilNeumann_ini[N_SoilNeumann_ini.index(iLayer*Ntot+ind)]
                                        #Diffusion
                                        if Dif_BC:
                                            temp=1.0E-04*heights[iLayer]*lengths[ind]/(thickness/2) #Wall surface to thickness ratio (cm)
                                            DF=temp*Diff_PW1 #"Diffusive flux" (cm^3/d) temp is the section to length ratio of the wall to junction path
                                            matrix_C2[iLayer*Ntot+ind,iLayer*Ntot+ind] -= DF
                                            rhs_C[iLayer*Ntot+ind] -= DF*cc_SoilNeumann_ini[N_SoilNeumann_ini.index(iLayer*Ntot+ind)]
                                    i+=1
                            for ind in Borderjunction:
                                for iLayer in range(AxialLayers):
                                    if ind+Ntot*iLayer in N_SoilNeumann_ini:
                                        #Convection
                                        Q=float(Q_soil[i])
                                        if Q<0.0: #Solute flowing out of the root
                                            matrix_C2[iLayer*Ntot+ind,iLayer*Ntot+ind]+=Q
                                        else: #Solute flowing into the root
                                            rhs_C[iLayer*Ntot+ind]-=Q*cc_SoilNeumann_ini[N_SoilNeumann_ini.index(iLayer*Ntot+ind)]
                                        #Diffusion
                                        if Dif_BC:
                                            temp=1.0E-04*heights[iLayer]*lengths[ind]/(thickness/2) #Wall surface to thickness ratio (cm)
                                            DF=temp*Diff_PW1 #"Diffusive flux" (cm^3/d) temp is the section to length ratio of the wall to junction path
                                            matrix_C2[iLayer*Ntot+ind,iLayer*Ntot+ind] -= DF
                                            rhs_C[iLayer*Ntot+ind] -= DF*cc_SoilNeumann_ini[N_SoilNeumann_ini.index(iLayer*Ntot+ind)]
                                    i+=1
                        else:
                            i=0
                            for ind in Borderwall:
                                for iLayer in range(AxialLayers):
                                    if ind+Ntot*iLayer in N_SoilNeumann_ini:
                                        #Convection
                                        Q=float(Q_soil[i])
                                        if Q<0.0: #Solute flowing out of the root
                                            if sparseM==1:
                                                matrix_C2[iLayer*Ntot+ind,iLayer*Ntot+ind]+=Q
                                            else:
                                                matrix_C[ind][ind]+=Q
                                        else: #Solute flowing into the root
                                            if sparseM==1:
                                                rhs_C[iLayer*Ntot+ind]-=Q*cc_SoilNeumann_ini[N_SoilNeumann_ini.index(iLayer*Ntot+ind)]
                                            else:
                                                rhs_C[ind]-=Q*cc_SoilNeumann_ini[N_SoilNeumann_ini.index(ind)]
                                        #Diffusion
                                        if Dif_BC:
                                            temp=1.0E-04*height*lengths[ind]/(thickness/2) #Wall surface to thickness ratio (cm)
                                            DF=temp*Diff_PW1 #"Diffusive flux" (cm^3/d) temp is the section to length ratio of the wall to junction path
                                            if sparseM==1:
                                                matrix_C2[iLayer*Ntot+ind,iLayer*Ntot+ind] -= DF
                                                rhs_C[iLayer*Ntot+ind] -= DF*cc_SoilNeumann_ini[N_SoilNeumann_ini.index(iLayer*Ntot+ind)]
                                            else:
                                                matrix_C[ind][ind] -= DF
                                                rhs_C[ind] -= DF*cc_SoilNeumann_ini[N_SoilNeumann_ini.index(ind)]
                                    i+=1
                            for ind in Borderjunction:
                                for iLayer in range(AxialLayers):
                                    if ind+Ntot*iLayer in N_SoilNeumann_ini:
                                        #Convection
                                        Q=float(Q_soil[i])
                                        if Q<0.0: #Solute flowing out of the root
                                            if sparseM==1:
                                                matrix_C2[iLayer*Ntot+ind,iLayer*Ntot+ind]+=Q
                                            else:
                                                matrix_C[ind][ind]+=Q
                                        else: #Solute flowing into the root
                                            if sparseM==1:
                                                rhs_C[iLayer*Ntot+ind]-=Q*cc_SoilNeumann_ini[N_SoilNeumann_ini.index(iLayer*Ntot+ind)]
                                            else:
                                                rhs_C[ind]-=Q*cc_SoilNeumann_ini[N_SoilNeumann_ini.index(ind)]
                                        #Diffusion
                                        if Dif_BC:
                                            temp=1.0E-04*height*lengths[ind]/(thickness/2) #Wall surface to thickness ratio (cm)
                                            DF=temp*Diff_PW1 #"Diffusive flux" (cm^3/d) temp is the section to length ratio of the wall to junction path
                                            if sparseM==1:
                                                matrix_C2[iLayer*Ntot+ind,iLayer*Ntot+ind] -= DF
                                                rhs_C[iLayer*Ntot+ind] -= DF*cc_SoilNeumann_ini[N_SoilNeumann_ini.index(iLayer*Ntot+ind)]
                                            else:
                                                matrix_C[ind][ind] -= DF
                                                rhs_C[i] -= DF*cc_SoilNeumann_ini[N_SoilNeumann_ini.index(ind)]
                                    i+=1
                    else:
                        i=0
                        for ind in Borderwall:
                            for iLayer in range(AxialLayers):
                                if ind+Ntot*iLayer in N_SoilNeumann_ini:
                                    Q=float(Q_soil[i])
                                    if Q<0.0:
                                        matrix_ApoC[ind][ind]+=Q
                                    else: #Solute flowing into the root
                                        rhs_ApoC[ind]-=Q*cc_SoilNeumann_ini[N_SoilNeumann_ini.index(ind)]
                                    #Diffusion
                                    if Dif_BC:
                                        temp=1.0E-04*height*lengths[ind]/(thickness/2) #Wall surface to thickness ratio (cm)
                                        DF=temp*Diff_PW1 #"Diffusive flux" (cm^3/d) temp is the section to length ratio of the wall to junction path
                                        matrix_ApoC[ind][ind] -= DF
                                        rhs_ApoC[i] -= DF*cc_SoilNeumann_ini[N_SoilNeumann_ini.index(ind)]
                                i+=1
                        for ind in Borderjunction:
                            for iLayer in range(AxialLayers):
                                if ind+Ntot*iLayer in N_SoilNeumann_ini:
                                    Q=float(Q_soil[i])
                                    if Q<0.0:
                                        matrix_ApoC[ind][ind]+=Q
                                    else: #Solute flowing into the root
                                        rhs_ApoC[ind]-=Q*cc_SoilNeumann_ini[N_SoilNeumann_ini.index(ind)]
                                    #Diffusion
                                    if Dif_BC:
                                        temp=1.0E-04*height*lengths[ind]/(thickness/2) #Wall surface to thickness ratio (cm)
                                        DF=temp*Diff_PW1 #"Diffusive flux" (cm^3/d) temp is the section to length ratio of the wall to junction path
                                        matrix_ApoC[ind][ind] -= DF
                                        rhs_ApoC[i] -= DF*cc_SoilNeumann_ini[N_SoilNeumann_ini.index(ind)]
                                i+=1
                
                #Transient BC
                if Transient_tracer:
                    if PileUp==2: #then rhs_C_transi is the same for all scenarios
                        rhs_C_transi = np.zeros((AxialLayers*Ntot),dtype=dp)
                        
                        #Solute production BC type
                        for cid in N_Prod_transi:
                            rhs_C_transi[cid] -= cc_Prod_transi[N_Prod_transi.index(cid)]*Volume_W2[cid] #units: mol/day
                                                
                        #Neumann solute BC
                        if XylOK>0:
                            i=0
                            for cid in listxyl:
                                Q_dist=float(Q_xyl[0][i]) #Positive upwards
                                Q_prox=float(Q_xyl[1][i])
                                if Q_prox<0: #Flow from shoot towards proximal xylem vessel
                                    if cid-NwallsJun in N_AxialNeumann_transi: #if cid not in Nodes_C_Dirichlet_ini:
                                        rhs_C_transi[-Ntot+cid] += Q_prox*cc_AxialNeumann_transi[N_AxialNeumann_transi.index(cid-NwallsJun)][1]
                                if Q_dist>0:#Flow from root tip towards distal xylem vessel
                                    if cid-NwallsJun in N_AxialNeumann_transi: #if cid not in Nodes_C_Dirichlet_ini:
                                        rhs_C_transi[cid] -= Q_dist*cc_AxialNeumann_transi[N_AxialNeumann_transi.index(cid-NwallsJun)][0] #possibly use of "" instead of i in case more than xylem vessels taken into account
                                i+=1

                        #Phloem BC for solute
                        i=0
                        for cid in listsieve:
                            if Q_sieve[1][i]<0: #Flow from shoot towards proximal phloem sieve tube
                                if cid-NwallsJun in N_AxialNeumann_transi: #if cid not in Nodes_C_Dirichlet_ini:
                                    if PileUp==1:
                                        print('No proximal solute flux BC when using the PileUp==1 mode')
                                        #rhs_C[-Ntot+cid][0] += Q_prox*Sym_cc_flows_ini[h][iMaturity][count-1][i][1]
                                    else:
                                        rhs_C_transi[-Ntot+cid] += Q_sieve[1][i]*cc_AxialNeumann_transi[N_AxialNeumann_transi.index(cid-NwallsJun)][1]
                            if Q_sieve[0][i]>0: #Flow from shoot towards proximal phloem sieve tube
                                if cid-NwallsJun in N_AxialNeumann_transi: #if cid not in Nodes_C_Dirichlet_ini:
                                    if PileUp==1:
                                        rhs_C_transi[cid] -= Q_sieve[0][i]*cc_AxialNeumann_transi[h][iMaturity][count-1][N_AxialNeumann_transi.index(cid-NwallsJun)][0]
                                    else:
                                        rhs_C_transi[cid] -= Q_sieve[0][i]*cc_AxialNeumann_transi[N_AxialNeumann_transi.index(cid-NwallsJun)][0]
                            i+=1
                        
                        i=0
                        for ind in Borderwall:
                            for iLayer in range(AxialLayers):
                                if ind+Ntot*iLayer in N_SoilNeumann_transi:
                                    #Convection
                                    Q=float(Q_soil[i])
                                    if Q>0.0: #Solute flowing into the root
                                        rhs_C_transi[iLayer*Ntot+ind]-=Q*cc_SoilNeumann_transi[N_SoilNeumann_transi.index(iLayer*Ntot+ind)]
                                    #Diffusion
                                    if Dif_BC:
                                        temp=1.0E-04*heights[iLayer]*lengths[ind]/(thickness/2) #Wall surface to thickness ratio (cm)
                                        DF=temp*Diff_PW1 #"Diffusive flux" (cm^3/d) temp is the section to length ratio of the wall to junction path
                                        rhs_C_transi[iLayer*Ntot+ind] -= DF*cc_SoilNeumann_transi[N_SoilNeumann_transi.index(iLayer*Ntot+ind)]
                                i+=1
                        for ind in Borderjunction:
                            for iLayer in range(AxialLayers):
                                if ind+Ntot*iLayer in N_SoilNeumann_transi:
                                    #Convection
                                    Q=float(Q_soil[i])
                                    if Q>0.0: #Solute flowing into the root
                                        rhs_C_transi[iLayer*Ntot+ind]-=Q*cc_SoilNeumann_transi[N_SoilNeumann_transi.index(iLayer*Ntot+ind)]
                                    #Diffusion
                                    if Dif_BC:
                                        temp=1.0E-04*heights[iLayer]*lengths[ind]/(thickness/2) #Wall surface to thickness ratio (cm)
                                        DF=temp*Diff_PW1 #"Diffusive flux" (cm^3/d) temp is the section to length ratio of the wall to junction path
                                        rhs_C_transi[iLayer*Ntot+ind] -= DF*cc_SoilNeumann_transi[N_SoilNeumann_transi.index(iLayer*Ntot+ind)]
                                i+=1
                                
                    elif PileUp==0: #then rhs_C_transi is the same for all scenarios
                        rhs_C_transi = np.zeros((AxialLayers*Ntot),dtype=dp)
                        
                        #Solute production BC type
                        for cid in N_Prod_transi:
                            rhs_C_transi[cid] -= cc_Prod_transi[N_Prod_transi.index(cid)]*Volume_W[cid] #units: mol/day
                        
                        #Neumann solute BC
                        if Barrier>0:
                            i=0
                            for cid in listxyl:
                                Q_dist=float(Q_xyl[0][i]) #Positive upwards
                                Q_prox=float(Q_xyl[1][i])
                                if Q_prox<0: #Flow from shoot towards proximal xylem vessel
                                    if cid-NwallsJun in N_AxialNeumann_transi: #if cid not in Nodes_C_Dirichlet_ini:
                                        rhs_C_transi[-Ntot+cid] += Q_prox*cc_AxialNeumann_transi[N_AxialNeumann_transi.index(cid-NwallsJun)][1]
                                if Q_dist>0:#Flow from root tip towards distal xylem vessel
                                    if cid-NwallsJun in N_AxialNeumann_transi: #if cid not in Nodes_C_Dirichlet_ini:
                                        rhs_C_transi[cid] -= Q_dist*cc_AxialNeumann_transi[N_AxialNeumann_transi.index(cid-NwallsJun)][0] #possibly use of "" instead of i in case more than xylem vessels taken into account
                                i+=1
                        
                        #Phloem BC for solute
                        i=0
                        for cid in listsieve:
                            if Q_sieve[1][i]<0: #Flow from shoot towards proximal phloem sieve tube
                                if cid-NwallsJun in N_AxialNeumann_transi: #if cid not in Nodes_C_Dirichlet_ini:
                                    if PileUp==1:
                                        print('No proximal solute flux BC when using the PileUp==1 mode')
                                        #rhs_C[-Ntot+cid][0] += Q_prox*Sym_cc_flows_ini[h][iMaturity][count-1][i][1]
                                    else:
                                        rhs_C_transi[-Ntot+cid] += Q_sieve[1][i]*cc_AxialNeumann_transi[N_AxialNeumann_transi.index(cid-NwallsJun)][1]
                            if Q_sieve[0][i]>0: #Flow from tip towards distal phloem sieve tube
                                if cid-NwallsJun in N_AxialNeumann_transi: #if cid not in Nodes_C_Dirichlet_ini:
                                    if PileUp==1:
                                        rhs_C_transi[cid] -= Q_sieve[0][i]*cc_AxialNeumann_transi[h][iMaturity][count-1][N_AxialNeumann_transi.index(cid-NwallsJun)][0]
                                    else:
                                        rhs_C_transi[cid] -= Q_sieve[0][i]*cc_AxialNeumann_transi[N_AxialNeumann_transi.index(cid-NwallsJun)][0]
                            i+=1
                        
                        i=0
                        for ind in Borderwall:
                            for iLayer in range(AxialLayers):
                                if ind+Ntot*iLayer in N_SoilNeumann_transi:
                                    #Convection
                                    Q=float(Q_soil[i])
                                    if Q>0.0: #Solute flowing into the root
                                        if sparseM==1:
                                            rhs_C_transi[iLayer*Ntot+ind]-=Q*cc_SoilNeumann_transi[N_SoilNeumann_transi.index(iLayer*Ntot+ind)]
                                        else:
                                            rhs_C_transi[ind]-=Q*cc_SoilNeumann_transi[N_SoilNeumann_transi.index(ind)]
                                    #Diffusion
                                    if Dif_BC:
                                        temp=1.0E-04*height*lengths[ind]/(thickness/2) #Wall surface to thickness ratio (cm)
                                        DF=temp*Diff_PW1 #"Diffusive flux" (cm^3/d) temp is the section to length ratio of the wall to junction path
                                        if sparseM==1:
                                            rhs_C_transi[iLayer*Ntot+ind] -= DF*cc_SoilNeumann_transi[N_SoilNeumann_transi.index(iLayer*Ntot+ind)]
                                        else:
                                            rhs_C_transi[ind] -= DF*cc_SoilNeumann_transi[N_SoilNeumann_transi.index(ind)]
                                i+=1
                        for ind in Borderjunction:
                            for iLayer in range(AxialLayers):
                                if ind+Ntot*iLayer in N_SoilNeumann_transi:
                                    #Convection
                                    Q=float(Q_soil[i])
                                    if Q>0.0: #Solute flowing into the root
                                        if sparseM==1:
                                            rhs_C_transi[iLayer*Ntot+ind]-=Q*cc_SoilNeumann_transi[N_SoilNeumann_transi.index(iLayer*Ntot+ind)]
                                        else:
                                            rhs_C_transi[ind]-=Q*cc_SoilNeumann_transi[N_SoilNeumann_transi.index(ind)]
                                    #Diffusion
                                    if Dif_BC:
                                        temp=1.0E-04*height*lengths[ind]/(thickness/2) #Wall surface to thickness ratio (cm)
                                        DF=temp*Diff_PW1 #"Diffusive flux" (cm^3/d) temp is the section to length ratio of the wall to junction path
                                        if sparseM==1:
                                            rhs_C_transi[iLayer*Ntot+ind] -= DF*cc_SoilNeumann_transi[N_SoilNeumann_transi.index(iLayer*Ntot+ind)]
                                        else:
                                            rhs_C_transi[ind] -= DF*cc_SoilNeumann_transi[N_SoilNeumann_transi.index(ind)]
                                i+=1
            
            t01 = time.perf_counter()
            print(t01-t00, 'seconds process time to list water fluxes and adjust matrix_C')
            
            ######################################################
            #Solving stationary & transient solute concentrations#
            ######################################################
            if Apo_Contagion==2 or Sym_Contagion==2: #Sym & Apo contagion
                if Apo_Contagion==2 and Sym_Contagion==2: #Sym & Apo contagion
                    if test_mass_balance==1:
                        for i in range(AxialLayers*Ntot):
                            temp=sum(matrix_C2[i])-rhs_C[i]
                            #matrix_C2[i,i]-=temp #Guarantees mass balance
                            print('mass_balance post BC',i,temp)#,'corrected',sum(matrix_C2[i])-rhs_C[i])
                    if not (Transient_tracer and Ini_cond==2):
                        #Setting stationary Dirichlet BC
                        Stored_matrix_C_ini=[]
                        i=0
                        for nid in N_Dirichlet_ini:
                            if sparseM==1:
                                Stored_matrix_C_ini.append(matrix_C2[nid,:].tocsr())
                                matrix_C2[nid,:]=zeros((1,size(matrix_C2,1))) #Deleting all row nid
                                matrix_C2[nid,nid]=1.0
                            else:
                                Stored_matrix_C_ini.append(matrix_C[nid,:])
                                matrix_C[nid,:]=zeros((1,size(matrix_C,1))) #Deleting all row nid
                                matrix_C[nid,nid]=1.0
                            rhs_C[nid]=cc_Dirichlet_ini[i] #N_Dirichlet_ini.index(nid) #1.0 #Concentration in source protoplasts defined in Geom
                            i+=1
                        
                        #Solving steady-state apoplastic & symplastic concentrations
                        print("Solving steady-state tracer flow & printing outputs")
                        t00 = time.perf_counter()
                        
                        soln_C = empty((Ntot*AxialLayers,1))
                        if sparseM==1:
                            if any(rhs_C): #There is a non-zero element in rhs_C
                                soln_C = spsolve(matrix_C2.tocsr(),rhs_C)
                            else:
                                soln_C = rhs_C
                            C_tot0 = float(np.dot(soln_C,Volume_W2))
                            if size(N_AxialNeumann_ini)>0:
                                Flux_div_C=-sum(rhs_C)-Degrad1*float(np.dot(Volume_W2,soln_C))
                                if not C_tot0==0.0:
                                    print('Steady-state tracer quantity:',C_tot0,'   Flux divergence (d^-1)',Flux_div_C,'   Flux balance error (%/s)',100*(Flux_div_C)/C_tot0/3600/24)
                                else:
                                    print('Steady-state tracer quantity:',C_tot0,'   Flux divergence (d^-1)',Flux_div_C)
                            else:
                                if len(N_Dirichlet_ini)<1000:
                                    Flux_div_C=-Degrad1*float(np.dot(Volume_W2,soln_C))
                                    for i_node in N_Dirichlet_ini:
                                        Flux_div_C+=float((sum(matrix_C2[:,i_node])-1.0)*rhs_C[i_node])
                                    if not C_tot0==0.0:
                                        print('Steady-state tracer quantity:',C_tot0,'   Flux divergence (d^-1)',Flux_div_C,'   Flux balance error (%/s)',100*(Flux_div_C)/C_tot0/3600/24)
                                    else:
                                        print('Steady-state tracer quantity:',C_tot0,'   Flux divergence (d^-1)',Flux_div_C)
                                else:
                                    print('Steady-state tracer quantity:',C_tot0)
                        else:
                            if any(rhs_C): #There is a non-zero element in rhs_C
                                soln_C = np.linalg.solve(matrix_C,rhs_C) #Solving the equation to get apoplastic relative concentrations
                            else:
                                soln_C = rhs_C
                            C_tot0 = float(np.dot(soln_C,Volume_W))
                            if size(N_AxialNeumann_ini)>0:
                                Flux_div_C=-sum(rhs_C)-Degrad1*float(np.dot(Volume_W2,soln_C))
                            else:
                                Flux_div_C=-Degrad1*float(np.dot(Volume_W2,soln_C))
                                for i_node in N_Dirichlet_ini:
                                    Flux_div_C+=(sum(matrix_C[:,i_node])-1.0)*rhs_C[i_node]
                                    error('Check column vector matrix_C[:,i_node]')
                            if not C_tot0==0.0:
                                print('Steady-state tracer quantity:',C_tot0,'   Flux divergence (d^-1)',Flux_div_C,'   Flux balance error (%/s)',100*(Flux_div_C)/C_tot0/3600/24)
                            else:
                                print('Steady-state tracer quantity:',C_tot0,'   Flux divergence (d^-1)',Flux_div_C)
                            
                        #Removing Dirichlet BC
                        i=0
                        for nid in N_Dirichlet_ini:
                            if sparseM==1:
                                matrix_C2[nid,:]=Stored_matrix_C_ini[i]
                            else:
                                matrix_C[nid,:]=Stored_matrix_C_ini[i]
                            i+=1
                    
                    ##Verification that computation was correct
                    #verif1=np.allclose(np.dot(matrix_C,soln_C),rhs_C)
                    
                    if Paraview==1 or ParTrack==1:
                        for Print_layer in Printrange:
                            iLayer=int(Print_layer.get("value"))
                            if iLayer<AxialLayers:
                                if PileUp==2:
                                    b=int(BarriersL[iLayer])
                                    he=heights[iLayer,0]
                                else:
                                    b=Barrier
                                    he=height
                                text_file = open(newpath+"Sym_Contagion_bottomb"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+"t00.pvtk", "w")
                                with open(newpath+"Sym_Contagion_bottomb"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+"t00.pvtk", "a") as myfile:
                                    myfile.write("# vtk DataFile Version 4.0 \n")
                                    myfile.write("Symplastic hormone concentration \n")
                                    myfile.write("ASCII \n")
                                    myfile.write(" \n")
                                    myfile.write("DATASET UNSTRUCTURED_GRID \n")
                                    myfile.write("POINTS "+str(len(ThickWalls))+" float \n")
                                    for ThickWallNode in ThickWalls:
                                        myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(iLayer*he/10) + " \n") #+height/200
                                    myfile.write(" \n")
                                    myfile.write("CELLS " + str(Ncells) + " " + str(int(Ncells+sum(nCell2ThickWalls))) + " \n") #The number of cells corresponds to the number of intercellular spaces
                                    for cid in range(Ncells):
                                        n=int(nCell2ThickWalls[cid]) #Total number of thick wall nodes around the protoplast
                                        Polygon=Cell2ThickWalls[cid][:n]
                                        ranking=list()
                                        ranking.append(int(Polygon[0]))
                                        ranking.append(ThickWalls[int(ranking[0])][5])
                                        ranking.append(ThickWalls[int(ranking[0])][6])
                                        for id1 in range(1,n):
                                            wid1=ThickWalls[int(ranking[id1])][5]
                                            wid2=ThickWalls[int(ranking[id1])][6]
                                            if wid1 not in ranking:
                                                ranking.append(wid1)
                                            if wid2 not in ranking:
                                                ranking.append(wid2)
                                        string=str(n)
                                        for id1 in ranking:
                                            string=string+" "+str(int(id1))
                                        myfile.write(string + " \n")
                                    myfile.write(" \n")
                                    myfile.write("CELL_TYPES " + str(Ncells) + " \n")
                                    for i in range(Ncells):
                                        myfile.write("6 \n") #Triangle-strip cell type
                                    myfile.write(" \n")
                                    myfile.write("POINT_DATA " + str(len(ThickWalls)) + " \n")
                                    myfile.write("SCALARS Hormone_Symplastic_Relative_Concentration_(-) float \n")
                                    myfile.write("LOOKUP_TABLE default \n")
                                    for ThickWallNode in ThickWalls:
                                        cellnumber1=ThickWallNode[2]-NwallsJun
                                        myfile.write(str(float(soln_C[iLayer*Ntot+int(cellnumber1+NwallsJun)])) + " \n")
                                myfile.close()
                                text_file.close()
                            
                                text_file = open(newpath+"Apo_Contagion_bottomb"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+"t00.pvtk", "w")
                                with open(newpath+"Apo_Contagion_bottomb"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+"t00.pvtk", "a") as myfile:
                                    myfile.write("# vtk DataFile Version 4.0 \n")
                                    myfile.write("Apoplastic hormone concentration \n")
                                    myfile.write("ASCII \n")
                                    myfile.write(" \n")
                                    myfile.write("DATASET UNSTRUCTURED_GRID \n")
                                    myfile.write("POINTS "+str(len(ThickWallsX))+" float \n")
                                    for ThickWallNodeX in ThickWallsX:
                                        myfile.write(str(ThickWallNodeX[1]) + " " + str(ThickWallNodeX[2]) + " " + str(iLayer*he/10) + " \n")
                                    myfile.write(" \n")
                                    myfile.write("CELLS " + str(int(NwallsJun+Nwalls-len(list_ghostwalls)*2-len(list_ghostjunctions))) + " " + str(int(2*Nwalls*5-len(list_ghostwalls)*10+sum(nWall2NewWallX[Nwalls:])+NwallsJun-Nwalls+2*len(Wall2NewWallX[Nwalls:])-nGhostJunction2Wall-len(list_ghostjunctions))) + " \n") #The number of cells corresponds to the number of lines in ThickWalls (if no ghost wall & junction)
                                    i=0
                                    for PolygonX in ThickWallPolygonX:
                                        if floor(i/2) not in list_ghostwalls:
                                            myfile.write("4 " + str(int(PolygonX[0])) + " " + str(int(PolygonX[1])) + " " + str(int(PolygonX[2])) + " " + str(int(PolygonX[3])) + " \n")
                                        i+=1
                                    j=Nwalls
                                    for PolygonX in Wall2NewWallX[Nwalls:]: #"junction" polygons
                                        #Would need to order them based on x or y position to make sure display fully covers the surface (but here we try a simpler not so good solution instead)
                                        if j not in list_ghostjunctions:
                                            string=str(int(nWall2NewWallX[j]+2)) #Added +2 so that the first and second nodes could be added again at the end (trying to fill the polygon better)
                                            for id1 in range(int(nWall2NewWallX[j])):
                                                string=string+" "+str(int(PolygonX[id1]))
                                            string=string+" "+str(int(PolygonX[0]))+" "+str(int(PolygonX[1])) #Adding the 1st and 2nd nodes again to the end
                                            myfile.write(string + " \n")
                                        j+=1
                                    myfile.write(" \n")
                                    myfile.write("CELL_TYPES " + str(NwallsJun+Nwalls-len(list_ghostwalls)*2-len(list_ghostjunctions)) + " \n")
                                    i=0
                                    for PolygonX in ThickWallPolygonX:
                                        if floor(i/2) not in list_ghostwalls:
                                            myfile.write("7 \n") #Polygon cell type (wall)
                                        i+=1
                                    j=Nwalls
                                    for PolygonX in Wall2NewWallX[Nwalls:]:
                                        if j not in list_ghostjunctions:
                                            myfile.write("6 \n") #Triangle-strip cell type (wall junction)
                                        j+=1
                                    myfile.write(" \n")
                                    myfile.write("POINT_DATA " + str(len(ThickWallsX)) + " \n")
                                    myfile.write("SCALARS Hormone_Apoplastic_Relative_Concentration_(-) float \n")
                                    myfile.write("LOOKUP_TABLE default \n")
                                    Newsoln_C=zeros((len(ThickWallsX),1))
                                    j=0
                                    for PolygonX in Wall2NewWallX:
                                        for id1 in range(int(nWall2NewWallX[j])):
                                            Newsoln_C[int(PolygonX[id1])]=soln_C[iLayer*Ntot+j]
                                        j+=1
                                    for i in range(len(ThickWallsX)):
                                        myfile.write(str(float(Newsoln_C[i])) + " \n")
                                myfile.close()
                                text_file.close()
                            else:
                                print('Cannot print Apo and Sym_contagion_bottom for iLayer',iLayer,'out of AxialLayers range',AxialLayers)
                            
                    if Transient_tracer:
                        #soln_C_print=zeros((Ntot,AxialLayers,Nprint+1),dtype=float)
                        #for i in range(Ntot):
                        #    soln_C_print[i,:,0]=soln_C[i:AxialLayers*Ntot:Ntot]
                        if PileUp==2:
                            b=9
                        else:
                            b=Barrier
                        if Ini_cond==2:
                            iresume="_resume"
                        else:
                            iresume=""
                        text_file = open(newpath+"Observation_nodes_cc"+str(b)+","+str(iMaturity)+"s"+str(count)+iresume+".txt", "w")
                        with open(newpath+"Observation_nodes_cc"+str(b)+","+str(iMaturity)+"s"+str(count)+iresume+".txt", "a") as myfile:
                            myfile.write("Observation points of tracer concentration"+" \n")
                            myfile.write("\t\tNodes (#) \n")
                            myfile.write("\t\t")
                            if int(Obs_node_range[0].get("id"))==-2: #All cells of the bottom layer
                                for node in Obs_nodes:
                                    myfile.write(str(node-NwallsJun) + "\t\t\t")
                            else:
                                for node in Obs_node_range:
                                    myfile.write(str(node.get("id")) + "\t\t\t")
                            myfile.write("\n")
                            myfile.write("\t\tLayers (#) \n")
                            myfile.write("\t\t")
                            if int(Obs_node_range[0].get("id"))==-2: #All cells of the bottom layer
                                for node in Obs_nodes:
                                    myfile.write(str(0) + "\t\t\t")
                            else:
                                for node in Obs_node_range:
                                    myfile.write(str(node.get("layer")) + "\t\t\t")
                            myfile.write("\n")
                            myfile.write("Time(s)\tConcentration(-) \n")
                            if Ini_cond==1:
                                myfile.write("0.000\t")
                            else:
                                myfile.write(str(np.round(t_resume*3600*24, 1)) + "\t")
                            if len(Obs_nodes)==1:
                                for cid in Obs_nodes:
                                    if cid<len(soln_C):
                                        if soln_C[cid]==0.0:
                                            myfile.write("0.000000\t")
                                        else:
                                            myfile.write(str('%.12E' % Decimal(soln_C[cid])) + "\t")
                                    else:
                                        print('Could not print concentration for observation node',cid,'out of range',len(soln_C))
                                myfile.write("\n")
                            else:
                                for cid in Obs_nodes:
                                    if cid<len(soln_C):
                                        if soln_C[cid]==0.0:
                                            myfile.write("0.000000\t")
                                        else:
                                            if Obs_node_range[0].get("precision")=="high":
                                                myfile.write(str('%.12E' % Decimal(soln_C[cid])) + "\t")
                                            else:
                                                myfile.write(str('%.2E' % Decimal(soln_C[cid])) + "\t")
                                    else:
                                        print('Could not print concentration for observation node',cid,'out of range',len(soln_C))
                                myfile.write("\n")
                        myfile.close()
                        #text_file.close()
                    if PileUp==1:
                        i=0
                        for cid in listxyl:
                            Sym_cc_flows_ini[h][iMaturity][count][N_AxialNeumann_ini.index(cid-NwallsJun)][0]=soln_C[cid]
                            i+=1
                        if Transient_tracer:
                            i=0
                            for cid in listxyl:
                                Pile_cc_time_series[count][i][0]=soln_C[cid] #Concentration at time 0 (steady-state) is recorded at all times (dt_tracer scale). Count==0 at BC, Count==1 in first slice.
                                i+=1
                    
                    if not (Transient_tracer and Ini_cond==2):
                        t01 = time.perf_counter()
                        print(t01-t00, 'seconds process time to solve steady-state tracer flow & write outputs')
                    
                    if Transient_tracer:
                        #Preparing xml saving of tracer transient flow
                        #Stationnary conditions first
                        if Ini_cond==2:
                            it_resume=int(Ini_path[-6:-4])
                        else:
                            it_resume=0
                        newtree = ET.parse(dir + 'cellsetdata/' + path_geom)
                        newroot = newtree.getroot()
                        for metadata in newroot.iter('metadata'):
                            times=ET.SubElement(metadata, 't_tracer')
                            layers=ET.SubElement(metadata, 'AxialLayers')
                        if Ini_cond==1:
                            times.set('value', str(0.0))
                        else:
                            times.set('value', str(t_resume))
                        layers.set('value', str(AxialLayers))
                        if PileUp==2:
                            for cell in newroot[1].iter('cell'):
                                cid=int(cell.get("id"))+NwallsJun
                                slices=ET.SubElement(cell, 'slices')
                                for iLayer in range(AxialLayers):
                                    slice_iL=ET.SubElement(slices, 'slice')
                                    slice_iL.set('z_prox', str(Cumheights[iLayer,0]+heights[iLayer,0]))#position of the proximal end of the root segment (initialization)
                                    slice_iL.set('cc1', str(float(soln_C[iLayer*Ntot+cid])))
                            for wall in newroot[2].iter('wall'):
                                wid=int(wall.get("id"))
                                slices=ET.SubElement(wall, 'slices')
                                for iLayer in range(AxialLayers):
                                    slice_iL=ET.SubElement(slices, 'slice')
                                    slice_iL.set('z_prox', str(Cumheights[iLayer,0]+heights[iLayer,0]))#position of the proximal end of the root segment (initialization)
                                    slice_iL.set('cc1', str(float(soln_C[iLayer*Ntot+wid])))
                            junctions=ET.SubElement(newroot, 'junctions')
                            for jid in range(Nwalls,NwallsJun):
                                junction=ET.SubElement(junctions, 'junction')
                                junction.set('id', str(jid))
                                slices=ET.SubElement(junction, 'slices')
                                for iLayer in range(AxialLayers):
                                    slice_iL=ET.SubElement(slices, 'slice')
                                    slice_iL.set('z_prox', str(Cumheights[iLayer,0]+heights[iLayer,0]))#position of the proximal end of the root segment (initialization)
                                    slice_iL.set('cc1', str(float(soln_C[iLayer*Ntot+jid])))
                        else:
                            for cell in newroot[1].iter('cell'):
                                cid=int(cell.get("id"))+NwallsJun
                                slices=ET.SubElement(cell, 'slices')
                                z_prox=height #position of the proximal end of the root segment (initialization)
                                for iLayer in range(AxialLayers):
                                    slice_iL=ET.SubElement(slices, 'slice')
                                    slice_iL.set('z_prox', str(z_prox))
                                    slice_iL.set('cc1', str(float(soln_C[iLayer*Ntot+cid])))
                                    z_prox+=height
                            for wall in newroot[2].iter('wall'):
                                wid=int(wall.get("id"))
                                slices=ET.SubElement(wall, 'slices')
                                z_prox=height #position of the proximal end of the root segment (initialization)
                                for iLayer in range(AxialLayers):
                                    slice_iL=ET.SubElement(slices, 'slice')
                                    slice_iL.set('z_prox', str(z_prox))
                                    slice_iL.set('cc1', str(float(soln_C[iLayer*Ntot+wid]))) 
                                    z_prox+=height
                            junctions=ET.SubElement(newroot, 'junctions')
                            for jid in range(Nwalls,NwallsJun):
                                junction=ET.SubElement(junctions, 'junction')
                                junction.set('id', str(jid))
                                slices=ET.SubElement(junction, 'slices')
                                z_prox=height #position of the proximal end of the root segment (initialization)
                                for iLayer in range(AxialLayers):
                                    slice_iL=ET.SubElement(slices, 'slice')
                                    slice_iL.set('z_prox', str(z_prox))
                                    slice_iL.set('cc1', str(float(soln_C[iLayer*Ntot+jid])))
                                    z_prox+=height
                        
                        if Ini_cond==2:
                            newtree.write(dir + Project + Ini_path[0:len(Ini_path)-4] + '_resume.xml')
                        else:
                            newtree.write(newpath+'3D_cc_b'+str(b)+','+str(iMaturity)+',t00.xml')
                        
                        #Simulation of transient tracer concentration
                        if Ini_cond==1:
                            t_tracer=0.0
                        else:
                            t_tracer=t_resume
                        soln_C_transi0 = soln_C
                        soln_VC_transi0 = soln_C_transi0*Volume_W2
                        if size(soln_VC_transi0)>size(soln_C_transi0):
                            error('Wrong multiplication type')
                        print('Initial tracer quantity:',float(sum(soln_VC_transi0)))
                        if min(soln_C)<0:
                            t_tracer=t_end
                            i=0
                            for conc in soln_C:
                                if float(conc)<0.0:
                                    print('Negative steady-state concentration:',float(conc),'Node:',i)
                                i+=1
                        
                        #Setting transient Dirichlet BC
                        Stored_matrix_C_transi=[]
                        i=0
                        for nid in N_Dirichlet_transi:
                            if sparseM==1:
                                Stored_matrix_C_transi.append(matrix_C2[nid,:].tocsr())
                                matrix_C2[nid,:]=zeros((1,size(matrix_C2,1))) #Deleting all row nid
                                matrix_C2[nid,nid]=1.0
                            else:
                                Stored_matrix_C_transi.append(matrix_C[nid,:])
                                matrix_C[nid,:]=zeros((1,size(matrix_C,1))) #Deleting all row nid
                                matrix_C[nid,nid]=1.0
                            rhs_C_transi[nid]=cc_Dirichlet_transi[i] #N_Dirichlet_ini.index(nid) #1.0 #Concentration in source protoplasts defined in Geom
                            i+=1
                        
                        if PileUp==0 or PileUp==2: #then rhs_C_transi is the same for all scenarios
                            #Calculation of "asymptotic concentration"
                            print("Solving asymptotic steady-state tracer flow under new transi BC")
                            t00 = time.perf_counter()
                            
                            soln_C_asymptot = empty((Ntot*AxialLayers,1))
                            if sparseM==1:
                                soln_C_asymptot = spsolve(matrix_C2.tocsr(),rhs_C_transi)
                            else:
                                soln_C_asymptot = np.linalg.solve(matrix_C,rhs_C_transi)
                            
                            t01 = time.perf_counter()
                            print(t01-t00, 'seconds process time to solve tracer flow')
                            
                            if min(soln_C_asymptot)<0:
                                t_tracer=t_end
                                i=0
                                for conc in soln_C_asymptot:
                                    if float(conc)<0.0:
                                        print('Negative asymptotic concentration:',float(conc),'Node:',i)
                                    i+=1
                        
                        
                        if switch_nt>0 and switch_dur>0.0:
                            dt_starterC=switch_dur/switch_nt
                            i=0
                            while i < switch_nt:
                                print("Solving transient tracer flow in smooth transition mode")
                                t00 = time.perf_counter()
                                if Num_meth==1: #Explicit scheme
                                    if sparseM==1:
                                        dCdt=np.divide(matrix_C2.tocsr().dot(soln_C_transi0)-(i+1)/switch_nt*rhs_C_transi-(switch_nt-i-1)/switch_nt*rhs_C,Volume_W2)
                                    else:
                                        dCdt=np.divide(np.dot(matrix_C,soln_C_transi0)-(i+1)/switch_nt*rhs_C_transi-(switch_nt-i-1)/switch_nt*rhs_C,Volume_W)
                                    soln_C_transi1=soln_C_transi0+dt_starterC*dCdt
                                elif Num_meth==3: #Implicit scheme 2
                                    if sparseM==1:
                                        A=spdiags(Volume_W2,0,AxialLayers*Ntot,AxialLayers*Ntot).tocsr()-dt_starterC/2*matrix_C2.tocsr()
                                        #A=identity(Ntot*AxialLayers,dtype=dp)*Volume_W2.T-dt_tracer/2*matrix_C2.tocsr()
                                        B=soln_C_transi0*Volume_W2+dt_starterC/2*matrix_C2.tocsr().dot(soln_C_transi0)-((i+1)/switch_nt*rhs_C_transi+(switch_nt-i-1)/switch_nt*rhs_C)*dt_starterC
                                        soln_C_transi1 = spsolve(A,B) #Trapezoidal rule
                                    else:
                                        A=np.identity(Ntot,dtype=dp)*Volume_W-dt_starterC/2*matrix_C
                                        B=soln_C_transi0*Volume_W+dt_starterC/2*np.dot(matrix_C,soln_C_transi0)-((i+1)/switch_nt*rhs_C_transi+(switch_nt-i-1)/switch_nt*rhs_C)*dt_starterC
                                        soln_C_transi1 = np.linalg.solve(A,B) #Trapezoidal rule
                                if min(soln_C_transi1)<-1E-1:
                                    error('Negative concentration, decrease time step')
                                elif min(soln_C_transi1)<0:
                                    for nid in range(0,Ntot*AxialLayers):
                                        if soln_C_transi1[nid]<0:
                                            print('Slightly negative concentration',nid)#soln_C_transi1[nid]=0.0
                                soln_C_transi0=soln_C_transi1
                                t01 = time.perf_counter()
                                i+=1
                                print(t01-t00, 'seconds process time, progress '+str(100*i/switch_nt)+'%')
                            t_tracer+=switch_dur
                        
                        #Time loop
                        it=1 #Print time index
                        idt=1 #Small time step index
                        for Print_time in Print_times:
                            t2=float(Print_time.get('value'))
                            while t_tracer+dt_tracer < t2:
                                if PileUp==1: #Then rhs_C_transi is redefined for each scenario(=slice) from input water concentration in the slice below
                                    rhs_C_transi = np.zeros((AxialLayers*Ntot),dtype=dp)
                                    
                                    #Dirichlet solute BC
                                    for i in N_Dirichlet_transi: 
                                        rhs_C_transi[i] = cc_Dirichlet_transi[N_Dirichlet_transi.index(i)]
                                    
                                    #Neumann axial solute BC
                                    if Barrier>0:
                                        i=0
                                        for cid in listxyl:
                                            Q_dist=float(Q_xyl[0][i]) #Positive upwards
                                            Q_prox=float(Q_xyl[1][i])
                                            if Q_prox<0: #Flow from shoot towards proximal xylem vessel
                                                print('No proximal solute flux BC when using the PileUp mode')
                                            if Q_dist>0:#Flow from root tip towards distal xylem vessel
                                                if cid-NwallsJun in N_AxialNeumann_transi: #if cid not in Nodes_C_Dirichlet_ini:
                                                    rhs_C_transi[cid] -= Q_dist*Pile_cc_time_series[count-1][N_AxialNeumann_transi.index(cid-NwallsJun)][idt-1]
                                            i+=1
                                    
                                    #Phloem BC for solute
                                    i=0
                                    for cid in listsieve:
                                        if Q_sieve[1][i]<0: #Flow from shoot towards proximal phloem sieve tube
                                            if cid-NwallsJun in N_AxialNeumann_transi: #if cid not in Nodes_C_Dirichlet_ini:
                                                print('No proximal solute flux BC when using the PileUp==1 mode')
                                        if Q_sieve[0][i]>0: #Flow from shoot towards proximal phloem sieve tube
                                            if cid-NwallsJun in N_AxialNeumann_transi: #if cid not in Nodes_C_Dirichlet_ini:
                                                rhs_C_transi[cid] -= Q_sieve[0][i]*Pile_cc_time_series[count-1][N_AxialNeumann_transi.index(cid-NwallsJun)][idt-1]
                                        i+=1
                                    
                                    #Calculation of "asymptotic concentration"
                                    soln_C_asymptot = empty((Ntot*AxialLayers,1))
                                    if sparseM==1:
                                        soln_C_asymptot = spsolve(matrix_C2.tocsr(),rhs_C_transi)
                                    else:
                                        soln_C_asymptot = np.linalg.solve(matrix_C,rhs_C_transi)
                                
                                #Calculating distance to "asymptotic concentrations" under new BC
                                dC_asymptot = soln_C_asymptot-soln_C_transi0
                                if any(soln_C_transi0==0.0):
                                    print('Relative distance from asymptot (cc diff) between:',min(dC_asymptot),'and',max(dC_asymptot))
                                else:
                                    temp=np.divide(abs(dC_asymptot),soln_C_transi0)
                                    if max(temp)<0.001: #relative differences between current state and asymptot lower than 0.1%
                                        t_tracer=t_end-dt_tracer
                                    else:
                                        print('Relative distance from asymptot (%) between:',100*min(abs(np.divide(dC_asymptot,soln_C_asymptot-soln_C))),'and',100*max(abs(np.divide(dC_asymptot,soln_C_asymptot-soln_C))))
                                
                                #Solving transient tracer concentration in slice "count" at time idt
                                if dt_adapt:
                                    dcc_OK=False
                                    while not dcc_OK:
                                        print("Solving transient tracer flow adaptive dt")
                                        t00 = time.perf_counter()
                                        if Num_meth==1: #Explicit scheme
                                            if sparseM==1:
                                                dCdt=np.divide(matrix_C2.tocsr().dot(soln_C_transi0)-rhs_C_transi,Volume_W2)
                                            else:
                                                dCdt=np.divide(np.dot(matrix_C,soln_C_transi0)-rhs_C_transi,Volume_W)
                                            soln_C_transi1=soln_C_transi0+dt_tracer*dCdt
                                        elif Num_meth==2: #Implicit scheme 1
                                            error('This type of implicit scheme is not accurate, please use numerical method 3')
                                        elif Num_meth==3: #Implicit scheme 2
                                            if sparseM==1:
                                                A=spdiags(Volume_W2,0,AxialLayers*Ntot,AxialLayers*Ntot).tocsr()-dt_tracer/2*matrix_C2.tocsr()
                                                #A=identity(Ntot*AxialLayers,dtype=dp)*Volume_W2.T-dt_tracer/2*matrix_C2.tocsr()
                                                B=soln_C_transi0*Volume_W2+dt_tracer/2*matrix_C2.tocsr().dot(soln_C_transi0)-rhs_C_transi*dt_tracer
                                                soln_C_transi1 = spsolve(A,B) #Trapezoidal rule
                                            else:
                                                A=np.identity(Ntot,dtype=dp)*Volume_W-dt_tracer/2*matrix_C
                                                B=soln_C_transi0*Volume_W+dt_tracer/2*np.dot(matrix_C,soln_C_transi0)-rhs_C_transi*dt_tracer
                                                soln_C_transi1 = np.linalg.solve(A,B) #Trapezoidal rule
                                        t01 = time.perf_counter()
                                        print(t01-t00, 'seconds process time to solve transient tracer flow')

                                        #dcc_max=max(-nanmin(nanmax((soln_C_transi1-soln_C_transi0)/soln_C_transi0)),nanmax((soln_C_transi1-soln_C_transi0)/soln_C_transi0))
                                        dcc_max=max(-min(soln_C_transi1-soln_C_transi0),max(soln_C_transi1-soln_C_transi0))
                                        if dcc_max>dcc_high:
                                            if dt_tracer/dt_factor<dt_min:
                                                if not dt_tracer==dt_min:
                                                    dt_tracer=dt_min
                                                    print('Concentration change ',dcc_max,' > ',dcc_high,'-> Time step decrease: dt_tracer = ',dt_tracer)
                                                else:
                                                    print('Concentration change ',dcc_max,' > ',dcc_high,' but dt_tracer == dt_min ',dt_tracer)
                                                    dcc_OK=True
                                            else:
                                                dt_tracer/=dt_factor
                                                print('Concentration change ',dcc_max,' > ',dcc_high,'-> Time step decrease: dt_tracer = ',dt_tracer)
                                        else:
                                            dcc_OK=True
                                            if dcc_max<dcc_low:
                                                if dt_tracer*dt_factor>dt_max:
                                                    if not dt_tracer==dt_max:
                                                        dt_tracer=dt_max
                                                        print('Concentration change ',dcc_max,' < ',dcc_low,'-> Time step increase: dt_tracer = ',dt_tracer)
                                                    else:
                                                        print('Concentration change ',dcc_max,' < ',dcc_low,' but dt_tracer == dt_max ',dt_tracer)
                                                else:
                                                    dt_tracer*=dt_factor
                                                    print('Concentration change ',dcc_max,' < ',dcc_low,'-> Time step increase: dt_tracer = ',dt_tracer)
                                else:    
                                    print("Solving transient tracer flow cst dt")
                                    t00 = time.perf_counter()
                                
                                    if Num_meth==1: #Explicit scheme
                                        if sparseM==1:
                                            dCdt=np.divide(matrix_C2.tocsr().dot(soln_C_transi0)-rhs_C_transi,Volume_W2)
                                        else:
                                            dCdt=np.divide(np.dot(matrix_C,soln_C_transi0)-rhs_C_transi,Volume_W)
                                        soln_C_transi1=soln_C_transi0+dt_tracer*dCdt
                                    elif Num_meth==2: #Implicit scheme 1
                                        error('This type of implicit scheme is not accurate, please use numerical method 3')
                                        #if sparseM==1:
                                        #    matrix_Cdt2V=dt_tracer/2.0*(matrix_C.T/Volume_W).T
                                        #    A=np.identity(Ntot,dtype=dp)-matrix_Cdt2V
                                        #    B=soln_C_transi0+np.dot(matrix_Cdt2V,soln_C_transi0)-np.divide(rhs_C_transi*dt_tracer,Volume_W)
                                        #    soln_C_transi1 = np.linalg.solve(A,B) #Trapezoidal rule
                                        #else:
                                        #    matrix_Cdt2V=dt_tracer/2.0*(matrix_C.T/Volume_W).T
                                        #    A=np.identity(Ntot,dtype=dp)-matrix_Cdt2V
                                        #    B=soln_C_transi0+np.dot(matrix_Cdt2V,soln_C_transi0)-np.divide(rhs_C_transi*dt_tracer,Volume_W)
                                        #    soln_C_transi1 = np.linalg.solve(A,B) #Trapezoidal rule
                                        ##Avoiding division by volume in increment
                                        ##soln_VC_transi1 = soln_VC_transi0 + dt_tracer*(np.dot(matrix_C,soln_C_transi0)-rhs_C_transi)
                                    elif Num_meth==3: #Implicit scheme 2
                                        if sparseM==1:
                                            A=spdiags(Volume_W2,0,AxialLayers*Ntot,AxialLayers*Ntot).tocsr()-dt_tracer/2*matrix_C2.tocsr()
                                            #A=identity(Ntot*AxialLayers,dtype=dp)*Volume_W2.T-dt_tracer/2*matrix_C2.tocsr()
                                            B=soln_C_transi0*Volume_W2+dt_tracer/2*matrix_C2.tocsr().dot(soln_C_transi0)-rhs_C_transi*dt_tracer
                                            soln_C_transi1 = spsolve(A,B) #Trapezoidal rule
                                        else:
                                            A=np.identity(Ntot,dtype=dp)*Volume_W-dt_tracer/2*matrix_C
                                            B=soln_C_transi0*Volume_W+dt_tracer/2*np.dot(matrix_C,soln_C_transi0)-rhs_C_transi*dt_tracer
                                            soln_C_transi1 = np.linalg.solve(A,B) #Trapezoidal rule
                                    t01 = time.perf_counter()
                                    print(t01-t00, 'seconds process time to solve transient tracer flow')
                                
                                if sparseM==1:    
                                    soln_VC_transi1=soln_C_transi1*Volume_W2
                                    if size(soln_VC_transi1)>size(soln_C_transi1):
                                        error('Wrong multiplication type')
                                else:
                                    soln_VC_transi1=soln_C_transi1*Volume_W
                                    if size(soln_VC_transi1)>size(soln_C_transi1):
                                        error('Wrong multiplication type')
                                C_tot0=float(sum(soln_VC_transi0))
                                delta_C_transi=float(sum(soln_VC_transi1))-C_tot0
                                if size(N_AxialNeumann_transi)>0:
                                    Flux_div_C=-float(sum(rhs_C_transi))-Degrad1*float(np.dot(Volume_W2,soln_C))
                                    if not C_tot0==0.0:
                                        print('Time (d):',float(t_tracer+dt_tracer),'   Tracer quantity and storage change:',C_tot0,delta_C_transi,'   Mass balance error (%)',100*(delta_C_transi-Flux_div_C*dt_tracer)/C_tot0)#,'   Xylem input & output:',float(-dt_tracer*sum(rhs_C_transi)),float(-dt_tracer*sum(Flow_xyl[0][1:Nxyl+1,[count]]*soln_C_transi0[Ntot:Ntot+Nxyl])),'   Soil surface output:',sum(soln_C_transi0[hstack((Borderwall,Borderjunction))]*Q_soil_neg)) #concatenate
                                    else:
                                        print('Time (d):',float(t_tracer+dt_tracer),'   Tracer quantity and storage change:',C_tot0,delta_C_transi)
                                else:
                                    if len(N_Dirichlet_transi)<1000:
                                        Flux_div_C=-Degrad1*float(np.dot(Volume_W2,soln_C))
                                        for i_node in N_Dirichlet_transi:
                                            Flux_div_C+=float((sum(matrix_C2[:,i_node])-1.0)*rhs_C_transi[i_node])
                                        if not C_tot0==0.0:
                                            print('Time (d):',float(t_tracer+dt_tracer),'   Tracer quantity and storage change:',C_tot0,delta_C_transi,'   Mass balance error (%)',100*(delta_C_transi-Flux_div_C*dt_tracer)/C_tot0)#,'   Xylem input & output:',float(-dt_tracer*sum(rhs_C_transi)),float(-dt_tracer*sum(Flow_xyl[0][1:Nxyl+1,[count]]*soln_C_transi0[Ntot:Ntot+Nxyl])),'   Soil surface output:',sum(soln_C_transi0[hstack((Borderwall,Borderjunction))]*Q_soil_neg)) #concatenate
                                        else:
                                            print('Time (d):',float(t_tracer+dt_tracer),'   Tracer quantity and storage change:',C_tot0,delta_C_transi)
                                    else:
                                        print('Time (d):',float(t_tracer+dt_tracer),'   Tracer quantity and storage change:',C_tot0,delta_C_transi)
                                if min(soln_C_transi1)<-0.1:
                                    t_tracer=t_end
                                    i=0
                                    for conc in soln_C_transi1:
                                        if float(conc)<0.0: # and i<2*Ntot
                                            print('Negative transient concentration:',float(conc),'Node:',i)
                                        i+=1
                                
                                #Writing observation point output
                                if idt%Freq_obs==0.0: #Remainder of division of time step number by frequency of recorded observations
                                    with open(newpath+"Observation_nodes_cc"+str(b)+","+str(iMaturity)+"s"+str(count)+iresume+".txt", "a") as myfile:
                                        myfile.write(str(np.round((t_tracer+dt_tracer)*3600*24, 1)) + "\t")
                                        if len(Obs_nodes)==1:
                                            for cid in Obs_nodes:
                                                if cid<len(soln_C_transi1):
                                                    if soln_C_transi1[cid]==0.0:
                                                        myfile.write("0.000000\t")
                                                    else:
                                                        if Obs_node_range[0].get("precision")=="high":
                                                            myfile.write(str('%.12E' % Decimal(soln_C_transi1[cid])) + "\t")
                                                        else:
                                                            myfile.write(str('%.2E' % Decimal(soln_C_transi1[cid])) + "\t")
                                                else:
                                                    print('Could not print concentration for observation node',cid, 'out of range',len(soln_C_transi1))
                                            myfile.write("\n")
                                        else:
                                            for cid in Obs_nodes:
                                                if cid<len(soln_C_transi1):
                                                    if soln_C_transi1[cid]==0.0:
                                                        myfile.write("0.000000\t")
                                                    else:
                                                        myfile.write(str('%.2E' % Decimal(soln_C_transi1[cid])) + "\t")
                                                else:
                                                    print('Could not print concentration for observation node',cid, 'out of range',len(soln_C_transi1))
                                            myfile.write("\n")
                                    myfile.close()
                                #text_file.close()
                                
                                #Saving xylem tracer concentration across dt
                                if PileUp==1:
                                    i=0
                                    for cid in listxyl:
                                        Pile_cc_time_series[count][i][idt] = soln_C_transi0[cid][0]
                                        i+=1
                                #Preparing inputs for next time step
                                soln_C_transi0=soln_C_transi1
                                soln_VC_transi0=soln_VC_transi1
                                #temp1=100.0*np.divide(soln_C_transi0-soln_C,soln_C)
                                #print(1,temp1.tolist().index(max(temp1)),max(temp1),temp1.tolist().index(min(temp1)),min(temp1))
                                t_tracer+=dt_tracer
                                idt+=1
                            #end of the loop "while t_tracer+dt_tracer < t2:"
                            
                            if PileUp==1: #Then rhs_C_transi is redefined for each scenario(=slice) from input water concentration in the slice below
                                rhs_C_transi = np.zeros((AxialLayers*Ntot))
                                
                                #Dirichlet solute BC
                                for i in N_Dirichlet_transi:
                                    rhs_C_transi[i] = cc_Dirichlet_transi[N_Dirichlet_transi.index(i)]
                                
                                #Neumann solute BC
                                if Barrier>0:
                                    i=0
                                    for cid in listxyl:
                                        Q_dist=float(Q_xyl[0][i]) #Positive upwards
                                        Q_prox=float(Q_xyl[1][i])
                                        if Q_prox<0: #Flow from shoot towards proximal xylem vessel
                                            print('No proximal solute flux BC when using the PileUp mode')
                                        if Q_dist>0:#Flow from root tip towards distal xylem vessel
                                            if cid-NwallsJun in N_AxialNeumann_transi: #if cid not in Nodes_C_Dirichlet_ini:
                                                rhs_C_transi[cid] -= Q_dist*Pile_cc_time_series[count-1][N_AxialNeumann_transi.index(cid-NwallsJun)][idt-1]
                                        i+=1
                                
                                #Phloem BC for solute
                                i=0
                                for cid in listsieve:
                                    if Q_sieve[1][i]<0: #Flow from shoot towards proximal phloem sieve tube
                                        if cid-NwallsJun in N_AxialNeumann_transi: #if cid not in Nodes_C_Dirichlet_ini:
                                            print('No proximal solute flux BC when using the PileUp==1 mode')
                                    if Q_sieve[0][i]>0: #Flow from tip towards distal phloem sieve tube
                                        if cid-NwallsJun in N_AxialNeumann_transi: #if cid not in Nodes_C_Dirichlet_ini:
                                            rhs_C_transi[cid] -= Q_sieve[0][i]*Pile_cc_time_series[count-1][N_AxialNeumann_transi.index(cid-NwallsJun)][idt-1]
                                    i+=1
                                
                            # #Check divergence of rhs_C over time
                            # temp=zeros((Ntot+Nxyl,1))
                            # for i in range(Ntot+Nxyl):
                            #     if not rhs_C[i]==0:
                            #         temp[i]=100.0*np.divide(rhs_C_transi[i]-rhs_C[i],rhs_C[i])
                            # print(2,max(temp))
                            #Solving transient tracer concentration in slice "count" in last bit before print time  
                            dt=t2-t_tracer
                            if Num_meth==1: #Explicit scheme
                                if sparseM==1:
                                    dCdt=np.divide(matrix_C2.tocsr().dot(soln_C_transi0)-rhs_C_transi,Volume_W2)
                                else:
                                    dCdt=np.divide(np.dot(matrix_C,soln_C_transi0)-rhs_C_transi,Volume_W)
                                soln_C_transi1=soln_C_transi0+dt*dCdt
                                #Avoiding division by volume in increment
                                #soln_VC_transi1 = soln_VC_transi0 + dt*(np.dot(matrix_C,soln_C_transi0)-rhs_C_transi)
                            elif Num_meth==2: #Implicit scheme 1
                                error('This type of implicit scheme is not accurate, please use numerical method 3')
                            #    matrix_Cdt2V=dt/2.0*(matrix_C.T/Volume_W).T
                            #    A=np.identity(Ntot,dtype=dp)-matrix_Cdt2V
                            #    B=soln_C_transi0+np.dot(matrix_Cdt2V,soln_C_transi0)-np.divide(rhs_C_transi*dt,Volume_W)
                            #    soln_C_transi1 = np.linalg.solve(A,B) #Trapezoidal rule
                            #    #soln_C_transi1 = np.linalg.solve(np.identity(Ntot+Nxyl,dtype=dp)-matrix_Cdt2V,soln_C_transi0+np.dot(matrix_Cdt2V,soln_C_transi0)-np.divide(rhs_C_transi,Volume_W)) #Trapezoidal rule
                            elif Num_meth==3: #Implicit scheme 2
                                if sparseM==1:
                                    A=spdiags(Volume_W2,0,AxialLayers*Ntot,AxialLayers*Ntot).tocsr()-dt/2*matrix_C2.tocsr()
                                    #A=identity(Ntot*AxialLayers,dtype=dp)*Volume_W2.T-dt_tracer/2*matrix_C2.tocsr()
                                    B=soln_C_transi0*Volume_W2+dt/2*matrix_C2.tocsr().dot(soln_C_transi0)-rhs_C_transi*dt
                                    soln_C_transi1 = spsolve(A,B) #Trapezoidal rule
                                else:
                                    A=np.identity(Ntot,dtype=dp)*Volume_W-dt/2*matrix_C
                                    B=soln_C_transi0*Volume_W+dt/2*np.dot(matrix_C,soln_C_transi0)-rhs_C_transi*dt
                                    soln_C_transi1 = np.linalg.solve(A,B) #Trapezoidal rule
                            if sparseM==1:
                                soln_VC_transi1=soln_C_transi1*Volume_W2
                                if size(soln_VC_transi1)>size(soln_C_transi1):
                                        error('Wrong multiplication type')
                            else:
                                soln_VC_transi1=soln_C_transi1*Volume_W
                                if size(soln_VC_transi1)>size(soln_C_transi1):
                                        error('Wrong multiplication type')
                            C_tot0=float(sum(soln_VC_transi0))
                            delta_C_transi=float(sum(soln_VC_transi1))-C_tot0
                            if size(N_AxialNeumann_transi)>0:
                                Flux_div_C=float(-sum(rhs_C_transi))-Degrad1*float(np.dot(Volume_W2,soln_C))
                                if not C_tot0==0.0:
                                    print('Time2(d):',float(t_tracer+dt),'   Tracer quantity and storage change:',C_tot0,delta_C_transi,'   Mass balance error (%)',100*(delta_C_transi-Flux_div_C*dt)/C_tot0)#,'Xylem input & output:',float(-dt*sum(rhs_C_transi)),float(-dt*sum(Flow_xyl[0][1:Nxyl+1,[count]]*soln_C_transi0[Ntot:Ntot+Nxyl])),'Soil surface output:',sum(soln_C_transi0[hstack((Borderwall,Borderjunction))]*Q_soil_neg))
                                else:
                                    print('Time2(d):',float(t_tracer+dt),'   Tracer quantity and storage change:',C_tot0,delta_C_transi)
                            else:
                                if len(N_Dirichlet_transi)<1000:
                                    Flux_div_C=-Degrad1*float(np.dot(Volume_W2,soln_C))
                                    for i_node in N_Dirichlet_transi:
                                        Flux_div_C+=float((sum(matrix_C2[:,i_node])-1.0)*rhs_C_transi[i_node])
                                    if not C_tot0==0.0:
                                        print('Time2(d):',float(t_tracer+dt),'   Tracer quantity and storage change:',C_tot0,delta_C_transi,'   Mass balance error (%)',100*(delta_C_transi-Flux_div_C*dt)/C_tot0)#,'Xylem input & output:',float(-dt*sum(rhs_C_transi)),float(-dt*sum(Flow_xyl[0][1:Nxyl+1,[count]]*soln_C_transi0[Ntot:Ntot+Nxyl])),'Soil surface output:',sum(soln_C_transi0[hstack((Borderwall,Borderjunction))]*Q_soil_neg))
                                    else:
                                        print('Time2(d):',float(t_tracer+dt),'   Tracer quantity and storage change:',C_tot0,delta_C_transi)
                                else:
                                    print('Time2(d):',float(t_tracer+dt),'   Tracer quantity and storage change:',C_tot0,delta_C_transi)
                            
                            #Writing observation point output
                            if idt%Freq_obs==0.0: #Remainder of division of time step number by frequency of recorded observations
                                with open(newpath+"Observation_nodes_cc"+str(b)+","+str(iMaturity)+"s"+str(count)+iresume+".txt", "a") as myfile:
                                    myfile.write(str(np.round(t2*3600*24, 1)) + "\t")
                                    if len(Obs_nodes)==1:
                                        for cid in Obs_nodes:
                                            if soln_C_transi1[cid]==0.0:
                                                myfile.write("0.000000\t")
                                            else:
                                                myfile.write(str('%.12E' % Decimal(soln_C_transi1[cid])) + "\t")
                                        myfile.write("\n")
                                    else:
                                        for cid in Obs_nodes:
                                            if soln_C_transi1[cid]==0.0:
                                                myfile.write("0.000000\t")
                                            else:
                                                if Obs_node_range[0].get("precision")=="high":
                                                    myfile.write(str('%.12E' % Decimal(soln_C_transi1[cid])) + "\t")
                                                else:
                                                    myfile.write(str('%.2E' % Decimal(soln_C_transi1[cid])) + "\t")
                                        myfile.write("\n")
                                myfile.close()
                            
                            #Saving xylem tracer concentration across dt
                            if PileUp==1:
                                i=0
                                for cid in listxyl:
                                    Pile_cc_time_series[count][i][idt] = soln_C_transi0[cid]
                                    Sym_cc_flows_transi[h][iMaturity][count][N_AxialNeumann_transi.index(cid-NwallsJun)][0][it]=soln_C_transi0[cid]
                                    i+=1
                            #Preparing inputs for next time step
                            soln_C_transi0=soln_C_transi1
                            soln_VC_transi0=soln_VC_transi1
                            #soln_C_transi0=np.divide(soln_VC_transi1,Volume_W)
                            
                            #temp1=100.0*np.divide(soln_C_transi0-soln_C,soln_C)
                            #print(2,temp1.tolist().index(max(temp1)),max(temp1),temp1.tolist().index(min(temp1)),min(temp1)) #find(
                            
                            t_tracer=t2
                            idt+1
                            
                            if Paraview==1 or ParTrack==1:
                                #Save data for .vtk and print xml output
                                for Print_layer in Printrange:
                                    iLayer=int(Print_layer.get("value"))
                                    if iLayer<AxialLayers:
                                        if PileUp==2:
                                            b=int(BarriersL[iLayer])
                                        else:
                                            b=Barrier
                                        if it+it_resume>9:
                                            fullpath=newpath+"Sym_Contagion_bottomb"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+"t"+str(it+it_resume)+".pvtk"
                                        else:
                                            fullpath=newpath+"Sym_Contagion_bottomb"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+"t0"+str(it+it_resume)+".pvtk"
                                        text_file = open(fullpath, "w")
                                        with open(fullpath, "a") as myfile:
                                            myfile.write("# vtk DataFile Version 4.0 \n")
                                            myfile.write("Symplastic hormone concentration \n")
                                            myfile.write("ASCII \n")
                                            myfile.write(" \n")
                                            myfile.write("DATASET UNSTRUCTURED_GRID \n")
                                            myfile.write("POINTS "+str(len(ThickWalls))+" float \n")
                                            for ThickWallNode in ThickWalls:
                                                myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(Cumheights[iLayer,0]/10) + " \n") #+height/200
                                            myfile.write(" \n")
                                            myfile.write("CELLS " + str(Ncells) + " " + str(int(Ncells+sum(nCell2ThickWalls))) + " \n") #The number of cells corresponds to the number of intercellular spaces
                                            for cid in range(Ncells):
                                                n=int(nCell2ThickWalls[cid]) #Total number of thick wall nodes around the protoplast
                                                Polygon=Cell2ThickWalls[cid][:n]
                                                ranking=list()
                                                ranking.append(int(Polygon[0]))
                                                ranking.append(ThickWalls[int(ranking[0])][5])
                                                ranking.append(ThickWalls[int(ranking[0])][6])
                                                for id1 in range(1,n):
                                                    wid1=ThickWalls[int(ranking[id1])][5]
                                                    wid2=ThickWalls[int(ranking[id1])][6]
                                                    if wid1 not in ranking:
                                                        ranking.append(wid1)
                                                    if wid2 not in ranking:
                                                        ranking.append(wid2)
                                                string=str(n)
                                                for id1 in ranking:
                                                    string=string+" "+str(int(id1))
                                                myfile.write(string + " \n")
                                            myfile.write(" \n")
                                            myfile.write("CELL_TYPES " + str(Ncells) + " \n")
                                            for i in range(Ncells):
                                                myfile.write("6 \n") #Triangle-strip cell type
                                            myfile.write(" \n")
                                            myfile.write("POINT_DATA " + str(len(ThickWalls)) + " \n")
                                            myfile.write("SCALARS Hormone_Symplastic_Relative_Concentration_(-) float \n")
                                            myfile.write("LOOKUP_TABLE default \n")
                                            for ThickWallNode in ThickWalls:
                                                cellnumber1=ThickWallNode[2]-NwallsJun
                                                myfile.write(str(float(soln_C_transi0[int(cellnumber1+NwallsJun)+iLayer*Ntot])) + " \n")
                                        myfile.close()
                                        text_file.close()
                                    
                                        if it+it_resume>9:
                                            fullpath=newpath+"Apo_Contagion_bottomb"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+"t"+str(it+it_resume)+".pvtk"
                                        else:
                                            fullpath=newpath+"Apo_Contagion_bottomb"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+"t0"+str(it+it_resume)+".pvtk"
                                        text_file = open(fullpath, "w")
                                        with open(fullpath, "a") as myfile:
                                            myfile.write("# vtk DataFile Version 4.0 \n")
                                            myfile.write("Apoplastic hormone concentration \n")
                                            myfile.write("ASCII \n")
                                            myfile.write(" \n")
                                            myfile.write("DATASET UNSTRUCTURED_GRID \n")
                                            myfile.write("POINTS "+str(len(ThickWallsX))+" float \n")
                                            for ThickWallNodeX in ThickWallsX:
                                                myfile.write(str(ThickWallNodeX[1]) + " " + str(ThickWallNodeX[2]) + " " + str(Cumheights[iLayer,0]/10) + " \n")
                                            myfile.write(" \n")
                                            myfile.write("CELLS " + str(int(NwallsJun+Nwalls-len(list_ghostwalls)*2-len(list_ghostjunctions))) + " " + str(int(2*Nwalls*5-len(list_ghostwalls)*10+sum(nWall2NewWallX[Nwalls:])+NwallsJun-Nwalls+2*len(Wall2NewWallX[Nwalls:])-nGhostJunction2Wall-len(list_ghostjunctions))) + " \n") #The number of cells corresponds to the number of lines in ThickWalls (if no ghost wall & junction)
                                            i=0
                                            for PolygonX in ThickWallPolygonX:
                                                if floor(i/2) not in list_ghostwalls:
                                                    myfile.write("4 " + str(int(PolygonX[0])) + " " + str(int(PolygonX[1])) + " " + str(int(PolygonX[2])) + " " + str(int(PolygonX[3])) + " \n")
                                                i+=1
                                            j=Nwalls
                                            for PolygonX in Wall2NewWallX[Nwalls:]: #"junction" polygons
                                                #Would need to order them based on x or y position to make sure display fully covers the surface (but here we try a simpler not so good solution instead)
                                                if j not in list_ghostjunctions:
                                                    string=str(int(nWall2NewWallX[j]+2)) #Added +2 so that the first and second nodes could be added again at the end (trying to fill the polygon better)
                                                    for id1 in range(int(nWall2NewWallX[j])):
                                                        string=string+" "+str(int(PolygonX[id1]))
                                                    string=string+" "+str(int(PolygonX[0]))+" "+str(int(PolygonX[1])) #Adding the 1st and 2nd nodes again to the end
                                                    myfile.write(string + " \n")
                                                j+=1
                                            myfile.write(" \n")
                                            myfile.write("CELL_TYPES " + str(NwallsJun+Nwalls-len(list_ghostwalls)*2-len(list_ghostjunctions)) + " \n")
                                            i=0
                                            for PolygonX in ThickWallPolygonX:
                                                if floor(i/2) not in list_ghostwalls:
                                                    myfile.write("7 \n") #Polygon cell type (wall)
                                                i+=1
                                            j=Nwalls
                                            for PolygonX in Wall2NewWallX[Nwalls:]:
                                                if j not in list_ghostjunctions:
                                                    myfile.write("6 \n") #Triangle-strip cell type (wall junction)
                                                j+=1
                                            myfile.write(" \n")
                                            myfile.write("POINT_DATA " + str(len(ThickWallsX)) + " \n")
                                            myfile.write("SCALARS Hormone_Apoplastic_Relative_Concentration_(-) float \n")
                                            myfile.write("LOOKUP_TABLE default \n")
                                            Newsoln_C=zeros((len(ThickWallsX),1))
                                            j=0
                                            for PolygonX in Wall2NewWallX:
                                                for id1 in range(int(nWall2NewWallX[j])):
                                                    Newsoln_C[int(PolygonX[id1])]=soln_C_transi0[j+iLayer*Ntot]
                                                j+=1
                                            for i in range(len(ThickWallsX)):
                                                myfile.write(str(float(Newsoln_C[i])) + " \n")
                                        myfile.close()
                                        text_file.close()
                                    
                            #Preparing xml saving of tracer transient flow
                            #Non-stationnary conditions
                            if PileUp==2:
                                b=9
                            else:
                                b=Barrier
                            newtree = ET.parse(dir + 'cellsetdata/' + path_geom)
                            newroot = newtree.getroot()
                            for metadata in newroot.iter('metadata'):
                                times=ET.SubElement(metadata, 't_tracer')
                                layers=ET.SubElement(metadata, 'AxialLayers')
                            times.set('value', str(t_tracer))
                            layers.set('value', str(AxialLayers))
                            if PileUp==2:
                                for cell in newroot[1].iter('cell'):
                                    cid=int(cell.get("id"))+NwallsJun
                                    slices=ET.SubElement(cell, 'slices')
                                    for iLayer in range(AxialLayers):
                                        slice_iL=ET.SubElement(slices, 'slice')
                                        slice_iL.set('z_prox', str(Cumheights[iLayer,0]+heights[iLayer,0]))#position of the proximal end of the root segment (initialization)
                                        slice_iL.set('cc1', str(float(soln_C_transi0[iLayer*Ntot+cid])))
                                for wall in newroot[2].iter('wall'):
                                    wid=int(wall.get("id"))
                                    slices=ET.SubElement(wall, 'slices')
                                    for iLayer in range(AxialLayers):
                                        slice_iL=ET.SubElement(slices, 'slice')
                                        slice_iL.set('z_prox', str(Cumheights[iLayer,0]+heights[iLayer,0]))#position of the proximal end of the root segment (initialization)
                                        slice_iL.set('cc1', str(float(soln_C_transi0[iLayer*Ntot+wid])))
                                junctions=ET.SubElement(newroot, 'junctions')
                                for jid in range(Nwalls,NwallsJun):
                                    junction=ET.SubElement(junctions, 'junction')
                                    junction.set('id', str(jid))
                                    slices=ET.SubElement(junction, 'slices')
                                    for iLayer in range(AxialLayers):
                                        slice_iL=ET.SubElement(slices, 'slice')
                                        slice_iL.set('z_prox', str(Cumheights[iLayer,0]+heights[iLayer,0]))#position of the proximal end of the root segment (initialization)
                                        slice_iL.set('cc1', str(float(soln_C_transi0[iLayer*Ntot+jid])))
                            else:
                                for cell in newroot[1].iter('cell'):
                                    cid=int(cell.get("id"))+NwallsJun
                                    slices=ET.SubElement(cell, 'slices')
                                    z_prox=height #position of the proximal end of the root segment (initialization)
                                    for iLayer in range(AxialLayers):
                                        slice_iL=ET.SubElement(slices, 'slice')
                                        slice_iL.set('z_prox', str(z_prox))
                                        slice_iL.set('cc1', str(float(soln_C_transi0[iLayer*Ntot+cid])))
                                        z_prox+=height
                                for wall in newroot[2].iter('wall'):
                                    wid=int(wall.get("id"))
                                    slices=ET.SubElement(wall, 'slices')
                                    z_prox=height #position of the proximal end of the root segment (initialization)
                                    for iLayer in range(AxialLayers):
                                        slice_iL=ET.SubElement(slices, 'slice')
                                        slice_iL.set('z_prox', str(z_prox))
                                        slice_iL.set('cc1', str(float(soln_C_transi0[iLayer*Ntot+wid])))
                                        z_prox+=height
                                junctions=ET.SubElement(newroot, 'junctions')
                                for jid in range(Nwalls,NwallsJun):
                                    junction=ET.SubElement(junctions, 'junction')
                                    junction.set('id', str(jid))
                                    slices=ET.SubElement(junction, 'slices')
                                    z_prox=height #position of the proximal end of the root segment (initialization)
                                    for iLayer in range(AxialLayers):
                                        slice_iL=ET.SubElement(slices, 'slice')
                                        slice_iL.set('z_prox', str(z_prox))
                                        slice_iL.set('cc1', str(float(soln_C_transi0[iLayer*Ntot+jid])))
                                        z_prox+=height
                            if it+it_resume>9:
                                newtree.write(newpath+'3D_cc_b'+str(b)+','+str(iMaturity)+',t'+str(it+it_resume)+'.xml')
                            else:
                                newtree.write(newpath+'3D_cc_b'+str(b)+','+str(iMaturity)+',t0'+str(it+it_resume)+'.xml')
                            it+=1
                        #end of loop on Print_times
                    
                    if Transient_tracer:
                        #myfile.close()
                        text_file.close()
                elif Apo_Contagion==2:
                    #Solving apoplastic concentrations
                    soln_ApoC = np.linalg.solve(matrix_ApoC,rhs_ApoC) #Solving the equation to get apoplastic & symplastic relative concentrations
                else: # Only Symplastic contagion
                    #Solving apoplastic concentrations
                    soln_SymC = np.linalg.solve(matrix_SymC,rhs_SymC) #Solving the equation to get symplastic relative concentrations
                
            #Resets matrix_C to geometrical factor values
            if Apo_Contagion==2:
                if Sym_Contagion==2: # Apo & Sym contagion
                    if matrix_analysis==1:
                        error('Checking content of matrices')
                    
                    if Transient_tracer: #Removing Dirichlet BC
                        i=0
                        for nid in N_Dirichlet_transi:
                            if sparseM==1:
                                matrix_C2[nid,:]=Stored_matrix_C_transi[i]
                            else:
                                matrix_C[nid,:]=Stored_matrix_C_transi[i]
                            i+=1
                    
                    for i,j,Fjw in Fjw_list:
                        for iLayer in range(AxialLayers):
                            if Fjw[iLayer]>0: #Flow from junction to wall
                                    if sparseM==1:
                                        matrix_C2[iLayer*Ntot+i,iLayer*Ntot+j] -= Fjw[iLayer] #Removing convective term
                                    else:
                                        matrix_C[i][j] -= Fjw[iLayer] #Removing convective term
                                    if sparseM==1:
                                        matrix_C2[iLayer*Ntot+j,iLayer*Ntot+j] += Fjw[iLayer] #Removing convective term
                                    else:
                                        matrix_C[j][j] += Fjw[iLayer] #Removing convective term
                            else: #Flow from wall to junction
                                    if sparseM==1:
                                        matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] -= Fjw[iLayer] #Removing convective term
                                    else:
                                        matrix_C[i][i] -= Fjw[iLayer] #Removing convective term
                                    if sparseM==1:
                                        matrix_C2[iLayer*Ntot+j,iLayer*Ntot+i] += Fjw[iLayer] #Removing convective term
                                    else:
                                        matrix_C[j][i] += Fjw[iLayer] #Removing convective term
                else: #Only Apo contagion
                    for i,j,Fjw in Fjw_list:
                        if Fjw>0: #Flow from junction to wall
                            if i not in N_Dirichlet_ini:
                                matrix_ApoC[i][j] -= Fjw #Removing convective term
                            if j not in N_Dirichlet_ini:
                                matrix_ApoC[j][j] += Fjw #Removing convective term
                        else: #Flow from wall to junction
                            if i not in N_Dirichlet_ini:
                                matrix_ApoC[i][i] -= Fjw #Removing convective term
                            if j not in N_Dirichlet_ini:
                                matrix_ApoC[j][i] += Fjw #Removing convective term
            
            if Sym_Contagion==2: #Convection across plasmodesmata
                if Apo_Contagion==2: #Apo & Sym Contagion
                    for i,j,Fcc in Fcc_list:
                        for iLayer in range(AxialLayers):
                            if Fcc[iLayer]>0: #Flow from j to i
                                    if sparseM==1:
                                        matrix_C2[iLayer*Ntot+j,iLayer*Ntot+j] += Fcc[iLayer] #Removing convective term
                                    else:
                                        matrix_C[j][j] += Fcc[iLayer] #Removing convective term
                                    if sparseM==1:
                                        matrix_C2[iLayer*Ntot+i,iLayer*Ntot+j] -= Fcc[iLayer] #Removing convective term
                                    else:
                                        matrix_C[i][j] -= Fcc[iLayer] #Removing convective term
                            else: #Flow from i to j
                                    if sparseM==1:
                                        matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] -= Fcc[iLayer] #Removing convective term
                                    else:
                                        matrix_C[i][i] -= Fcc[iLayer] #Removing convective term
                                    if sparseM==1:
                                        matrix_C2[iLayer*Ntot+j,iLayer*Ntot+i] += Fcc[iLayer] #Removing convective term
                                    else:
                                        matrix_C[j][i] += Fcc[iLayer] #Removing convective term
                else: #Only Sym contagion
                    for i,j,Fcc in Fcc_list:
                        if Fcc>0: #Flow from j to i
                            if j not in N_Dirichlet_ini:
                                matrix_SymC[j-NwallsJun][j-NwallsJun] += Fcc #Removing convective term
                            if ind_o_cell not in N_Dirichlet_ini:
                                matrix_SymC[i-NwallsJun][j-NwallsJun] -= Fcc #Removing convective term
                        else: #Flow from i to j
                            if i not in N_Dirichlet_ini:
                                matrix_SymC[i-NwallsJun][i-NwallsJun] -= Fcc #Removing convective term
                            if j not in N_Dirichlet_ini:
                                matrix_SymC[j-NwallsJun][i-NwallsJun] += Fcc #Removing convective term
            
            if Apo_Contagion==2 and Sym_Contagion==2:
                for i,j,Fcw,s in Fcw_list:
                    Fcw=-Fcw #Attention, -Fcw was saved
                    for iLayer in range(AxialLayers):
                        if Fcw[iLayer]>0: #Flow from wall to protoplast
                                if D2O1==1:#Solute that moves across membranes like water 
                                    if sparseM==1:
                                        matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] += Fcw[iLayer]
                                    else:
                                        matrix_C[i][i] += Fcw[iLayer] #Removing convective term
                                else: #Solute that moves across membranes independently of water (the membrane is possibly not one) 
                                    if sparseM==1:
                                        matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] += Fcw[iLayer]*(1-s) #Removing convective term
                                    else:
                                        matrix_C[i][i] += Fcw[iLayer]*(1-s) #Removing convective term
                                if D2O1==1:#Solute that moves across membranes like water 
                                    if sparseM==1:
                                        matrix_C2[iLayer*Ntot+j,iLayer*Ntot+i] -= Fcw[iLayer] #Removing convective term
                                    else:
                                        matrix_C[j][i] -= Fcw[iLayer] #Removing convective term
                                else: #Solute that moves across membranes independently of water (the membrane is possibly not one) 
                                    if sparseM==1:
                                        matrix_C2[iLayer*Ntot+j,iLayer*Ntot+i] -= Fcw[iLayer]*(1-s) #Removing convective term
                                    else:
                                        matrix_C[j][i] -= Fcw[iLayer]*(1-s) #Removing convective term
                        else: #Flow from protoplast to wall
                                if D2O1==1:#Solute that moves across membranes like water 
                                    if sparseM==1:
                                        matrix_C2[iLayer*Ntot+j,iLayer*Ntot+j] -= Fcw[iLayer] #Removing convective term
                                    else:
                                        matrix_C[j][j] -= Fcw[iLayer] #Removing convective term
                                else: #Solute that moves across membranes independently of water (the membrane is possibly not one) 
                                    if sparseM==1:
                                        matrix_C2[iLayer*Ntot+j,iLayer*Ntot+j] -= Fcw[iLayer]*(1-s) #Removing convective term
                                    else:
                                        matrix_C[j][j] -= Fcw[iLayer]*(1-s) #Removing convective term
                                if D2O1==1:#Solute that moves across membranes like water 
                                    if sparseM==1:
                                        matrix_C2[iLayer*Ntot+i,iLayer*Ntot+j] += Fcw[iLayer] #Removing convective term
                                    else:
                                        matrix_C[i][j] += Fcw[iLayer] #Removing convective term
                                else: #Solute that moves across membranes independently of water (the membrane is possibly not one) 
                                    if sparseM==1:
                                        matrix_C2[iLayer*Ntot+i,iLayer*Ntot+j] += Fcw[iLayer]*(1-s) #Removing convective term
                                    else:
                                        matrix_C[i][j] += Fcw[iLayer]*(1-s) #Removing convective term
            
            #Removal of ElongNeumann BC type, accounting for "elongation sink"
            if Apo_Contagion==2 and Sym_Contagion==2:
                if sparseM==1:
                    for iLayer in range(AxialLayers):
                        matrix_C2[iLayer*Ntot:(iLayer+1)*Ntot,iLayer*Ntot:(iLayer+1)*Ntot] += spdiags(rhs_e[iLayer*Ntot:(iLayer+1)*Ntot].T, 0, Ntot, Ntot).toarray()
                else:
                    matrix_C += diag(rhs_e[:,0])
            
            #Removal of SoilNeumann BC type
            if Apo_Contagion==2:
                if Sym_Contagion==2:
                    i=0
                    for ind in Borderwall:
                        for iLayer in range(AxialLayers):
                            if ind+Ntot*iLayer in N_SoilNeumann_ini:
                                Q=float(Q_soil[i]) #float(Q_soil[iLayer*len(Borderwall)+i])
                                if Q<0.0: #Solute flowing out of the root
                                    if sparseM==1:
                                        matrix_C2[iLayer*Ntot+ind,iLayer*Ntot+ind]-=Q #Removing convective term
                                    else:
                                        matrix_C[ind][ind]-=Q #Removing convective term
                            i+=1
                    for ind in Borderjunction:
                        for iLayer in range(AxialLayers):
                            if ind+Ntot*iLayer in N_SoilNeumann_ini:
                                Q=float(Q_soil[i]) #float(Q_soil[iLayer*len(Borderwall)+i])
                                if Q<0.0: #Solute flowing out of the root
                                    if sparseM==1:
                                        matrix_C2[iLayer*Ntot+ind,iLayer*Ntot+ind]-=Q #Removing convective term
                                    else:
                                        matrix_C[ind][ind]-=Q #Removing convective term
                            i+=1
                else:
                    i=0
                    for ind in Borderwall:
                        for iLayer in range(AxialLayers):
                            if ind+Ntot*iLayer in N_SoilNeumann_ini:
                                Q=float(Q_soil[i]) #float(Q_soil[iLayer*len(Borderwall)+i])
                                if Q<0.0:
                                    matrix_ApoC[ind][ind]-=Q #Removing convective term
                            i+=1
                    for ind in Borderjunction:
                        for iLayer in range(AxialLayers):
                            if ind+Ntot*iLayer in N_SoilNeumann_ini:
                                Q=float(Q_soil[i]) #float(Q_soil[iLayer*len(Borderwall)+i])
                                if Q<0.0:
                                    matrix_ApoC[ind][ind]-=Q #Removing convective term
                            i+=1
            
            #Removal of axial convective terms
            if Sym_Contagion==2 and Apo_Contagion==2:
                #Removal of axial solute convection terms, except at axial BCs
                i=0 #node index (0 to Ntot-1)
                for Faxial in Faxial_sym_apo_list:
                    for iLayer in range(AxialLayers-1):
                        Q=Faxial[iLayer] #Positive when pressure i (distal) > pressure i+Ntot (proximal)
                        if Q>0: #Flow from i upwards
                            if sparseM==1:
                                    matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] += Q
                                    matrix_C2[(iLayer+1)*Ntot+i,iLayer*Ntot+i] -= Q
                            else:
                                    matrix_C[iLayer*Ntot+i][iLayer*Ntot+i] += Q
                                    matrix_C[(iLayer+1)*Ntot+i][iLayer*Ntot+i] -= Q
                        else: #Flow from i+Ntot downwards
                            if sparseM==1:
                                    matrix_C2[(iLayer+1)*Ntot+i,(iLayer+1)*Ntot+i] -= Q  #changed sign 12/2/2019
                                    matrix_C2[iLayer*Ntot+i,(iLayer+1)*Ntot+i] += Q #changed sign 12/2/2019
                            else:
                                    matrix_C[(iLayer+1)*Ntot+i][(iLayer+1)*Ntot+i] -= Q  #changed sign 12/2/2019
                                    matrix_C[iLayer*Ntot+i][(iLayer+1)*Ntot+i] += Q #changed sign 12/2/2019
                    i+=1
                
                i=NwallsJun #cell node index (NwallsJun to Ntot-1)
                for Faxial in Faxial_TM_list:
                    for iLayer in range(AxialLayers-1):
                        Q=Faxial[iLayer] #Positive when pressure i (distal) > pressure i+Ntot (proximal)
                        if Q>0: #Flow from i upwards
                            if sparseM==1:
                                if D2O1==1:#Solute that moves across membranes like water 
                                    matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] += Q
                                    matrix_C2[(iLayer+1)*Ntot+i,iLayer*Ntot+i] -= Q
                                else:#Solute that moves across membranes independently of water (the membrane is possibly not one) 
                                    if PileUp==2:
                                        matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] += Q*(1-s_cells[i-NwallsJun,MaturL[iLayer]])
                                        matrix_C2[(iLayer+1)*Ntot+i,iLayer*Ntot+i] -= Q*(1-s_cells[i-NwallsJun,MaturL[iLayer]])
                                    else:
                                        matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] += Q*(1-s_cells[i-NwallsJun])
                                        matrix_C2[(iLayer+1)*Ntot+i,iLayer*Ntot+i] -= Q*(1-s_cells[i-NwallsJun])
                            else:
                                    matrix_C[iLayer*Ntot+i][iLayer*Ntot+i] += Q
                                    matrix_C[(iLayer+1)*Ntot+i][iLayer*Ntot+i] -= Q
                        else: #Flow from i+Ntot downwards
                            if sparseM==1:
                                if D2O1==1:#Solute that moves across membranes like water 
                                    matrix_C2[(iLayer+1)*Ntot+i,(iLayer+1)*Ntot+i] -= Q  #changed sign 12/2/2019
                                    matrix_C2[iLayer*Ntot+i,(iLayer+1)*Ntot+i] += Q #changed sign 12/2/2019
                                else:
                                    if PileUp==2:
                                        matrix_C2[(iLayer+1)*Ntot+i,(iLayer+1)*Ntot+i] -= Q*(1-s_cells[i-NwallsJun,MaturL[iLayer]])  #changed sign 12/2/2019
                                        matrix_C2[iLayer*Ntot+i,(iLayer+1)*Ntot+i] += Q*(1-s_cells[i-NwallsJun,MaturL[iLayer]]) #changed sign 12/2/2019
                                    else:
                                        matrix_C2[(iLayer+1)*Ntot+i,(iLayer+1)*Ntot+i] -= Q*(1-s_cells[i-NwallsJun])  #changed sign 12/2/2019
                                        matrix_C2[iLayer*Ntot+i,(iLayer+1)*Ntot+i] += Q*(1-s_cells[i-NwallsJun]) #changed sign 12/2/2019
                            else:
                                    matrix_C[(iLayer+1)*Ntot+i][(iLayer+1)*Ntot+i] -= Q  #changed sign 12/2/2019
                                    matrix_C[iLayer*Ntot+i][(iLayer+1)*Ntot+i] += Q #changed sign 12/2/2019
                    i+=1
                #for Faxial in Faxial_list:
                #    for iLayer in range(AxialLayers-1):
                #        Q=Faxial[iLayer] #Positive when pressure i (distal) > pressure i+Ntot (proximal)
                #        if Q>0: #Flow from i upwards
                #            if sparseM==1:
                #                    matrix_C2[iLayer*Ntot+i,iLayer*Ntot+i] += Q
                #                    matrix_C2[(iLayer+1)*Ntot+i,iLayer*Ntot+i] -= Q
                #            else:
                #                    matrix_C[iLayer*Ntot+i][iLayer*Ntot+i] += Q
                #                    matrix_C[(iLayer+1)*Ntot+i][iLayer*Ntot+i] -= Q
                #        else: #Flow from i+Ntot downwards
                #            if sparseM==1:
                #                    matrix_C2[(iLayer+1)*Ntot+i,(iLayer+1)*Ntot+i] -= Q  #changed sign 12/2/2019
                #                    matrix_C2[iLayer*Ntot+i,(iLayer+1)*Ntot+i] += Q #changed sign 12/2/2019
                #            else:
                #                    matrix_C[(iLayer+1)*Ntot+i][(iLayer+1)*Ntot+i] -= Q  #changed sign 12/2/2019
                #                    matrix_C[iLayer*Ntot+i][(iLayer+1)*Ntot+i] += Q #changed sign 12/2/2019
                #    i+=1
                
                #Removal of AxialNeumann BC
                if XylOK>0:
                    i=0
                    for cid in listxyl:
                        Q_dist=float(Q_xyl[0][i]) #Positive upwards
                        Q_prox=float(Q_xyl[1][i])
                        if Q_prox>0: #Flow from proximal xylem vessel toward shoot
                            if sparseM==1:
                                    matrix_C2[-Ntot+cid,-Ntot+cid] += Q_prox
                            else:
                                    matrix_C[-Ntot+cid][-Ntot+cid] += Q_prox
                        if Q_dist<0: #Flow from distal xylem vessel toward root tip
                            if sparseM==1:
                                    matrix_C2[cid,cid] -= Q_dist
                            else:
                                    matrix_C[cid][cid] -= Q_dist
                        i+=1
                
                #Removal of Phloem BC for solute
                i=0
                for cid in listsieve:
                    if Q_sieve[1][i]>0: #Proximal phloem BC shootward
                        if sparseM==1:
                            matrix_C2[-Ntot+cid,-Ntot+cid] += Q_sieve[1][i]
                        else:
                            matrix_C[-Ntot+cid][-Ntot+cid] += Q_sieve[1][i]
                    if Q_sieve[0][i]<0: #Distal phloem BC tipward
                        if sparseM==1:
                            matrix_C2[cid,cid] -= Q_sieve[0][i]
                        else:
                            matrix_C[cid][cid] -= Q_sieve[0][i]
                    i+=1                
            
            ####################################
            ## Creates .vtk file for Paraview ##
            ####################################
            
            if Paraview==1 or ParTrack==1:
                if Sym_Contagion==1:
                    Sym_Zombies=[]
                    for source in Sym_source_ini_range:
                        Sym_Zombies.append(int(source.get("id")))
                    iZombie=0
                    while not iZombie == size(Sym_Zombies):
                        itemp=0
                        for cid in Cell_connec[int(Sym_Zombies[iZombie])][0:int(nCell_connec[int(Sym_Zombies[iZombie])])]:
                            if Cell_connec_flow[int(Sym_Zombies[iZombie])][itemp] == -1 and (cid not in Sym_Zombies): #Infection
                                if cid in Sym_Immune:
                                    print(cid,': "You shall not pass!"')
                                else:
                                    Sym_Zombies.append(cid)
                                    print(cid,': "Aaargh!"      Zombie count:', size(Sym_Zombies)+1)
                            itemp+=1
                        iZombie+=1
                    print('End of the propagation. Survivor count:', Ncells-size(Sym_Zombies)-1)
                    for cid in Sym_Target:
                        if cid in Sym_Zombies:
                            print('Target '+ str(cid) +' down. XXX')
                        else:
                            print('Target '+ str(cid) +' missed!')
                    if Sym_Target[0] in Sym_Zombies:
                        if Sym_Target[1] in Sym_Zombies:
                            Hydropatterning[iMaturity][count]=0 #Both targets reached
                        else:
                            Hydropatterning[iMaturity][count]=1 #Target1 reached only
                    elif Sym_Target[1] in Sym_Zombies:
                        Hydropatterning[iMaturity][count]=2 #Target2 reached only
                    else:
                        Hydropatterning[iMaturity][count]=-1 #Not target reached
                    
                    
                    for Print_layer in Printrange:
                        iLayer=int(Print_layer.get("value"))
                        if PileUp==2:
                            b=int(BarriersL[iLayer])
                        else:
                            b=Barrier
                        text_file = open(newpath+"Sym_Contagion_bottomb"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "w")
                        with open(newpath+"Sym_Contagion_bottomb"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "a") as myfile:
                            myfile.write("# vtk DataFile Version 4.0 \n")
                            myfile.write("Contaminated symplastic space geometry \n")
                            myfile.write("ASCII \n")
                            myfile.write(" \n")
                            myfile.write("DATASET UNSTRUCTURED_GRID \n")
                            myfile.write("POINTS "+str(len(ThickWalls))+" float \n")
                            for ThickWallNode in ThickWalls:
                                myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(Cumheights[iLayer,0]/10) + " \n") #height/200
                            myfile.write(" \n")
                            myfile.write("CELLS " + str(len(Sym_Zombies)) + " " + str(int(len(Sym_Zombies)+sum(nCell2ThickWalls[Sym_Zombies]))) + " \n") #The number of cells corresponds to the number of intercellular spaces
                            Sym_Contagion_order=zeros((Ncells,1))
                            temp=0
                            for cid in Sym_Zombies:
                                n=int(nCell2ThickWalls[cid]) #Total number of thick wall nodes around the protoplast
                                Polygon=Cell2ThickWalls[cid][:n]
                                ranking=list()
                                ranking.append(int(Polygon[0]))
                                ranking.append(ThickWalls[int(ranking[0])][5])
                                ranking.append(ThickWalls[int(ranking[0])][6])
                                for id1 in range(1,n):
                                    wid1=ThickWalls[int(ranking[id1])][5]
                                    wid2=ThickWalls[int(ranking[id1])][6]
                                    if wid1 not in ranking:
                                        ranking.append(wid1)
                                    if wid2 not in ranking:
                                        ranking.append(wid2)
                                string=str(n)
                                for id1 in ranking:
                                    string=string+" "+str(int(id1))
                                myfile.write(string + " \n")
                                Sym_Contagion_order[cid]=temp
                                temp+=1
                            myfile.write(" \n")
                            myfile.write("CELL_TYPES " + str(len(Sym_Zombies)) + " \n")
                            for i in range(len(Sym_Zombies)):
                                myfile.write("6 \n") #Triangle-strip cell type
                            myfile.write(" \n")
                            myfile.write("POINT_DATA " + str(len(ThickWalls)) + " \n")
                            myfile.write("SCALARS Sym_Contagion_order_(#) float \n")
                            myfile.write("LOOKUP_TABLE default \n")
                            for ThickWallNode in ThickWalls:
                                cellnumber1=ThickWallNode[2]-NwallsJun
                                myfile.write(str(int(Sym_Contagion_order[int(cellnumber1)])) + " \n") #Flow rate from wall (non junction) to cell    min(sath1,max(satl1,  ))
                        myfile.close()
                        text_file.close()
                
                if Apo_Contagion==1:
                    Apo_w_Zombies=N_Dirichlet_ini
                    iZombie=0
                    while not iZombie == size(Apo_w_Zombies):
                        id1=Apo_w_Zombies[iZombie]
                        for id2 in Apo_connec_flow[id1][0:nApo_connec_flow[id1]]:
                            if id2 not in Apo_w_Zombies: #Infection
                                if id2 in Apo_w_Immune:
                                    print(id2,': "You shall not pass!"')
                                else:
                                    Apo_w_Zombies.append(id2)
                                    print(id2,': "Aaargh!"      Zombie count:', size(Apo_w_Zombies))
                        iZombie+=1
                    print('End of the propagation. Survivor count:', NwallsJun-size(Apo_w_Zombies))
                    temp=0
                    for wid in Apo_w_Target:
                        if wid in Apo_w_Zombies:
                            temp+=1
                            print('Target '+ str(wid) +' down. XXX')
                        else:
                            print('Target '+ str(wid) +' missed!')
                    Hydrotropism[iMaturity][count]=float(temp)/size(Apo_w_Target) #0: No apoplastic target reached; 1: All apoplastic targets reached
                    
                    for Print_layer in Printrange:
                        iLayer=int(Print_layer.get("value"))
                        if PileUp==2:
                            b=int(BarriersL[iLayer])
                        else:
                            b=Barrier
                        text_file = open(newpath+"Apo_Contagion_bottomb"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "w")
                        with open(newpath+"Apo_Contagion_bottomb"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "a") as myfile:
                            myfile.write("# vtk DataFile Version 4.0 \n")
                            myfile.write("Contaminated Apoplastic space geometry \n")
                            myfile.write("ASCII \n")
                            myfile.write(" \n")
                            myfile.write("DATASET UNSTRUCTURED_GRID \n")
                            myfile.write("POINTS "+str(len(ThickWallsX))+" float \n")
                            for ThickWallNodeX in ThickWallsX:
                                myfile.write(str(ThickWallNodeX[1]) + " " + str(ThickWallNodeX[2]) + " " + str(Cumheights[iLayer,0]/10) + " \n")
                            myfile.write(" \n")
                            myfile.write("CELLS " + str(int(NwallsJun+Nwalls-len(list_ghostwalls)*2-len(list_ghostjunctions))) + " " + str(int(2*Nwalls*5-len(list_ghostwalls)*10+sum(nWall2NewWallX[Nwalls:])+NwallsJun-Nwalls+2*len(Wall2NewWallX[Nwalls:])-nGhostJunction2Wall-len(list_ghostjunctions))) + " \n") #The number of cells corresponds to the number of lines in ThickWalls (if no ghost wall & junction)
                            i=0
                            for PolygonX in ThickWallPolygonX:
                                if floor(i/2) not in list_ghostwalls:
                                    myfile.write("4 " + str(int(PolygonX[0])) + " " + str(int(PolygonX[1])) + " " + str(int(PolygonX[2])) + " " + str(int(PolygonX[3])) + " \n")
                                i+=1
                            j=Nwalls
                            for PolygonX in Wall2NewWallX[Nwalls:]: #"junction" polygons
                                #Would need to order them based on x or y position to make sure display fully covers the surface (but here we try a simpler not so good solution instead)
                                if j not in list_ghostjunctions:
                                    string=str(int(nWall2NewWallX[j]+2)) #Added +2 so that the first and second nodes could be added again at the end (trying to fill the polygon better)
                                    for id1 in range(int(nWall2NewWallX[j])):
                                        string=string+" "+str(int(PolygonX[id1]))
                                    string=string+" "+str(int(PolygonX[0]))+" "+str(int(PolygonX[1])) #Adding the 1st and 2nd nodes again to the end
                                    myfile.write(string + " \n")
                                j+=1
                            myfile.write(" \n")
                            myfile.write("CELL_TYPES " + str(NwallsJun+Nwalls-len(list_ghostwalls)*2-len(list_ghostjunctions)) + " \n")
                            i=0
                            for PolygonX in ThickWallPolygonX:
                                if floor(i/2) not in list_ghostwalls:
                                    myfile.write("7 \n") #Polygon cell type (wall)
                                i+=1
                            j=Nwalls
                            for PolygonX in Wall2NewWallX[Nwalls:]:
                                if j not in list_ghostjunctions:
                                    myfile.write("6 \n") #Triangle-strip cell type (wall junction)
                                j+=1
                            myfile.write(" \n")
                            myfile.write("POINT_DATA " + str(len(ThickWallsX)) + " \n")
                            myfile.write("SCALARS Apo_Contagion_order_(#) float \n")
                            myfile.write("LOOKUP_TABLE default \n")
                            Apo_Contagion_order=zeros((NwallsJun,1))+int(len(Apo_w_Zombies)*1.6)
                            temp=0
                            for wid in Apo_w_Zombies:
                                Apo_Contagion_order[wid]=temp
                                temp+=1
                            NewApo_Contagion_order=zeros((len(ThickWallsX),1))
                            j=0
                            for PolygonX in Wall2NewWallX:
                                for id1 in range(int(nWall2NewWallX[j])):
                                    NewApo_Contagion_order[int(PolygonX[id1])]=Apo_Contagion_order[j]
                                j+=1
                            for i in range(len(ThickWallsX)):
                                myfile.write(str(float(NewApo_Contagion_order[i])) + " \n")
                        myfile.close()
                        text_file.close()
                    
            if Paraview==1:
                if ParaviewWP==1: #2D visualization of walls pressure potentials
                    for Print_layer in Printrange:
                        iLayer=int(Print_layer.get("value"))
                        if PileUp==2:
                            b=int(BarriersL[iLayer])
                        else:
                            b=Barrier
                        text_file = open(newpath+"Walls2Db"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "w")
                        with open(newpath+"Walls2Db"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "a") as myfile:
                            myfile.write("# vtk DataFile Version 4.0 \n")     #("Purchase Amount: %s" % TotalAmount)
                            myfile.write("Wall geometry 2D \n")
                            myfile.write("ASCII \n")
                            myfile.write(" \n")
                            myfile.write("DATASET UNSTRUCTURED_GRID \n")
                            myfile.write("POINTS "+str(len(G.node))+" float \n")
                            for node in G:
                                myfile.write(str(float(position[node][0])) + " " + str(float(position[node][1])) + " " + str(Cumheights[iLayer,0]/10) + " \n")
                            myfile.write(" \n")
                            myfile.write("CELLS " + str(Nwalls*2-len(list_ghostwalls)*2) + " " + str(Nwalls*6-len(list_ghostwalls)*6) + " \n") #len(G.node)
                            for node, edges in G.adjacency_iter():
                                i=indice[node]
                                if i not in list_ghostwalls:
                                    for neighbour, eattr in edges.items(): #Loop on connections (edges)
                                        j=indice[neighbour]
                                        if j>i and eattr['path']=='wall':
                                            #print(nx.get_node_attributes(edges,'path'))
                                            myfile.write(str(2) + " " + str(i) + " " + str(j) + " \n")
                            myfile.write(" \n")
                            myfile.write("CELL_TYPES " + str(Nwalls*2-len(list_ghostwalls)*2) + " \n") #The number of nodes corresponds to the number of wall to wall connections.... to be checked, might not be generality
                            for node, edges in G.adjacency_iter():
                                i=indice[node]
                                if i not in list_ghostwalls:
                                    for neighbour, eattr in edges.items(): #Loop on connections (edges)
                                        j=indice[neighbour]
                                        if j>i and eattr['path']=='wall':
                                            #print(nx.get_node_attributes(edges,'path'))
                                            myfile.write(str(3) + " \n") #Line cell type
                            myfile.write(" \n")
                            myfile.write("POINT_DATA " + str(len(G.node)) + " \n")
                            myfile.write("SCALARS Wall_pressure float \n")
                            myfile.write("LOOKUP_TABLE default \n")
                            for node in G:
                                myfile.write(str(float(soln[iLayer*Ntot+node])) + " \n") #Line cell type      min(sath0,max(satl0,   ))
                        myfile.close()
                        text_file.close()
                
                if ParaviewWP==1 and ParaviewCP: #2D visualization of walls & cells osmotic potentials
                    for Print_layer in Printrange:
                        iLayer=int(Print_layer.get("value"))
                        if PileUp==2:
                            b=int(BarriersL[iLayer])
                            iM=int(MaturL[iLayer])
                        else:
                            b=Barrier
                            iM=iMaturity
                        text_file = open(newpath+"WallsOsAndCellsOs2Db"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "w")
                        with open(newpath+"WallsOsAndCellsOs2Db"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "a") as myfile:
                            myfile.write("# vtk DataFile Version 4.0 \n")     #("Purchase Amount: %s" % TotalAmount)
                            myfile.write("Wall geometry 2D \n")
                            myfile.write("ASCII \n")
                            myfile.write(" \n")
                            myfile.write("DATASET UNSTRUCTURED_GRID \n")
                            myfile.write("POINTS "+str(len(G.node))+" float \n")
                            for node in G:
                                myfile.write(str(float(position[node][0])) + " " + str(float(position[node][1])) + " " + str(Cumheights[iLayer,0]/10) + " \n")
                            myfile.write(" \n")                                     
                            myfile.write("CELLS " + str(Nwalls*2-len(list_ghostwalls)*2+Ncells) + " " + str(Nwalls*6-len(list_ghostwalls)*6+Ncells*2) + " \n") #
                            for node, edges in G.adjacency_iter():
                                i=indice[node]
                                if i not in list_ghostwalls:
                                    for neighbour, eattr in edges.items(): #Loop on connections (edges)
                                        j=indice[neighbour]
                                        if j>i and eattr['path']=='wall':
                                            #print(nx.get_node_attributes(edges,'path'))
                                            myfile.write(str(2) + " " + str(i) + " " + str(j) + " \n")
                                if i>=NwallsJun: #Cell node
                                    myfile.write("1 " + str(i) + " \n")
                            myfile.write(" \n")
                            myfile.write("CELL_TYPES " + str(Nwalls*2-len(list_ghostwalls)*2+Ncells) + " \n") #
                            for node, edges in G.adjacency_iter():
                                i=indice[node]
                                if i not in list_ghostwalls:
                                    for neighbour, eattr in edges.items(): #Loop on connections (edges)
                                        j=indice[neighbour]
                                        if j>i and eattr['path']=='wall':
                                            #print(nx.get_node_attributes(edges,'path'))
                                            myfile.write(str(3) + " \n") #Line cell type
                                if i>=NwallsJun: #Cell node
                                    myfile.write("1 \n")
                            myfile.write(" \n")
                            myfile.write("POINT_DATA " + str(len(G.node)) + " \n")
                            myfile.write("SCALARS Wall_and_Cell_osmotic_pot float \n")
                            myfile.write("LOOKUP_TABLE default \n")
                            for node, edges in G.adjacency_iter():
                                i=indice[node] #Node ID number
                                if i<Nwalls: #Wall node
                                    if Ana_Os and AxialLayers>1:
                                        myfile.write(str(float(Os_walls[i,iM])) + " \n")
                                    else:
                                        myfile.write(str(float(Os_walls[i,0])) + " \n")
                                elif i<NwallsJun: #Junction node
                                    myfile.write(str(float(0.0)) + " \n")
                                else: #Cell node
                                    if PileUp==2:
                                        myfile.write(str(float(Os_cells[i-NwallsJun,iM])) + " \n")
                                    else:
                                        myfile.write(str(float(Os_cells[i-NwallsJun,0])) + " \n")
                        myfile.close()
                        text_file.close()
                    

                
                if ParaviewWP==1 and ParaviewCP==1: #2D visualization of walls & cells water potentials
                    for Print_layer in Printrange:
                        iLayer=int(Print_layer.get("value"))
                        if PileUp==2:
                            b=int(BarriersL[iLayer])
                        else:
                            b=Barrier
                        text_file = open(newpath+"WallsAndCells2Db"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "w")
                        with open(newpath+"WallsAndCells2Db"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "a") as myfile:
                            myfile.write("# vtk DataFile Version 4.0 \n")     #("Purchase Amount: %s" % TotalAmount)
                            myfile.write("Water potential distribution in cells and walls 2D \n")
                            myfile.write("ASCII \n")
                            myfile.write(" \n")
                            myfile.write("DATASET UNSTRUCTURED_GRID \n")
                            myfile.write("POINTS "+str(len(G.node))+" float \n")
                            for node in G:
                                myfile.write(str(float(position[node][0])) + " " + str(float(position[node][1])) + " " + str(Cumheights[iLayer,0]/10) + " \n")
                            myfile.write(" \n")
                            myfile.write("CELLS " + str(Nwalls*2-len(list_ghostwalls)*2+Ncells) + " " + str(Nwalls*6-len(list_ghostwalls)*6+Ncells*2) + " \n") #
                            for node, edges in G.adjacency_iter():
                                i=indice[node]
                                if i not in list_ghostwalls:
                                    for neighbour, eattr in edges.items(): #Loop on connections (edges)
                                        j=indice[neighbour]
                                        if j>i and eattr['path']=='wall':
                                            #print(nx.get_node_attributes(edges,'path'))
                                            myfile.write(str(2) + " " + str(i) + " " + str(j) + " \n")
                                if i>=NwallsJun: #Cell node
                                    myfile.write("1 " + str(i) + " \n")
                            myfile.write(" \n")
                            myfile.write("CELL_TYPES " + str(Nwalls*2-len(list_ghostwalls)*2+Ncells) + " \n") #
                            for node, edges in G.adjacency_iter():
                                i=indice[node]
                                if i not in list_ghostwalls:
                                    for neighbour, eattr in edges.items(): #Loop on connections (edges)
                                        j=indice[neighbour]
                                        if j>i and eattr['path']=='wall':
                                            #print(nx.get_node_attributes(edges,'path'))
                                            myfile.write(str(3) + " \n") #Line cell type
                                if i>=NwallsJun: #Cell node
                                    myfile.write("1 \n")
                            myfile.write(" \n")
                            myfile.write("POINT_DATA " + str(len(G.node)) + " \n")
                            myfile.write("SCALARS pressure float \n")
                            myfile.write("LOOKUP_TABLE default \n")
                            for node in G:
                                myfile.write(str(float(soln[iLayer*Ntot+node])) + " \n") #Line cell type
                        myfile.close()
                        text_file.close()
                
                if ParaviewCP==1: #2D visualization of cells water potentials
                    for Print_layer in Printrange:
                        iLayer=int(Print_layer.get("value"))
                        if PileUp==2:
                            b=int(BarriersL[iLayer])
                        else:
                            b=Barrier
                        text_file = open(newpath+"Cells2Db"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "w")
                        with open(newpath+"Cells2Db"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "a") as myfile:
                            myfile.write("# vtk DataFile Version 4.0 \n")     #("Purchase Amount: %s" % TotalAmount)
                            myfile.write("Pressure potential distribution in cells 2D \n")
                            myfile.write("ASCII \n")
                            myfile.write(" \n")
                            myfile.write("DATASET UNSTRUCTURED_GRID \n")
                            myfile.write("POINTS "+str(len(G.node))+" float \n")
                            for node in G:
                                myfile.write(str(float(position[node][0])) + " " + str(float(position[node][1])) + " " + str(Cumheights[iLayer,0]/10) + " \n")
                            myfile.write(" \n")
                            myfile.write("CELLS " + str(Ncells) + " " + str(Ncells*2) + " \n") #
                            for node, edges in G.adjacency_iter():
                                i=indice[node]
                                if i>=NwallsJun: #Cell node
                                    myfile.write("1 " + str(i) + " \n")
                            myfile.write(" \n")
                            myfile.write("CELL_TYPES " + str(Ncells) + " \n") #
                            for node, edges in G.adjacency_iter():
                                i=indice[node]
                                if i>=NwallsJun: #Cell node
                                    myfile.write("1 \n")
                            myfile.write(" \n")
                            myfile.write("POINT_DATA " + str(len(G.node)) + " \n")
                            myfile.write("SCALARS Cell_pressure float \n")
                            myfile.write("LOOKUP_TABLE default \n")
                            for node in G:
                                myfile.write(str(float(soln[iLayer*Ntot+node])) + " \n") #Line cell type      min(sath01,max(satl01,   ))
                        myfile.close()
                        text_file.close()
                    
                
                if ParaviewMF==1: #3D visualization of membrane fluxes
                    for Print_layer in Printrange:
                        iLayer=int(Print_layer.get("value"))
                        if PileUp==2:
                            b=int(BarriersL[iLayer])
                        else:
                            b=Barrier
                        text_file = open(newpath+"Membranes3Db"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "w")
                        with open(newpath+"Membranes3Db"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "a") as myfile:
                            myfile.write("# vtk DataFile Version 4.0 \n")
                            myfile.write("Membranes geometry 3D \n")
                            myfile.write("ASCII \n")
                            myfile.write(" \n")
                            myfile.write("DATASET UNSTRUCTURED_GRID \n")
                            myfile.write("POINTS "+str(len(ThickWalls)*2)+" float \n")
                            for ThickWallNode in ThickWalls:
                                myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(Cumheights[iLayer,0]/10) + " \n")
                            if PileUp==2 or AxialLayers>1:
                                for ThickWallNode in ThickWalls:
                                    if iLayer<AxialLayers-1:
                                        myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(Cumheights[iLayer+1,0]/10) + " \n")
                                    else:
                                        myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str((Cumheights[iLayer,0]+Cumheights[iLayer,0]-Cumheights[iLayer-1,0])/10) + " \n")
                            else:
                                for ThickWallNode in ThickWalls:
                                    myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(height/10) + " \n")
                            myfile.write(" \n")
                            myfile.write("CELLS " + str(len(ThickWalls)-len(list_ghostwalls)*4) + " " + str(len(ThickWalls)*5-len(list_ghostwalls)*20) + " \n") #The number of cells corresponds to the number of lines in ThickWalls
                            for ThickWallNode in ThickWalls:
                                if ThickWallNode[1]>=Nwalls: #wall that is a junction
                                    if ThickWalls[int(ThickWallNode[5])][1] not in list_ghostwalls:
                                        myfile.write("4 " + str(int(ThickWallNode[0])) + " " + str(int(ThickWallNode[5])) + " " + str(int(ThickWallNode[5])+len(ThickWalls)) + " " + str(int(ThickWallNode[0])+len(ThickWalls)) + " \n") #All points were repeated twice (once at z=0 and once at z=height), so adding len(ThickWalls) is the same point at z=height
                                    if ThickWalls[int(ThickWallNode[6])][1] not in list_ghostwalls:
                                        myfile.write("4 " + str(int(ThickWallNode[0])) + " " + str(int(ThickWallNode[6])) + " " + str(int(ThickWallNode[6])+len(ThickWalls)) + " " + str(int(ThickWallNode[0])+len(ThickWalls)) + " \n")
                            myfile.write(" \n")
                            myfile.write("CELL_TYPES " + str(len(ThickWalls)-len(list_ghostwalls)*4) + " \n")
                            for ThickWallNode in ThickWalls:
                                if ThickWallNode[1]>=Nwalls: #wall that is a junction
                                    if ThickWalls[int(ThickWallNode[5])][1] not in list_ghostwalls:
                                        myfile.write("9 \n") #Quad cell type
                                    if ThickWalls[int(ThickWallNode[6])][1] not in list_ghostwalls:
                                        myfile.write("9 \n") #Quad cell type
                            myfile.write(" \n")
                            myfile.write("POINT_DATA " + str(len(ThickWalls)*2) + " \n")
                            myfile.write("SCALARS TM_flux_(m/s) float \n")
                            myfile.write("LOOKUP_TABLE default \n")
                            for ThickWallNode in ThickWalls:
                                if ThickWallNode[0]<len(MembraneFlowDensity):
                                    myfile.write(str(float(MembraneFlowDensity[int(ThickWallNode[0])][iLayer])/sperd/cmperm) + " \n") #Flow rate from wall (non junction) to cell   min(sath1,max(satl1,  ))
                                else:
                                    myfile.write(str(float((MembraneFlowDensity[int(ThickWallNode[5])][iLayer]+MembraneFlowDensity[int(ThickWallNode[6])][iLayer])/2)/sperd/cmperm) + " \n") #Flow rate from junction wall to cell is the average of the 2 neighbouring wall flow rates   min(sath1,max(satl1,  ))
                            for ThickWallNode in ThickWalls:
                                if ThickWallNode[0]<len(MembraneFlowDensity):
                                    myfile.write(str(float(MembraneFlowDensity[int(ThickWallNode[0])][iLayer])/sperd/cmperm) + " \n") #Flow rate from wall (non junction) to cell   min(sath1,max(satl1,  ))
                                else:
                                    myfile.write(str(float((MembraneFlowDensity[int(ThickWallNode[5])][iLayer]+MembraneFlowDensity[int(ThickWallNode[6])][iLayer])/2)/sperd/cmperm) + " \n") #Flow rate from junction wall to cell is the average of the 2 neighbouring wall flow rates   min(sath1,max(satl1,  ))
                        myfile.close()
                        text_file.close()
                
                if ParaviewWF==1: #Wall flow density data
                    #maxWallFlowDensity=0.0
                    #for ir in range(int(len(WallFlowDensity))):
                    #    maxWallFlowDensity=max(maxWallFlowDensity,max(abs(WallFlowDensity[ir][2:])))
                    #sath2=maxWallFlowDensity*color_threshold #(1-(1-color_threshold)/2)
                    #satl2=0.0
                    for Print_layer in Printrange:
                        iLayer=int(Print_layer.get("value"))
                        if PileUp==2:
                            b=int(BarriersL[iLayer])
                        else:
                            b=Barrier
                        text_file = open(newpath+"WallsThick3D_bottomb"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "w")
                        with open(newpath+"WallsThick3D_bottomb"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "a") as myfile:
                            myfile.write("# vtk DataFile Version 4.0 \n")
                            myfile.write("Wall geometry 3D including thickness bottom \n")
                            myfile.write("ASCII \n")
                            myfile.write(" \n")
                            myfile.write("DATASET UNSTRUCTURED_GRID \n")
                            myfile.write("POINTS "+str(len(ThickWallsX))+" float \n")
                            for ThickWallNodeX in ThickWallsX:
                                myfile.write(str(ThickWallNodeX[1]) + " " + str(ThickWallNodeX[2]) + " " + str(Cumheights[iLayer,0]/10) + " \n")
                            myfile.write(" \n")
                            myfile.write("CELLS " + str(int(NwallsJun+Nwalls-len(list_ghostwalls)*2-len(list_ghostjunctions))) + " " + str(int(2*Nwalls*5-len(list_ghostwalls)*10+sum(nWall2NewWallX[Nwalls:])+NwallsJun-Nwalls+2*len(Wall2NewWallX[Nwalls:])-nGhostJunction2Wall-len(list_ghostjunctions))) + " \n") #The number of cells corresponds to the number of lines in ThickWalls (if no ghost wall & junction)
                            i=0
                            for PolygonX in ThickWallPolygonX:
                                if floor(i/2) not in list_ghostwalls:
                                    myfile.write("4 " + str(int(PolygonX[0])) + " " + str(int(PolygonX[1])) + " " + str(int(PolygonX[2])) + " " + str(int(PolygonX[3])) + " \n")
                                i+=1
                            j=Nwalls
                            for PolygonX in Wall2NewWallX[Nwalls:]: #"junction" polygons
                                #Would need to order them based on x or y position to make sure display fully covers the surface (but here we try a simpler not so good solution instead)
                                if j not in list_ghostjunctions:
                                    string=str(int(nWall2NewWallX[j]+2)) #Added +2 so that the first and second nodes could be added again at the end (trying to fill the polygon better)
                                    for id1 in range(int(nWall2NewWallX[j])):
                                        string=string+" "+str(int(PolygonX[id1]))
                                    string=string+" "+str(int(PolygonX[0]))+" "+str(int(PolygonX[1])) #Adding the 1st and 2nd nodes again to the end
                                    myfile.write(string + " \n")
                                j+=1
                            myfile.write(" \n")
                            myfile.write("CELL_TYPES " + str(NwallsJun+Nwalls-len(list_ghostwalls)*2-len(list_ghostjunctions)) + " \n")
                            i=0
                            for PolygonX in ThickWallPolygonX:
                                if floor(i/2) not in list_ghostwalls:
                                    myfile.write("7 \n") #Polygon cell type (wall)
                                i+=1
                            j=Nwalls
                            for PolygonX in Wall2NewWallX[Nwalls:]:
                                if j not in list_ghostjunctions:
                                    myfile.write("6 \n") #Triangle-strip cell type (wall junction)
                                j+=1
                            myfile.write(" \n")
                            myfile.write("POINT_DATA " + str(len(ThickWallsX)) + " \n")
                            myfile.write("SCALARS Apo_flux_(m/s) float \n")
                            myfile.write("LOOKUP_TABLE default \n")
                            NewWallFlowDensity=zeros((len(ThickWallsX),2))
                            i=0
                            for PolygonX in ThickWallPolygonX:
                                for id1 in range(4):
                                    if abs(float(WallFlowDensity[i][2][iLayer]))>min(NewWallFlowDensity[int(PolygonX[id1])]):
                                        NewWallFlowDensity[int(PolygonX[id1])][0]=max(NewWallFlowDensity[int(PolygonX[id1])])
                                        NewWallFlowDensity[int(PolygonX[id1])][1]=abs(float(WallFlowDensity[i][2][iLayer]))
                                i+=1
                            for i in range(len(ThickWallsX)):
                                myfile.write(str(float(mean(NewWallFlowDensity[i]))/sperd/cmperm) + " \n")  # min(sath2,  )
                        myfile.close()
                        text_file.close()
                    
                        text_file = open(newpath+"WallsThick3Dcos_bottomb"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "w")
                        with open(newpath+"WallsThick3Dcos_bottomb"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "a") as myfile:
                            myfile.write("# vtk DataFile Version 4.0 \n")
                            myfile.write("Wall geometry 3D including thickness bottom \n")
                            myfile.write("ASCII \n")
                            myfile.write(" \n")
                            myfile.write("DATASET UNSTRUCTURED_GRID \n")
                            myfile.write("POINTS "+str(len(ThickWallsX))+" float \n")
                            for ThickWallNodeX in ThickWallsX:
                                myfile.write(str(ThickWallNodeX[1]) + " " + str(ThickWallNodeX[2]) + " " + str(Cumheights[iLayer,0]/10) + " \n")
                            myfile.write(" \n")
                            myfile.write("CELLS " + str(int(NwallsJun+Nwalls-len(list_ghostwalls)*2-len(list_ghostjunctions))) + " " + str(int(2*Nwalls*5-len(list_ghostwalls)*10+sum(nWall2NewWallX[Nwalls:])+NwallsJun-Nwalls+2*len(Wall2NewWallX[Nwalls:])-nGhostJunction2Wall-len(list_ghostjunctions))) + " \n") #The number of cells corresponds to the number of lines in ThickWalls (if no ghost wall & junction)
                            i=0
                            for PolygonX in ThickWallPolygonX:
                                if floor(i/2) not in list_ghostwalls:
                                    myfile.write("4 " + str(int(PolygonX[0])) + " " + str(int(PolygonX[1])) + " " + str(int(PolygonX[2])) + " " + str(int(PolygonX[3])) + " \n")
                                i+=1
                            j=Nwalls
                            for PolygonX in Wall2NewWallX[Nwalls:]: #"junction" polygons
                                #Would need to order them based on x or y position to make sure display fully covers the surface (but here we try a simpler not so good solution instead)
                                if j not in list_ghostjunctions:
                                    string=str(int(nWall2NewWallX[j]+2)) #Added +2 so that the first and second nodes could be added again at the end (trying to fill the polygon better)
                                    for id1 in range(int(nWall2NewWallX[j])):
                                        string=string+" "+str(int(PolygonX[id1]))
                                    string=string+" "+str(int(PolygonX[0]))+" "+str(int(PolygonX[1])) #Adding the 1st and 2nd nodes again to the end
                                    myfile.write(string + " \n")
                                j+=1
                            myfile.write(" \n")
                            myfile.write("CELL_TYPES " + str(NwallsJun+Nwalls-len(list_ghostwalls)*2-len(list_ghostjunctions)) + " \n")
                            i=0
                            for PolygonX in ThickWallPolygonX:
                                if floor(i/2) not in list_ghostwalls:
                                    myfile.write("7 \n") #Polygon cell type (wall)
                                i+=1
                            j=Nwalls
                            for PolygonX in Wall2NewWallX[Nwalls:]:
                                if j not in list_ghostjunctions:
                                    myfile.write("6 \n") #Triangle-strip cell type (wall junction)
                                j+=1
                            myfile.write(" \n")
                            myfile.write("POINT_DATA " + str(len(ThickWallsX)) + " \n")
                            myfile.write("SCALARS Apo_flux_cosine_(m/s) float \n")
                            myfile.write("LOOKUP_TABLE default \n")
                            NewWallFlowDensity_cos=zeros((len(ThickWallsX),2))
                            i=0
                            for PolygonX in ThickWallPolygonX:
                                for id1 in range(4):
                                    if abs(float(WallFlowDensity_cos[i][2][iLayer]))>min(abs(NewWallFlowDensity_cos[int(PolygonX[id1])])):
                                        #Horizontal component of the flux
                                        if abs(NewWallFlowDensity_cos[int(PolygonX[id1])][1])>abs(NewWallFlowDensity_cos[int(PolygonX[id1])][0]): #Keeping the most extreme value
                                            NewWallFlowDensity_cos[int(PolygonX[id1])][0]=NewWallFlowDensity_cos[int(PolygonX[id1])][1]
                                        NewWallFlowDensity_cos[int(PolygonX[id1])][1]=float(WallFlowDensity_cos[i][2][iLayer])
                                i+=1
                            for i in range(len(ThickWallsX)):
                                myfile.write(str(float(mean(NewWallFlowDensity_cos[i]))/sperd/cmperm) + " \n")  # min(sath2,  )
                        myfile.close()
                        text_file.close()
                
                        if b>0:
                            text_file = open(newpath+"InterC3D_bottomb"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "w")
                            with open(newpath+"InterC3D_bottomb"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "a") as myfile:
                                myfile.write("# vtk DataFile Version 4.0 \n")
                                myfile.write("Intercellular space geometry 3D \n")
                                myfile.write("ASCII \n")
                                myfile.write(" \n")
                                myfile.write("DATASET UNSTRUCTURED_GRID \n")
                                myfile.write("POINTS "+str(len(ThickWalls))+" float \n")
                                for ThickWallNode in ThickWalls:
                                    myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(Cumheights[iLayer,0]/10) + " \n") #+height/200
                                myfile.write(" \n")
                                myfile.write("CELLS " + str(len(InterCid)) + " " + str(int(len(InterCid)+sum(nCell2ThickWalls[InterCid]))) + " \n") #The number of cells corresponds to the number of intercellular spaces
                                InterCFlowDensity=zeros((Ncells,1))
                                for cid in InterCid:
                                    n=int(nCell2ThickWalls[cid]) #Total number of thick wall nodes around the protoplast
                                    Polygon=Cell2ThickWalls[cid][:n]
                                    ranking=list()
                                    ranking.append(int(Polygon[0]))
                                    ranking.append(ThickWalls[int(ranking[0])][5])
                                    ranking.append(ThickWalls[int(ranking[0])][6])
                                    for id1 in range(1,n):
                                        wid1=ThickWalls[int(ranking[id1])][5]
                                        wid2=ThickWalls[int(ranking[id1])][6]
                                        if wid1 not in ranking:
                                            ranking.append(wid1)
                                        if wid2 not in ranking:
                                            ranking.append(wid2)
                                    string=str(n)
                                    for id1 in ranking:
                                        string=string+" "+str(int(id1))
                                    myfile.write(string + " \n")
                                    for twpid in Polygon[:int(n/2)]: #The first half of nodes are wall nodes actually connected to cells
                                        InterCFlowDensity[cid]+=abs(MembraneFlowDensity[int(twpid)][iLayer])/n #Mean absolute flow density calculation
                                myfile.write(" \n")
                                myfile.write("CELL_TYPES " + str(len(InterCid)) + " \n")
                                for i in range(len(InterCid)):
                                    myfile.write("6 \n") #Triangle-strip cell type
                                myfile.write(" \n")
                                myfile.write("POINT_DATA " + str(len(ThickWalls)) + " \n")
                                myfile.write("SCALARS Apo_flux_(m/s) float \n")
                                myfile.write("LOOKUP_TABLE default \n")
                                for ThickWallNode in ThickWalls:
                                    cellnumber1=ThickWallNode[2]-NwallsJun
                                    myfile.write(str(float(InterCFlowDensity[int(cellnumber1)])/sperd/cmperm) + " \n") #Flow rate from wall (non junction) to cell    min(sath1,max(satl1,  ))
                            myfile.close()
                            text_file.close()
                
                
                if ParaviewPF==1: #Plasmodesmata flow density data disks
                    
                    for Print_layer in Printrange:
                        iLayer=int(Print_layer.get("value"))
                        if PileUp==2:
                            b=int(BarriersL[iLayer])
                        else:
                            b=Barrier
                        text_file = open(newpath+"Plasmodesm3Db"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "w")
                        with open(newpath+"Plasmodesm3Db"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "a") as myfile:
                            myfile.write("# vtk DataFile Version 4.0 \n")
                            myfile.write("PD flux disks 3D \n")
                            myfile.write("ASCII \n")
                            myfile.write(" \n")
                            myfile.write("DATASET UNSTRUCTURED_GRID \n")
                            myfile.write("POINTS "+str(len(PlasmodesmFlowDensity)*12)+" float \n")
                            for ThickWallNode in ThickWalls:
                                if ThickWallNode[1]<Nwalls: #selection of new walls (not new junctions)
                                    if ThickWallNode[7]==0: #new walls that are not at the interface with soil or xylem, where there is no plasmodesmata   #if G.node[int(ThickWallNode[1])]['borderlink']==0
                                        #calculate the XY slope between the two neighbouring new junctions
                                        twpid1=int(ThickWallNode[5])
                                        twpid2=int(ThickWallNode[6])
                                        if not ThickWalls[twpid1][3]==ThickWalls[twpid2][3]: #Otherwise we'll get a division by 0 error
                                            slopeNJ=(ThickWalls[twpid1][4]-ThickWalls[twpid2][4])/(ThickWalls[twpid1][3]-ThickWalls[twpid2][3]) #slope of the line connecting the new junction nodes neighbouring the new wall
                                        else:
                                            slopeNJ=inf
                                        x0=ThickWallNode[3]
                                        y0=ThickWallNode[4]
                                        z0=radiusPlasmodesm_disp*3
                                        #Calculate the horizontal distance between XY0 and the cell center, compare it with the distance between the mean position of the new junctions. If the latter is closer to the cell center, it becomes the new XY0 to make sur the disk is visible
                                        xC=position[int(ThickWallNode[2])][0]
                                        yC=position[int(ThickWallNode[2])][1]
                                        xNJ=(ThickWalls[twpid1][3]+ThickWalls[twpid2][3])/2.0
                                        yNJ=(ThickWalls[twpid1][4]+ThickWalls[twpid2][4])/2.0
                                        if sqrt(square(x0-xC)+square(y0-yC)) > sqrt(square(xNJ-xC)+square(yNJ-yC)):
                                            x0=xNJ
                                            y0=yNJ
                                        for i in range(12):
                                            x=x0+cos(arctan(slopeNJ))*radiusPlasmodesm_disp*cos(int(i)*pi/6.0)
                                            y=y0+sin(arctan(slopeNJ))*radiusPlasmodesm_disp*cos(int(i)*pi/6.0)
                                            z=z0+radiusPlasmodesm_disp*sin(int(i)*pi/6.0)
                                            myfile.write(str(x) + " " + str(y) + " " + str(Cumheights[iLayer,0]/10+z) + " \n")
                                else:
                                    break #interrupts the for loop in case we reached the new junction nodes
                            myfile.write(" \n")
                            myfile.write("CELLS " + str(len(PlasmodesmFlowDensity)) + " " + str(len(PlasmodesmFlowDensity)*13) + " \n") #The number of cells corresponds to the number of lines in ThickWalls
                            for i in range(len(PlasmodesmFlowDensity)):
                                if PlasmodesmFlowDensity[i][iLayer]==0:
                                    myfile.write("12 " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " " + str(i*12+0) + " \n")
                                else:
                                    myfile.write("12 " + str(i*12+0) + " " + str(i*12+1) + " " + str(i*12+2) + " " + str(i*12+3) + " " + str(i*12+4) + " " + str(i*12+5) + " " + str(i*12+6) + " " + str(i*12+7) + " " + str(i*12+8) + " " + str(i*12+9) + " " + str(i*12+10) + " " + str(i*12+11) + " \n")
                            myfile.write(" \n")
                            myfile.write("CELL_TYPES " + str(len(PlasmodesmFlowDensity)) + " \n")
                            for i in range(len(PlasmodesmFlowDensity)):
                                myfile.write("7 \n") #Polygon cell type 
                            myfile.write(" \n")
                            myfile.write("POINT_DATA " + str(len(PlasmodesmFlowDensity)*12) + " \n")
                            myfile.write("SCALARS PD_Flux_(m/s) float \n")
                            myfile.write("LOOKUP_TABLE default \n")
                            for i in range(len(PlasmodesmFlowDensity)):
                                for j in range(12):
                                    myfile.write(str(float(PlasmodesmFlowDensity[i][iLayer])/sperd/cmperm) + " \n") #min(sath3,max(satl3, ))
                        myfile.close()
                        text_file.close()
                
                
                if ParaviewMF==1 and ParaviewPF==1: #Membranes and plasmodesms in the same file
                    for Print_layer in Printrange:
                        iLayer=int(Print_layer.get("value"))
                        if PileUp==2:
                            b=int(BarriersL[iLayer])
                        else:
                            b=Barrier
                        text_file = open(newpath+"Membranes_n_plasmodesm3Db"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "w")
                        with open(newpath+"Membranes_n_plasmodesm3Db"+str(b)+","+str(iMaturity)+"z"+str(iLayer)+"s"+str(count)+".pvtk", "a") as myfile:
                            myfile.write("# vtk DataFile Version 4.0 \n")
                            myfile.write("Membranes geometry and plasmodesm disks 3D \n")
                            myfile.write("ASCII \n")
                            myfile.write(" \n")
                            myfile.write("DATASET UNSTRUCTURED_GRID \n")
                            myfile.write("POINTS "+str(len(ThickWalls)*2+len(PlasmodesmFlowDensity)*12)+" float \n")
                            for ThickWallNode in ThickWalls:
                                myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(Cumheights[iLayer,0]/10) + " \n")
                            if PileUp==2 or AxialLayers>1:
                                for ThickWallNode in ThickWalls:
                                    if iLayer<AxialLayers-1:
                                        myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(Cumheights[iLayer+1,0]/10) + " \n")
                                    else:
                                        myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str((Cumheights[iLayer,0]+Cumheights[iLayer,0]-Cumheights[iLayer-1,0])/10) + " \n")
                            else:
                                myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(height/10) + " \n")
                            for ThickWallNode in ThickWalls:
                                if ThickWallNode[1]<Nwalls: #selection of new walls (not new junctions)
                                    if ThickWallNode[7]==0: #new walls that are not at the interface with soil or xylem, where there is no plasmodesmata   #if G.node[int(ThickWallNode[1])]['borderlink']==0
                                        #calculate the XY slope between the two neighbouring new junctions
                                        twpid1=int(ThickWallNode[5])
                                        twpid2=int(ThickWallNode[6])
                                        if not ThickWalls[twpid1][3]==ThickWalls[twpid2][3]: #Otherwise we'll get a division by 0 error
                                            slopeNJ=(ThickWalls[twpid1][4]-ThickWalls[twpid2][4])/(ThickWalls[twpid1][3]-ThickWalls[twpid2][3]) #slope of the line connecting the new junction nodes neighbouring the new wall
                                        else:
                                            slopeNJ=inf
                                        x0=ThickWallNode[3]
                                        y0=ThickWallNode[4]
                                        z0=radiusPlasmodesm_disp*3
                                        #Calculate the horizontal distance between XY0 and the cell center, compare it with the distance between the mean position of the new junctions. If the latter is closer to the cell center, it becomes the new XY0 to make sur the disk is visible
                                        xC=position[int(ThickWallNode[2])][0]
                                        yC=position[int(ThickWallNode[2])][1]
                                        xNJ=(ThickWalls[twpid1][3]+ThickWalls[twpid2][3])/2.0
                                        yNJ=(ThickWalls[twpid1][4]+ThickWalls[twpid2][4])/2.0
                                        if sqrt(square(x0-xC)+square(y0-yC)) > sqrt(square(xNJ-xC)+square(yNJ-yC)):
                                            x0=xNJ
                                            y0=yNJ
                                        for i in range(12):
                                            x=x0+cos(arctan(slopeNJ))*radiusPlasmodesm_disp*cos(int(i)*pi/6.0)
                                            y=y0+sin(arctan(slopeNJ))*radiusPlasmodesm_disp*cos(int(i)*pi/6.0)
                                            z=z0+radiusPlasmodesm_disp*sin(int(i)*pi/6.0)
                                            myfile.write(str(x) + " " + str(y) + " " + str(Cumheights[iLayer,0]/10+z) + " \n")
                                else:
                                    break #interrupts the for loop in case we reached the new junction nodes
                            myfile.write(" \n")
                            myfile.write("CELLS " + str(len(ThickWalls)-len(list_ghostwalls)*4+len(PlasmodesmFlowDensity)) + " " + str(len(ThickWalls)*5-len(list_ghostwalls)*20+len(PlasmodesmFlowDensity)*13) + " \n") #The number of cells corresponds to the number of lines in ThickWalls
                            for ThickWallNode in ThickWalls:
                                if ThickWallNode[1]>=Nwalls: #wall that is a junction
                                    if ThickWalls[int(ThickWallNode[5])][1] not in list_ghostwalls:
                                        myfile.write("4 " + str(int(ThickWallNode[0])) + " " + str(int(ThickWallNode[5])) + " " + str(int(ThickWallNode[5])+len(ThickWalls)) + " " + str(int(ThickWallNode[0])+len(ThickWalls)) + " \n")
                                    if ThickWalls[int(ThickWallNode[6])][1] not in list_ghostwalls:
                                        myfile.write("4 " + str(int(ThickWallNode[0])) + " " + str(int(ThickWallNode[6])) + " " + str(int(ThickWallNode[6])+len(ThickWalls)) + " " + str(int(ThickWallNode[0])+len(ThickWalls)) + " \n")
                            for i in range(len(PlasmodesmFlowDensity)):
                                if PlasmodesmFlowDensity[i][iLayer]==0:
                                    myfile.write("12 " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+0+len(ThickWalls)*2) + " \n")
                                else:
                                    myfile.write("12 " + str(i*12+0+len(ThickWalls)*2) + " " + str(i*12+1+len(ThickWalls)*2) + " " + str(i*12+2+len(ThickWalls)*2) + " " + str(i*12+3+len(ThickWalls)*2) + " " + str(i*12+4+len(ThickWalls)*2) + " " + str(i*12+5+len(ThickWalls)*2) + " " + str(i*12+6+len(ThickWalls)*2) + " " + str(i*12+7+len(ThickWalls)*2) + " " + str(i*12+8+len(ThickWalls)*2) + " " + str(i*12+9+len(ThickWalls)*2) + " " + str(i*12+10+len(ThickWalls)*2) + " " + str(i*12+11+len(ThickWalls)*2) + " \n")
                            myfile.write(" \n")
                            myfile.write("CELL_TYPES " + str(len(ThickWalls)-len(list_ghostwalls)*4+len(PlasmodesmFlowDensity)) + " \n")
                            for ThickWallNode in ThickWalls:
                                if ThickWallNode[1]>=Nwalls: #wall that is a junction
                                    if ThickWalls[int(ThickWallNode[5])][1] not in list_ghostwalls:
                                        myfile.write("9 \n") #Quad cell type
                                    if ThickWalls[int(ThickWallNode[6])][1] not in list_ghostwalls:
                                        myfile.write("9 \n") #Quad cell type
                            for i in range(len(PlasmodesmFlowDensity)):
                                myfile.write("7 \n") #Polygon cell type 
                            myfile.write(" \n")
                            myfile.write("POINT_DATA " + str(len(ThickWalls)*2+len(PlasmodesmFlowDensity)*12) + " \n")
                            myfile.write("SCALARS TM_n_PD_flux_(m/s) float \n")
                            myfile.write("LOOKUP_TABLE default \n")
                            for ThickWallNode in ThickWalls:
                                if ThickWallNode[0]<len(MembraneFlowDensity):
                                    myfile.write(str(float(MembraneFlowDensity[int(ThickWallNode[0])][iLayer])/sperd/cmperm) + " \n") #Flow rate from wall (non junction) to cell
                                else:
                                    myfile.write(str(float((MembraneFlowDensity[int(ThickWallNode[5])][iLayer]+MembraneFlowDensity[int(ThickWallNode[6])][iLayer])/2)/sperd/cmperm) + " \n") #Flow rate from junction wall to cell is the average of the 2 neighbouring wall flow rates
                            for ThickWallNode in ThickWalls:
                                if ThickWallNode[0]<len(MembraneFlowDensity):
                                    myfile.write(str(float(MembraneFlowDensity[int(ThickWallNode[0])][iLayer])/sperd/cmperm) + " \n") #Flow rate from wall (non junction) to cell
                                else:
                                    myfile.write(str(float((MembraneFlowDensity[int(ThickWallNode[5])][iLayer]+MembraneFlowDensity[int(ThickWallNode[6])][iLayer])/2)/sperd/cmperm) + " \n") #Flow rate from junction wall to cell is the average of the 2 neighbouring wall flow rates
                            for i in range(len(PlasmodesmFlowDensity)):
                                for j in range(12):
                                    myfile.write(str(float(PlasmodesmFlowDensity[i][iLayer])/sperd/cmperm) + " \n")
                        myfile.close()
                        text_file.close()
                    
        #Maturity loop (break if PileUp==2?)
        if ParaviewPF==1 and PD_freq_source==2: # or Sym_Contagion==2
            text_file = open(newpath+"PD_frequency_bottom.pvtk", "w")
            with open(newpath+"PD_frequency_bottom.pvtk", "a") as myfile:
                myfile.write("# vtk DataFile Version 4.0 \n")
                myfile.write("Plasmodesmata frequency per membrane surface area \n")
                myfile.write("ASCII \n")
                myfile.write(" \n")
                myfile.write("DATASET UNSTRUCTURED_GRID \n")
                if UniXwalls==1:
                    myfile.write("POINTS "+str(len(ThickWallsX)+len(ThickWallPolygonX)*4)+" float \n")
                    for ThickWallNodeX in ThickWallsX:
                        myfile.write(str(ThickWallNodeX[1]) + " " + str(ThickWallNodeX[2]) + " " + str(Cumheights[iLayer,0]/10) + " \n")
                    for PolygonX in ThickWallPolygonX:
                        for NodeX in PolygonX:
                            myfile.write(str(ThickWallsX[int(NodeX)][1]) + " " + str(ThickWallsX[int(NodeX)][2]) + " " + str(Cumheights[iLayer,0]/10) + " \n")                        
                    myfile.write(" \n")
                    myfile.write("CELLS " + str(int(NwallsJun+Nwalls-len(list_ghostwalls)*2-len(list_ghostjunctions))) + " " + str(int(2*Nwalls*5-len(list_ghostwalls)*10+sum(nWall2NewWallX[Nwalls:])+NwallsJun-Nwalls+2*len(Wall2NewWallX[Nwalls:])-nGhostJunction2Wall-len(list_ghostjunctions))) + " \n") #The number of cells corresponds to the number of lines in ThickWalls (if no ghost wall & junction)
                    i=len(ThickWallsX)
                    for PolygonX in ThickWallPolygonX:
                        if floor((i-len(ThickWallsX))/4) not in list_ghostwalls:
                            myfile.write("4 " + str(i) + " " + str(i+1) + " " + str(i+2) + " " + str(i+3) + " \n")
                        i+=4
                    j=Nwalls
                    for PolygonX in Wall2NewWallX[Nwalls:]: #"junction" polygons
                        #Would need to order them based on x or y position to make sure display fully covers the surface (but here we try a simpler not so good solution instead)
                        if j not in list_ghostjunctions:
                            string=str(int(nWall2NewWallX[j]+2)) #Added +2 so that the first and second nodes could be added again at the end (trying to fill the polygon better)
                            for id1 in range(int(nWall2NewWallX[j])):
                                string=string+" "+str(int(PolygonX[id1]))
                            string=string+" "+str(int(PolygonX[0]))+" "+str(int(PolygonX[1])) #Adding the 1st and 2nd nodes again to the end
                            myfile.write(string + " \n")
                        j+=1
                    myfile.write(" \n")
                    myfile.write("CELL_TYPES " + str(NwallsJun+Nwalls-len(list_ghostwalls)*2-len(list_ghostjunctions)) + " \n")
                    i=0
                    for PolygonX in ThickWallPolygonX:
                        if floor(i/2) not in list_ghostwalls:
                            myfile.write("7 \n") #Polygon cell type (wall)
                        i+=1
                    j=Nwalls
                    for PolygonX in Wall2NewWallX[Nwalls:]:
                        if j not in list_ghostjunctions:
                            myfile.write("6 \n") #Triangle-strip cell type (wall junction)
                        j+=1
                    myfile.write(" \n")
                    myfile.write("POINT_DATA " + str(len(ThickWallsX)+len(ThickWallPolygonX)*4) + " \n")
                    myfile.write("SCALARS PD_frequency(#_per_square_micron) float \n")
                    myfile.write("LOOKUP_TABLE default \n")
                    NewWall_PD=zeros((len(ThickWallsX),1))
                    #j=0
                    #for PolygonX in Wall2NewWallX:
                    #    for id1 in range(int(nWall2NewWallX[j])):
                    #        if j<Nwalls:
                    #            NewWall_PD[int(PolygonX[id1])]=float(Wall_PD[j])
                    #        else:
                    #            NewWall_PD[int(PolygonX[id1])]=mean(Wall_PD[Junction2Wall[j-Nwalls]])
                    #    j+=1
                    wid=0
                    for PolygonX in ThickWallPolygonX:
                        for id1 in range(4):
                            NewWall_PD[int(PolygonX[id1])]+=float(Wall_PD[int(floor(wid/2))])/2
                        wid+=1
                    for i in range(len(ThickWallsX)):
                        myfile.write(str(float(NewWall_PD[i])) + " \n")
                    for it1 in range(0,Nwalls):
                        for it2 in range(0,8):
                            myfile.write(str(float(Wall_PD[it1])) + " \n") 
                else:
                    myfile.write("POINTS "+str(len(ThickWallsX))+" float \n")
                    for ThickWallNodeX in ThickWallsX:
                        myfile.write(str(ThickWallNodeX[1]) + " " + str(ThickWallNodeX[2]) + " " + str(Cumheights[iLayer,0]/10) + " \n")
                    myfile.write(" \n")
                    myfile.write("CELLS " + str(int(NwallsJun+Nwalls-len(list_ghostwalls)*2-len(list_ghostjunctions))) + " " + str(int(2*Nwalls*5-len(list_ghostwalls)*10+sum(nWall2NewWallX[Nwalls:])+NwallsJun-Nwalls+2*len(Wall2NewWallX[Nwalls:])-nGhostJunction2Wall-len(list_ghostjunctions))) + " \n") #The number of cells corresponds to the number of lines in ThickWalls (if no ghost wall & junction)
                    i=0
                    for PolygonX in ThickWallPolygonX:
                        if floor(i/2) not in list_ghostwalls:
                            myfile.write("4 " + str(int(PolygonX[0])) + " " + str(int(PolygonX[1])) + " " + str(int(PolygonX[2])) + " " + str(int(PolygonX[3])) + " \n")
                        i+=1
                    j=Nwalls
                    for PolygonX in Wall2NewWallX[Nwalls:]: #"junction" polygons
                        #Would need to order them based on x or y position to make sure display fully covers the surface (but here we try a simpler not so good solution instead)
                        if j not in list_ghostjunctions:
                            string=str(int(nWall2NewWallX[j]+2)) #Added +2 so that the first and second nodes could be added again at the end (trying to fill the polygon better)
                            for id1 in range(int(nWall2NewWallX[j])):
                                string=string+" "+str(int(PolygonX[id1]))
                            string=string+" "+str(int(PolygonX[0]))+" "+str(int(PolygonX[1])) #Adding the 1st and 2nd nodes again to the end
                            myfile.write(string + " \n")
                        j+=1
                    myfile.write(" \n")
                    myfile.write("CELL_TYPES " + str(NwallsJun+Nwalls-len(list_ghostwalls)*2-len(list_ghostjunctions)) + " \n")
                    i=0
                    for PolygonX in ThickWallPolygonX:
                        if floor(i/2) not in list_ghostwalls:
                            myfile.write("7 \n") #Polygon cell type (wall)
                        i+=1
                    j=Nwalls
                    for PolygonX in Wall2NewWallX[Nwalls:]:
                        if j not in list_ghostjunctions:
                            myfile.write("6 \n") #Triangle-strip cell type (wall junction)
                        j+=1
                    myfile.write(" \n")
                    myfile.write("POINT_DATA " + str(len(ThickWallsX)) + " \n")
                    myfile.write("SCALARS PD_frequency(#_per_square_micron) float \n")
                    myfile.write("LOOKUP_TABLE default \n")
                    NewWall_PD=zeros((len(ThickWallsX),1))
                    #j=0
                    #for PolygonX in Wall2NewWallX:
                    #    for id1 in range(int(nWall2NewWallX[j])):
                    #        if j<Nwalls:
                    #            NewWall_PD[int(PolygonX[id1])]=float(Wall_PD[j])
                    #        else:
                    #            NewWall_PD[int(PolygonX[id1])]=mean(Wall_PD[Junction2Wall[j-Nwalls]])
                    #    j+=1
                    wid=0
                    for PolygonX in ThickWallPolygonX:
                        for id1 in range(4):
                            NewWall_PD[int(PolygonX[id1])]+=float(Wall_PD[int(floor(wid/2))])/2
                        wid+=1
                    for i in range(len(ThickWallsX)):
                        myfile.write(str(float(NewWall_PD[i])) + " \n")
            myfile.close()
            text_file.close()
    
    #write down kr_tot and Uptake distributions in matrices
    #iMaturity=-1
    TopLayer=0
    for iMaturity in range(Nmaturity_loop):
        if PileUp==2:
            b=9
            AxialLayers=TotLayers #Total number of stacked cross-sections across maturity levels
        else:
            b=int(Maturityrange[iMaturity].get("Barrier"))
            AxialLayers=int(Maturityrange[iMaturity].get("Nlayers")) #Total number of stacked cross-sections at this maturity level
        #iMaturity+=1
        TopLayer+=AxialLayers
        text_file = open(newpath+"Macro_prop_"+str(b)+","+str(iMaturity)+".txt", "w")
        with open(newpath+"Macro_prop_"+str(b)+","+str(iMaturity)+".txt", "a") as myfile:
            myfile.write("Macroscopic root radial hydraulic properties, apoplastic barrier "+str(b)+","+str(iMaturity)+" \n")
            myfile.write("\n")
            myfile.write(str(Nscenarios-1)+" scenarios \n")
            myfile.write("\n")
            myfile.write("Stack height: "+str((Totheight)*1.0E-04)+" cm \n")
            myfile.write("\n")
            myfile.write("Number of zones: "+str(len(Nlayers))+" \n")
            myfile.write("\n")
            temp=str(Nlayers)
            #if len(Nlayers)>1:
            myfile.write("Number of layers: "+temp[1:-1]+" \n")
            #else:
            #    myfile.write("Number of layers: "+temp+" \n")
            myfile.write("\n")
            if PileUp==2:
                temp=str(Barriers)
                myfile.write("Type of layers: "+temp[1:-1]+" \n")
            else:
                myfile.write("Type of layers: "+str(b)+" \n")
            myfile.write("\n")
            myfile.write("Cross-section perimeter: "+str(perimeter[0])+" cm \n")
            myfile.write("\n")
            myfile.write("Xylem specific axial conductance: "+str(K_xyl_spec)+" cm^4/hPa/d \n")
            myfile.write("\n")
            myfile.write("Cross-section radial conductivity: "+str(kr_tot[iMaturity][0])+" cm/hPa/d \n")
            myfile.write("\n")
            myfile.write("Number of radial discretization boxes: \n")
            r_discret_txt=' '.join(map(str, r_discret.T)) 
            myfile.write(r_discret_txt[1:21]+" \n")
            myfile.write("\n")
            myfile.write("Radial distance from stele centre (microns): \n")
            for j in Layer_dist2:
                myfile.write(str(float(j))+" \n")
            myfile.write("\n")
            myfile.write("Standard Transmembrane uptake Fractions (%): \n")
            for j in range(int(r_discret[0])):
                if AxialLayers==1:
                    myfile.write(str(float(STFlayer_plus[j,TopLayer-AxialLayers:TopLayer]*100))+" \n")
                else:
                    temp=str(list(STFlayer_plus[j,TopLayer-AxialLayers:TopLayer]*100))
                    myfile.write(temp[1:-1]+" \n")
            myfile.write("\n")
            myfile.write("Standard Transmembrane release Fractions (%): \n")
            for j in range(int(r_discret[0])):
                if AxialLayers==1:
                    myfile.write(str(float(STFlayer_minus[j,TopLayer-AxialLayers:TopLayer]*100))+" \n")
                else:
                    temp=str(list(STFlayer_minus[j,TopLayer-AxialLayers:TopLayer]*100))
                    myfile.write(temp[1:-1]+" \n")
            for i in range(1,Nscenarios):
                myfile.write("\n")
                myfile.write("\n")
                myfile.write("Scenario "+str(i)+" \n")
                myfile.write("\n")
                myfile.write("h_x: "+str(Psi_xyl[1][iMaturity][i])+" hPa \n")
                myfile.write("\n")
                myfile.write("h_s: "+str(Psi_soil[0][i])+" to "+str(Psi_soil[1][i])+" hPa \n")
                myfile.write("\n")
                myfile.write("h_p: "+str(Psi_sieve[1][iMaturity][i])+" hPa \n")
                myfile.write("\n")
                #if AxialLayers>1:
                #    temp=str(list(Os_apo_eq[:,1,i]))
                #    myfile.write("O_apo_stele_eq: "+temp[1:-1]+" hPa \n")
                #    temp=str(list(Os_sym_eq[:,1,i]))
                #    myfile.write("O_sym_stele_eq: "+temp[1:-1]+" hPa \n")
                #    temp=str(list(Os_apo_eq[:,0,i]))
                #    myfile.write("O_apo_cortex_eq: "+temp[1:-1]+" hPa \n")
                #    temp=str(list(Os_sym_eq[:,0,i]))
                #    myfile.write("O_sym_cortex_eq: "+temp[1:-1]+" hPa \n")
                #else:
                myfile.write("O_apo_stele_eq: "+str(Os_apo_eq[iMaturity,1,i])+" hPa \n")
                myfile.write("O_sym_stele_eq: "+str(Os_sym_eq[iMaturity,1,i])+" hPa \n")
                myfile.write("O_apo_cortex_eq: "+str(Os_apo_eq[iMaturity,0,i])+" hPa \n")
                myfile.write("O_sym_cortex_eq: "+str(Os_sym_eq[iMaturity,0,i])+" hPa \n")
                myfile.write("\n")
                myfile.write("O_x: "+str(Os_xyl[0][i])+" to "+str(Os_xyl[1][i])+" hPa \n")
                myfile.write("\n")
                myfile.write("O_s: "+str(Os_soil[0][i])+" to "+str(Os_soil[1][i])+" hPa \n")
                myfile.write("\n")
                myfile.write("O_p: "+str(Os_sieve[0][i])+" hPa \n")
                myfile.write("\n")
                if PileUp==2:
                    myfile.write("Xcontact: "+str(Xcontacts)+" microns \n")
                else:
                    myfile.write("Xcontact: "+str(Xcontact)+" microns \n")
                myfile.write("\n")
                if b==0:
                    myfile.write("Elong_cell: "+str(Elong_cell[0][i])+" cm/d \n")
                    myfile.write("\n")
                    myfile.write("Elong_cell_side_diff: "+str(Elong_cell_side_diff[0][i])+" cm/d \n")
                    myfile.write("\n")
                else:
                    myfile.write("Elong_cell: "+str(0.0)+" cm/d \n")
                    myfile.write("\n")
                    myfile.write("Elong_cell_side_diff: "+str(0.0)+" cm/d \n")
                    myfile.write("\n")
                myfile.write("kw: "+str(kw)+" cm^2/hPa/d \n")
                myfile.write("\n")
                myfile.write("Kpl: "+str(Kpl)+" cm^3/hPa/d \n")
                myfile.write("\n")
                myfile.write("kAQP: "+str(kaqp_cortex)+" cm/hPa/d \n")
                myfile.write("\n")
                myfile.write("s_hetero: "+str(s_hetero[0][count])+" \n")
                myfile.write("\n")
                myfile.write("s_factor: "+str(s_factor[0][count])+" \n")
                myfile.write("\n")
                myfile.write("Os_hetero: "+str(Os_hetero[0][count])+" \n")
                myfile.write("\n")
                myfile.write("Os_cortex: "+str(Os_cortex[0][count])+" hPa \n")
                myfile.write("\n")
                myfile.write("q_tot: "+str(Q_tot[iMaturity][i]/(Totheight)/1.0E-04)+" cm^2/d \n")
                myfile.write("\n")
                myfile.write("Stele, cortex, and epidermis uptake distribution cm^3/d: \n")
                for j in range(int(r_discret[0])):
                    if AxialLayers==1:
                        myfile.write(str(float(UptakeLayer_plus[j,TopLayer-AxialLayers:TopLayer,i]))+" \n")
                    else:
                        temp=str(list(UptakeLayer_plus[j,TopLayer-AxialLayers:TopLayer,i]))
                        myfile.write(temp[1:-1]+" \n")
                myfile.write("\n")
                myfile.write("Stele, cortex, and epidermis release distribution cm^3/d: \n")
                for j in range(int(r_discret[0])):
                    if AxialLayers==1:
                        myfile.write(str(float(UptakeLayer_minus[j,TopLayer-AxialLayers:TopLayer,i]))+" \n")
                    else:
                        temp=str(list(UptakeLayer_minus[j,TopLayer-AxialLayers:TopLayer,i]))
                        myfile.write(temp[1:-1]+" \n")
                myfile.write("\n")
                myfile.write("Xylem uptake distribution cm^3/d: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(float(Q_xyl_layer[j][iMaturity][i]))+" \n")
                myfile.write("\n")
                myfile.write("Phloem uptake distribution cm^3/d: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(float(Q_sieve_layer[j][iMaturity][i]))+" \n")
                myfile.write("\n")
                myfile.write("Elongation flow convergence distribution cm^3/d: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(float(Q_elong_layer[j][iMaturity][i]))+" \n")
                myfile.write("\n")
                myfile.write("Cell layers pressure potentials: \n")
                for j in range(int(r_discret[0])):
                    if AxialLayers==1:
                        myfile.write(str(float(PsiCellLayer[j,TopLayer-AxialLayers:TopLayer,i]))+" \n")
                    else:
                        temp=str(list(PsiCellLayer[j,TopLayer-AxialLayers:TopLayer,i]))
                        myfile.write(temp[1:-1]+" \n")
                myfile.write("\n")
                myfile.write("Cell layers osmotic potentials: \n")
                if PileUp==2:
                    for j in range(int(r_discret[0])):
                        temp=str(list(OsCellLayer[j,:,i]))
                        myfile.write(temp[1:-1]+" \n")
                else:
                    for j in range(int(r_discret[0])):
                        myfile.write(str(float(OsCellLayer[j][iMaturity][i]))+" \n")
                myfile.write("\n")
                myfile.write("Wall layers pressure potentials: \n")
                for j in range(int(r_discret[0])):
                    if NWallLayer[j][iMaturity][i]>1:
                        if AxialLayers>1:
                            temp=str(list(PsiWallLayer[j,TopLayer-AxialLayers:TopLayer,i]/NWallLayer[j,iMaturity,i]))
                            myfile.write(temp[1:-1]+" \n")
                        else:
                            myfile.write(str(float(PsiWallLayer[j,TopLayer-AxialLayers:TopLayer,i]/NWallLayer[j,iMaturity,i]))+" \n")
                    else:
                        if AxialLayers>1:
                            temp=str(list(PsiWallLayer[j,TopLayer-AxialLayers:TopLayer,i]))
                            myfile.write(temp[1:-1]+" \n")
                        else:
                            myfile.write(str(float(PsiWallLayer[j,TopLayer-AxialLayers:TopLayer,i]))+" \n")
                myfile.write("\n")
                myfile.write("Wall layers osmotic potentials: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(float(OsWallLayer[j][iMaturity][i]))+" \n")
        myfile.close()
        text_file.close()
    
    if Sym_Contagion == 1: #write down results of the hydropatterning study
        #iMaturity=-1
        for iMaturity in range(Nmaturity_loop):
            if PileUp==2:
                b=9
            else:
                b=int(Maturityrange[iMaturity].get("Barrier"))
            #iMaturity+=1
            text_file = open(newpath+"Hydropatterning_"+str(b)+","+str(iMaturity)+".txt", "w")
            with open(newpath+"Hydropatterning_"+str(b)+","+str(iMaturity)+".txt", "a") as myfile:
                myfile.write("Is there symplastic mass flow from source to target cells? Apoplastic barrier "+str(b)+","+str(iMaturity)+" \n")
                myfile.write("\n")
                myfile.write(str(Nscenarios-1)+" scenarios \n")
                myfile.write("\n")
                myfile.write("Template: "+path_geom+" \n")
                myfile.write("\n")
                myfile.write("Source cell: "+str(N_Dirichlet_ini)+" \n")
                myfile.write("\n")
                myfile.write("Target cells: "+str(Sym_Target)+" \n")
                myfile.write("\n")
                myfile.write("Immune cells: "+str(Sym_Immune)+" \n")
                myfile.write("\n")
                myfile.write("Stack height: "+str((Totheight)*1.0E-04)+" cm \n")
                myfile.write("\n")
                myfile.write("Cross-section perimeter: "+str(perimeter[0])+" cm \n")
                myfile.write("\n")
                myfile.write("Xcontact: "+str(Xcontact)+" microns \n")
                myfile.write("\n")
                myfile.write("kw: "+str(kw)+" cm^2/hPa/d \n")
                myfile.write("\n")
                myfile.write("Kpl: "+str(Kpl)+" cm^3/hPa/d \n")
                myfile.write("\n")
                myfile.write("kAQP: "+str(kaqp_cortex)+" cm/hPa/d \n")
                myfile.write("\n")
                if b==0:
                    myfile.write("Cell elongation rate: "+str(Elong_cell)+" cm/d \n")
                else: #No elongation after formation of the Casparian strip
                    myfile.write("Cell elongation rate: "+str(0.0)+" cm/d \n")
                myfile.write("\n")
                for i in range(1,Nscenarios):
                    myfile.write("\n")
                    myfile.write("\n")
                    myfile.write("Scenario "+str(i)+" \n")
                    myfile.write("\n")
                    myfile.write("Expected hydropatterining response (1: Wet-side XPP; -1 to 0: Unclear; 2: Dry-side XPP) \n")
                    myfile.write("Hydropat.: "+str(int(Hydropatterning[iMaturity][i]))+" \n")
                    myfile.write("\n")
                    myfile.write("h_x: "+str(Psi_xyl[1][iMaturity][i])+" hPa, h_s: "+str(Psi_soil[0][i])+" to "+str(Psi_soil[1][i])+" hPa, h_p: "+str(Psi_sieve[1][iMaturity][i])+" hPa \n")
                    myfile.write("\n")
                    myfile.write("O_x: "+str(Os_xyl[0][i])+" to "+str(Os_xyl[0][i])+" hPa, O_s: "+str(Os_soil[0][i])+" to "+str(Os_soil[1][i])+" hPa, O_p: "+str(Os_sieve[0][i])+" hPa \n")
                    myfile.write("\n")
                    myfile.write("Os_cortex: "+str(Os_cortex[0][count])+" hPa, Os_hetero: "+str(Os_hetero[0][count])+", s_hetero: "+str(s_hetero[0][count])+", s_factor: "+str(s_factor[0][count])+" \n")
                    myfile.write("\n")
                    myfile.write("q_tot: "+str(Q_tot[iMaturity][i]/(Totheight)/1.0E-04)+" cm^2/d \n")
                    myfile.write("\n")
            myfile.close()
            text_file.close()
    
    if Apo_Contagion == 1: #write down results of the hydrotropism study
        #iMaturity=-1
        for iMaturity in range(Nmaturity_loop):
            b=int(Maturityrange[iMaturity].get("Barrier"))
            #iMaturity+=1
            text_file = open(newpath+"Hydrotropism_"+str(b)+","+str(iMaturity)+".txt", "w")
            with open(newpath+"Hydrotropism_"+str(b)+","+str(iMaturity)+".txt", "a") as myfile:
                myfile.write("Is there apoplastic mass flow from source to target cells? Apoplastic barrier "+str(b)+","+str(iMaturity)+" \n")
                myfile.write("\n")
                myfile.write(str(Nscenarios-1)+" scenarios \n")
                myfile.write("\n")
                myfile.write("Template: "+path_geom+" \n")
                myfile.write("\n")
                myfile.write("Source cell: "+str(N_Dirichlet_ini)+" \n")
                myfile.write("\n")
                myfile.write("Target cells: "+str(Apo_Target)+" \n")
                myfile.write("\n")
                myfile.write("Immune cells: "+str(Apo_Immune)+" \n")
                myfile.write("\n")
                myfile.write("Stack height: "+str((Totheight)*1.0E-04)+" cm \n")
                myfile.write("\n")
                myfile.write("Cross-section perimeter: "+str(perimeter[0])+" cm \n")
                myfile.write("\n")
                myfile.write("Xcontact: "+str(Xcontact)+" microns \n")
                myfile.write("\n")
                myfile.write("kw: "+str(kw)+" cm^2/hPa/d \n")
                myfile.write("\n")
                myfile.write("Kpl: "+str(Kpl)+" cm^3/hPa/d \n")
                myfile.write("\n")
                myfile.write("kAQP: "+str(kaqp_cortex)+" cm/hPa/d \n")
                myfile.write("\n")
                if b==0:
                    myfile.write("Cell elongation rate: "+str(Elong_cell)+" cm/d \n")
                else: #No elongation after formation of the Casparian strip
                    myfile.write("Cell elongation rate: "+str(0.0)+" cm/d \n")
                myfile.write("\n")
                for i in range(1,Nscenarios):
                    myfile.write("\n")
                    myfile.write("\n")
                    myfile.write("Scenario "+str(i)+" \n")
                    myfile.write("\n")
                    myfile.write("Expected hydrotropism response (1: All cell walls reached by ABA; 0: No target walls reached by ABA) \n")
                    myfile.write("Hydropat.: "+str(int(Hydrotropism[iMaturity][i]))+" \n")
                    myfile.write("\n")
                    myfile.write("h_x: "+str(Psi_xyl[1][iMaturity][i])+" hPa, h_s: "+str(Psi_soil[0][i])+" to "+str(Psi_soil[1][i])+" hPa, h_p: "+str(Psi_sieve[1][iMaturity][i])+" hPa \n")
                    myfile.write("\n")
                    myfile.write("O_x: "+str(Os_xyl[0][i])+" to "+str(Os_xyl[0][i])+" hPa, O_s: "+str(Os_soil[0][i])+" to "+str(Os_soil[1][i])+" hPa, O_p: "+str(Os_sieve[0][i])+" hPa \n")
                    myfile.write("\n")
                    myfile.write("Os_cortex: "+str(Os_cortex[0][count])+" hPa, Os_hetero: "+str(Os_hetero[0][count])+", s_hetero: "+str(s_hetero[0][count])+", s_factor: "+str(s_factor[0][count])+" \n")
                    myfile.write("\n")
                    myfile.write("q_tot: "+str(Q_tot[iMaturity][i]/(Totheight)/1.0E-04)+" cm^2/d \n")
                    myfile.write("\n")
            myfile.close()
            text_file.close()
        
    if InterC_perim_search==1:
        text_file = open(newpath+"Cortical_cell_perimeters.txt", "w")
        with open(newpath+"Cortical_cell_perimeters.txt", "a") as myfile:
            for j in range(10):
                myfile.write(str(cortex_cellperimeters[j][:])+" \n")
        myfile.close()
        text_file.close()
    
    

if Paraview==1:
    rank1=28
    rank2=28
    NselecC=sum(logical_and(Cell_rank>=rank1,Cell_rank<=rank2))
    text_file = open(newpath+"Cell_rank_"+str(rank1)+"_to_"+str(rank2)+"_bottom.pvtk", "w")
    with open(newpath+"Cell_rank_"+str(rank1)+"_to_"+str(rank2)+"_bottom.pvtk", "a") as myfile:
        myfile.write("# vtk DataFile Version 4.0 \n")
        myfile.write("Intercellular space geometry 3D \n")
        myfile.write("ASCII \n")
        myfile.write(" \n")
        myfile.write("DATASET UNSTRUCTURED_GRID \n")
        myfile.write("POINTS "+str(len(ThickWalls))+" float \n")
        for ThickWallNode in ThickWalls:
            myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(Cumheights[iLayer,0]/10) + " \n") #+height/200
        myfile.write(" \n")
        myfile.write("CELLS " + str(NselecC) + " " + str(int(NselecC+sum(nCell2ThickWalls[logical_and(Cell_rank>=rank1,Cell_rank<=rank2)]))) + " \n") #The number of cells corresponds to the number of intercellular spaces
        selecCFlowDensity=zeros((Ncells,1))
        for cid in range(Ncells):
            if logical_and(Cell_rank[cid]>=rank1,Cell_rank[cid]<=rank2):
                n=int(nCell2ThickWalls[cid]) #Total number of thick wall nodes around the protoplast
                Polygon=Cell2ThickWalls[cid][:n]
                ranking=list()
                ranking.append(int(Polygon[0]))
                ranking.append(ThickWalls[int(ranking[0])][5])
                ranking.append(ThickWalls[int(ranking[0])][6])
                for id1 in range(1,n):
                    wid1=ThickWalls[int(ranking[id1])][5]
                    wid2=ThickWalls[int(ranking[id1])][6]
                    if wid1 not in ranking:
                        ranking.append(wid1)
                    if wid2 not in ranking:
                        ranking.append(wid2)
                string=str(n)
                for id1 in ranking:
                    string=string+" "+str(int(id1))
                myfile.write(string + " \n")
                for twpid in Polygon[:int(n/2)]: #The first half of nodes are wall nodes actually connected to cells
                    selecCFlowDensity[cid]+=abs(MembraneFlowDensity[int(twpid)][iLayer])/n #Mean absolute flow density calculation
        myfile.write(" \n")
        myfile.write("CELL_TYPES " + str(NselecC) + " \n")
        for i in range(NselecC):
            myfile.write("6 \n") #Triangle-strip cell type
        myfile.write(" \n")
        myfile.write("POINT_DATA " + str(len(ThickWalls)) + " \n")
        myfile.write("SCALARS Apo_flux_(m/s) float \n")
        myfile.write("LOOKUP_TABLE default \n")
        for ThickWallNode in ThickWalls:
            cellnumber1=ThickWallNode[2]-NwallsJun
            myfile.write(str(float(selecCFlowDensity[int(cellnumber1)])/sperd/cmperm) + " \n") #Flow rate from wall (non junction) to cell    min(sath1,max(satl1,  ))
    myfile.close()
    text_file.close()
    
    
    #text_file = open(newpath+"PD_flow_rates.txt", "w")
    #with open(newpath+"Cortical_cell_perimeters.txt", "a") as myfile:
    #    for j in range(10):
    #        myfile.write(str(cortex_cellperimeters[j][:])+" \n")
    #myfile.close()
    #text_file.close()
t1 = time.perf_counter()
print(t1-t0, "seconds process time")