# -*- coding: utf-8 -*-

#Directory
dir='C:/Users/heymansad/Documents/GitHub/granar_examples/'

#Project
Project='Projects/granar/'#BBSRC/'#'Projects/

#Inputs
Gen='Maize_General.xml'#'Arabido1_General.xml' #'MilletLR3_General.xml' #
Geom='Root_example.xml'#'Arabido4_Geometry_BBSRC.xml' #'Maize2_Geometry.xml' #''MilletLR3_Geometry.xml' #'Wheat1_Nodal_Geometry_aerenchyma.xml' #'Maize1_Geometry.xml' #
Hydr='Maize_Hydraulics.xml' #'Arabido1_Hydraulics_ERC.xml' #'MilletLR3_Hydraulics.xml' #'Test_Hydraulics.xml' #
BC='Maize_BC_kr.xml' #'Arabido4_BC_BBSRC2.xml' #'Arabido1_BC_Emily.xml' #'Arabido3_BC_BBSRC.xml' #'Maize_BC_SoluteAna_krOsmo.xml'#'Maize_BC_OSxyl_hetero.xml' #'Arabido1_BC_Emily.xml' #'BC_Test.xml' #'Maize_BC_Plant_phys.xml'
Horm='Maize_Hormones_Carriers.xml'

#Libraries
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

import scipy.linalg as slin #Linear algebra functions

import pylab #Found in the package pyqt
from pylab import *  # for plotting

import networkx as nx #NetworkX is a Python language software package for the creation,
                      #manipulation, and study of the structure, dynamics, and functions of complex networks.
from lxml import etree #Tree element analysis module

import sys, os # On importe le module os qui dispose de variables 
               # et de fonctions utiles pour dialoguer avec votre 
               # système d'exploitation

#Import General data
#print('Importing geometrical data')
OS=etree.parse(dir + Project + 'in/' + Gen).getroot().xpath('OS')[0].get("value")
Output_path=etree.parse(dir + Project + 'in/' + Gen).getroot().xpath('Output')[0].get("path")
Paraview=int(etree.parse(dir + Project + 'in/' + Gen).getroot().xpath('Paraview')[0].get("value"))
ParaviewWF=int(etree.parse(dir + Project + 'in/' + Gen).getroot().xpath('Paraview')[0].get("WallFlux"))
ParaviewMF=int(etree.parse(dir + Project + 'in/' + Gen).getroot().xpath('Paraview')[0].get("MembraneFlux"))
ParaviewPF=int(etree.parse(dir + Project + 'in/' + Gen).getroot().xpath('Paraview')[0].get("PlasmodesmataFlux"))
ParaviewWP=int(etree.parse(dir + Project + 'in/' + Gen).getroot().xpath('Paraview')[0].get("WallPot"))
ParaviewCP=int(etree.parse(dir + Project + 'in/' + Gen).getroot().xpath('Paraview')[0].get("CellPot"))
ParTrack=int(etree.parse(dir + Project + 'in/' + Gen).getroot().xpath('ParTrack')[0].get("value"))
Sym_Contagion=int(etree.parse(dir + Project + 'in/' + Gen).getroot().xpath('Sym_Contagion')[0].get("value"))
Apo_Contagion=int(etree.parse(dir + Project + 'in/' + Gen).getroot().xpath('Apo_Contagion')[0].get("value"))
color_threshold=float(etree.parse(dir + Project + 'in/' + Gen).getroot().xpath('color_threshold')[0].get("value"))
thickness_disp=float(etree.parse(dir + Project + 'in/' + Gen).getroot().xpath('thickness_disp')[0].get("value"))
thicknessJunction_disp=float(etree.parse(dir + Project + 'in/' + Gen).getroot().xpath('thicknessJunction_disp')[0].get("value"))
radiusPlasmodesm_disp=float(etree.parse(dir + Project + 'in/' + Gen).getroot().xpath('radiusPlasmodesm_disp')[0].get("value"))

#Import Geometrical data
Plant=etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('Plant')[0].get("value")
path=etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('path')[0].get("value")
im_scale=float(etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('im_scale')[0].get("value"))
Maturityrange=etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('Maturityrange/Maturity')
passage_cell_range=etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('passage_cell_range/passage_cell')
aerenchyma_range=etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('aerenchyma_range/aerenchyma')
passage_cell_ID=[]
for passage_cell in passage_cell_range:
    passage_cell_ID.append(int(passage_cell.get("id")))
PPP=list()
InterCid=list() #Aerenchyma is classified as intercellular space
for aerenchyma in aerenchyma_range:
    InterCid.append(int(aerenchyma.get("id"))) #Cell id starting at 0
InterC_perim1=float(etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('InterC_perim1')[0].get("value"))
InterC_perim2=float(etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('InterC_perim2')[0].get("value"))
InterC_perim3=float(etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('InterC_perim3')[0].get("value"))
InterC_perim4=float(etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('InterC_perim4')[0].get("value"))
InterC_perim5=float(etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('InterC_perim5')[0].get("value"))
kInterC=float(etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('kInterC')[0].get("value"))
cell_per_layer=zeros((2,1))
cell_per_layer[0][0]=float(etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('cell_per_layer')[0].get("cortex"))
cell_per_layer[1][0]=float(etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('cell_per_layer')[0].get("stele"))
diffusion_length=zeros((2,1))
diffusion_length[0][0]=float(etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('diffusion_length')[0].get("cortex"))
diffusion_length[1][0]=float(etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('diffusion_length')[0].get("stele"))
thickness=float(etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('thickness')[0].get("value")) #micron
PD_section=float(etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('PD_section')[0].get("value")) #micron^2
Xylem_pieces=False
if float(etree.parse(dir + Project + 'in/' + Geom).getroot().xpath('Xylem_pieces')[0].get("flag"))==1:
    Xylem_pieces=True

#Import hormone properties
Degrad1=float(etree.parse(dir + Project + 'in/' + Horm).getroot().xpath('Hormone_movement/Degradation_constant_H1')[0].get("value")) #Hormone 1 degradation constant (mol degraded / mol-day)
Diff_PD1=float(etree.parse(dir + Project + 'in/' + Horm).getroot().xpath('Hormone_movement/Diffusivity_PD_H1')[0].get("value")) #Hormone 1 diffusivity constant (cm^2/day)
Diff_PW1=float(etree.parse(dir + Project + 'in/' + Horm).getroot().xpath('Hormone_movement/Diffusivity_PW_H1')[0].get("value")) #Hormone 1 diffusivity constant (cm^2/day)
D2O1=int(etree.parse(dir + Project + 'in/' + Horm).getroot().xpath('Hormone_movement/H1_D2O')[0].get("flag")) #Hormone 1 diffusivity constant (cm^2/day)
Active_transport_range=etree.parse(dir + Project + 'in/' + Horm).getroot().xpath('Hormone_active_transport/carrier_range/carrier')
Sym_source_range=etree.parse(dir + Project + 'in/' + Horm).getroot().xpath('Sym_Contagion/source_range/source')
Sym_Zombie0=[]
for source in Sym_source_range:
    Sym_Zombie0.append(int(source.get("id")))
Sym_cc=[]
for source in Sym_source_range:
    Sym_cc.append(float(source.get("concentration")))
Sym_target_range=etree.parse(dir + Project + 'in/' + Horm).getroot().xpath('Sym_Contagion/target_range/target')
Sym_Target=[]
for target in Sym_target_range:
    Sym_Target.append(int(target.get("id")))
Sym_immune_range=etree.parse(dir + Project + 'in/' + Horm).getroot().xpath('Sym_Contagion/immune_range/immune')
Sym_Immune=[]
for immune in Sym_immune_range:
    Sym_Immune.append(int(immune.get("id")))
Apo_source_range=etree.parse(dir + Project + 'in/' + Horm).getroot().xpath('Apo_Contagion/source_range/source')
Apo_Zombie0=[]
for source in Apo_source_range:
    Apo_Zombie0.append(int(source.get("id")))
Apo_cc=[]
for source in Apo_source_range:
    Apo_cc.append(float(source.get("concentration")))
Apo_target_range=etree.parse(dir + Project + 'in/' + Horm).getroot().xpath('Apo_Contagion/target_range/target')
Apo_Target=[]
for target in Apo_target_range:
    Apo_Target.append(int(target.get("id")))
Apo_immune_range=etree.parse(dir + Project + 'in/' + Horm).getroot().xpath('Apo_Contagion/immune_range/immune')
Apo_Immune=[]
for immune in Apo_immune_range:
    Apo_Immune.append(int(immune.get("id")))
contact_range=etree.parse(dir + Project + 'in/' + Horm).getroot().xpath('Contactrange/Contact')
Contact=[]
for contact in contact_range:
    Contact.append(int(contact.get("id")))

#Import cellset data
tree = etree.parse(dir + 'cellsetdata/' + path) #Maize_Charles\\Maize_pip_cross4.xml') # #Parse literally decrypts the tree element data         SteleOK_high.xml
rootelt = tree.getroot()
Cell2Wall_loop = rootelt.xpath('cells/cell/walls') #Cell2Wall_loop contains cell wall groups info (one group by cell), searched by xpath ("Smart" element identifier)

#Set path
points = rootelt.xpath('walls/wall/points') #points contains the wall elements attributes
Walls_loop = rootelt.xpath('cells/cell/walls/wall') #Walls_loop contains the individual cell to wall associations
Cells_loop=rootelt.xpath('cells/cell') #Cells_loop contains the cell attributes
newpath=dir+Project+Output_path+Plant+'/'
#print('Outputs in '+newpath)
if not os.path.exists(newpath):
    os.makedirs(newpath)

#Initializes network structure
G = nx.Graph() #Full network

#Creates wall & junction nodes
#print('Creating network nodes')
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
                print('Warning null wall segment length! wid:',wid)
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

NwallsJun=Nwalls+jid
position=nx.get_node_attributes(G,'position') #Nodes XY positions (micrometers)

#Junction nodes are pointwise by definition so their length is null, except for junctions at root surface, which are attributed a quarter of the length of each surface neighbouring wall for radial transport 
lengths=nx.get_node_attributes(G,'length') #Walls lengths (micrometers)

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

#Get X and Y for Cell nodes and cell nodes
listsieve=[]
listxyl=[]
listxylwalls=[]
Apo_w_Zombies0=[]
Apo_w_cc=[]
Apo_w_Target=[]
Apo_w_Immune=[]
for w in Cell2Wall_loop: #Loop on cells. Cell2Wall_loop contains cell wall groups info (one group by cell)
	totx=0.0 #Summing up cell walls X positions
	toty=0.0 #Summing up cell walls Y positions
	cellnumber1 = int(w.getparent().get("id")) #Cell ID number
	cgroup=int(w.getparent().get("group")) #Cell type (1=Exodermis;2=epidermis;3=endodermis;4=cortex;5=stele;16=pericycle)
	div=float(len(w)) #Total number of walls around the current cell
	for r in w: #w points to the cell walls around the current cell
	    wid= int(r.get("id")) #Wall ID number
	    totx += position[wid][0] #Contains the walls average X positions
	    toty += position[wid][1] #Contains the walls average Y positions
	finalx=totx/div #Average cell X position (from the average position of its walls)
	finaly=toty/div #Average cell Y position (from the average position of its walls)
	G.add_node(NwallsJun + cellnumber1, indice=(NwallsJun) + cellnumber1, type="cell", position = (finalx,finaly), cgroup=cgroup) #Adding cell nodes  borderlink=0,
	if cgroup==11 or cgroup==23: #Phloem sieve tube
	    listsieve.append(NwallsJun+cellnumber1)
	elif cgroup==13 or cgroup==19 or cgroup==20: #Xylem vessel
	    listxyl.append(NwallsJun+cellnumber1)
	    for r in w: #w points to the cell walls around the current cell
	        wid= int(r.get("id")) #Wall ID number
	        listxylwalls.append(wid) #ghost walls crossing xylem vessels will appear twice
	if Apo_Contagion:
	    if cellnumber1 in Apo_Zombie0:
	        cc=Apo_cc[Apo_Zombie0.index(cellnumber1)]
	        for r in w: #w points to the cell walls around the current cell
	            wid= int(r.get("id")) #Wall ID number
	            if wid not in Apo_w_Zombies0:
	                Apo_w_Zombies0.append(wid)
	                Apo_w_cc.append(cc)
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

#add Edges
#print('Creating network connections')
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
        if dist_cell==0:
            print(dist_cell,cellnumber1)
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

Cell_connec=-ones((Ncells,25),dtype=int) #Connected cells for further ranking
nCell_connec=zeros((Ncells,1),dtype=int) #Quantity of cell to cell connectionsC:\Users\heymansad
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
x_grav/=n_cell_endo
y_grav/=n_cell_endo

#Calculation of a and b parameter for cortex AQP activity radial distribution
#print('Identifying cell layers')
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
#<group id="26" name="CompanionCell" /> is companion cell

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
	elif celltype==26: #Companion cell in new Cellset version
	    celltype=12
	Cell_rank[cellnumber1]=celltype #Later on, cell types 4 and 5 will be updated to account for their ranking within the cortex / stele
	x_cell=position[NwallsJun + cellnumber1][0] #Cell position (micrometers)
	y_cell=position[NwallsJun + cellnumber1][1]
	dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	Layer_dist[celltype]+=dist
	nLayer[celltype]+=1
	if celltype==13:
	    xyl_dist.append([dist])
xyl80_dist=percentile(xyl_dist, 80)

if nLayer[16]==0: #If there is no labelled pericycle
    stele_connec_rank=3 #Endodermis connected to stele cells
else:
    stele_connec_rank=16 #Pericycle connected to stele cells
if nLayer[1]==0: #If there is no labelled exodermis
    outercortex_connec_rank=2 #Cortex connected to epidermis cells
else:
    outercortex_connec_rank=1 #Cortex connected to exodermis cells
#rank_cellperimeters_in=linspace(nan,nan,100) #
#rank_cellperimeters_out=linspace(nan,nan,100)
listprotosieve=[]
mincid=99999
for w in Cell2Wall_loop: #Loop on cells. Cell2Wall_loop contains cell wall groups info (one group by cell)
	cellnumber1 = int(w.getparent().get("id")) #Cell ID number
	#celltype=G.node[NwallsJun + cellnumber1]['cgroup']
	celltype=Cell_rank[cellnumber1] #Cell types 4 and 5 updated to account for their ranking within the cortex / stele
	if celltype==4: #Cortex
	    if any(Cell_rank[Cell_connec[cellnumber1][0:nCell_connec[cellnumber1][0]]]==3): #Cell to cell connection with endodermis
	        Cell_rank[cellnumber1]=40
	        x_cell=position[NwallsJun + cellnumber1][0] #Cell position (micrometers)
	        y_cell=position[NwallsJun + cellnumber1][1]
	        dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	        Layer_dist[40]+=dist
	        nLayer[40]+=1
	        #rank_cellperimeters_in[int(nLayer[40]-1)]=cellperimeter[cellnumber1]
	        if cellperimeter[cellnumber1]<InterC_perim1:
	            InterCid.append(cellnumber1) #Cell id starting at 0
	    elif any(Cell_rank[Cell_connec[cellnumber1][0:nCell_connec[cellnumber1][0]]]==outercortex_connec_rank): #Cell to cell connection with exodermis
	        Cell_rank[cellnumber1]=49
	        x_cell=position[NwallsJun + cellnumber1][0] #Cell position (micrometers)
	        y_cell=position[NwallsJun + cellnumber1][1]
	        dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	        Layer_dist[49]+=dist
	        nLayer[49]+=1
	        #rank_cellperimeters_out[int(nLayer[49]-1)]=cellperimeter[cellnumber1]
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

#cortex_cellperimeters_in=rank_cellperimeters_in #Inner part of the cortex (close to endodermis)
#cortex_cellperimeters_out=rank_cellperimeters_out #Outer part of cortex
for i in range(12):
    #rank_cellperimeters_in=linspace(nan,nan,100)
    #rank_cellperimeters_out=linspace(nan,nan,100)
    for w in Cell2Wall_loop: #Loop on cells. Cell2Wall_loop contains cell wall groups info (one group by cell)
        cellnumber1 = int(w.getparent().get("id")) #Cell ID number
        celltype=Cell_rank[cellnumber1] #Cell types 4 and 5 updated to account for their ranking within the cortex / stele
        if celltype==4 and i<4: #Cortex
	   # if i<4: #Within 3 layers of cortical sides
	        if any(Cell_rank[Cell_connec[cellnumber1][0:nCell_connec[cellnumber1][0]]]==(40+i)): #Cell to cell connection with endodermis
	            Cell_rank[cellnumber1]=41+i
	            x_cell=position[NwallsJun + cellnumber1][0] #Cell position (micrometers)
	            y_cell=position[NwallsJun + cellnumber1][1]
	            dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	            Layer_dist[41+i]+=dist
	            nLayer[41+i]+=1
	            #rank_cellperimeters_in[int(nLayer[41+i]-1)]=cellperimeter[cellnumber1]
	            if i==0 and cellperimeter[cellnumber1]<InterC_perim2:
	                InterCid.append(cellnumber1)
	            elif i==1 and cellperimeter[cellnumber1]<InterC_perim3:
	                InterCid.append(cellnumber1)
	            elif i==2 and cellperimeter[cellnumber1]<InterC_perim4:
	                InterCid.append(cellnumber1)
	            elif i>2 and cellperimeter[cellnumber1]<InterC_perim5:
	                InterCid.append(cellnumber1)
	        elif any(Cell_rank[Cell_connec[cellnumber1][0:nCell_connec[cellnumber1][0]]]==(49-i)): #Cell to cell connection with exodermis
	            Cell_rank[cellnumber1]=48-i
	            x_cell=position[NwallsJun + cellnumber1][0] #Cell position (micrometers)
	            y_cell=position[NwallsJun + cellnumber1][1]
	            dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
	            Layer_dist[48-i]+=dist
	            nLayer[48-i]+=1
	            #rank_cellperimeters_out[int(nLayer[48-i]-1)]=cellperimeter[cellnumber1]
	            if cellperimeter[cellnumber1]<InterC_perim5:
	                InterCid.append(cellnumber1)
        elif celltype==5 or celltype==11 or celltype==12 or celltype==13: #Stele
            if i<10:
                if any(Cell_rank[Cell_connec[cellnumber1][0:nCell_connec[cellnumber1][0]]]==(50+i)): #Cell to cell connection with pericycle
                    Cell_rank[cellnumber1]=51+i
                    x_cell=position[NwallsJun + cellnumber1][0] #Cell position (micrometers)
                    y_cell=position[NwallsJun + cellnumber1][1]
                    dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
                    Layer_dist[51+i]+=dist
                    nLayer[51+i]+=1
            else: #No more than 11 stele cell layers
                Cell_rank[cellnumber1]=61
                x_cell=position[NwallsJun + cellnumber1][0] #Cell position (micrometers)
                y_cell=position[NwallsJun + cellnumber1][1]
                dist=hypot(x_cell-x_grav,y_cell-y_grav) #(micrometers)
                Layer_dist[61]+=dist
                nLayer[61]+=1
#    if i<4:
        #cortex_cellperimeters_in=vstack((cortex_cellperimeters_in,rank_cellperimeters_in))
        #cortex_cellperimeters_out=vstack((rank_cellperimeters_out,cortex_cellperimeters_out))
#cortex_cellperimeters=vstack((cortex_cellperimeters_in,cortex_cellperimeters_out))
InterCid=InterCid[1:]

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
            for neighboor, eattr in edges.items(): #Loop on connections (edges)
                if eattr['path'] == "plasmodesmata" and (G.node[indice[neighboor]]['cgroup']==11 or G.node[indice[neighboor]]['cgroup']==23): #Plasmodesmata connection  #eattr is the edge attribute (i.e. connection type)
                    PPP.append(i-NwallsJun)
        elif G.node[i]['cgroup']==outercortex_connec_rank or G.node[i]['cgroup']==4 or G.node[i]['cgroup']==3: #exodermis or cortex or endodermis (or epidermis if there is no exodermis)
            if i-NwallsJun not in InterCid: #The loop focuses on exo, cortex and endodermis cells that are not intercellular spaces
                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                    if eattr['path'] == "plasmodesmata": #Plasmodesmata connection  #eattr is the edge attribute (i.e. connection type)
                        j = (indice[neighboor]) #neighbouring node number
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
    rank2row[16]=j #Pericycle
    Layer_dist2=vstack((Layer_dist2,Layer_dist[16]))
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

i1=40
nLayer_ref=nLayer[i1]
Layer_dist_ref=Layer_dist[i1]
rank2row[i1]=j
Layer_dist2=vstack((Layer_dist2,Layer_dist[i1]))
j+=1 #number of the row that was just added
i1=41
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
		if G.node[NwallsJun + cellnumber1]['cgroup']==4: #Cortex
		    dmax_cortex=max(dmax_cortex,dist_grav[wid])
		    dmin_cortex=min(dmin_cortex,dist_grav[wid])
		elif G.node[NwallsJun + cellnumber1]['cgroup']==2: #Epidermis
		    davg_epi+=dist_grav[wid]
		    temp+=1.0
davg_epi/=temp #Last step of averaging (note that we take both inner and outer membranes into account in the averaging)
perimeter=2*pi*davg_epi*1.0E-04 #(cm)

#Calculating cell borders accounting for wall thickness
if Paraview==1 or ParTrack==1 or Apo_Contagion>0 or Sym_Contagion>0:
    print('Preparing geometrical properties for Paraview')
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
    Cell2ThickWalls=empty((Ncells,32))
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
            for neighboor, eattr in edges.items(): #Loop on connections (edges)
                cid=int(indice[neighboor])
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
                elif G.node[neighboor]['type']=='apo': #Node j is a junction
                    Wall2Junction[wid][int(nWall2Junction[wid])]=indice[neighboor]
                    nWall2Junction[wid]+=1
            cid1=Wall2Cell[wid][0]
            cid2=Wall2Cell[wid][1]
            rank1=Cell_rank[int(cid1-NwallsJun)]
            row1=rank2row[int(rank1)]
            if not isnan(cid2):
                rank2=Cell_rank[int(cid2-NwallsJun)]
                row2=rank2row[int(rank2)]
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
    L_diff=((abs(float(Layer_dist[3]-Layer_dist[2])*1.0E-04),abs(float(Layer_dist[3] - xyl80_dist)*1.0E-04))) #Diffusion lengths (cm)
    for node, edges in G.adjacency_iter():
        wid=indice[node]
        if wid<Nwalls: #wall that is not a junction (connected to a cell)
            for neighboor, eattr in edges.items(): #Loop on connections (edges)
                cid=int(indice[neighboor])
                if G.node[cid]['type']=='cell':
                    Wall2Cell[wid][int(nWall2Cell[wid])]=cid
                    nWall2Cell[wid]+=1
                elif G.node[neighboor]['type']=='apo': #Node j is a junction
                    Wall2Junction[wid][int(nWall2Junction[wid])]=indice[neighboor]
                    nWall2Junction[wid]+=1
            cid1=Wall2Cell[wid][0]
            cid2=Wall2Cell[wid][1]
            rank1=Cell_rank[int(cid1-NwallsJun)]
            row1=rank2row[int(rank1)]
            if not isnan(cid2):
                rank2=Cell_rank[int(cid2-NwallsJun)]
                row2=rank2row[int(rank2)]
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

#temp=0
#for node, edges in G.adjacency_iter():
#    for neighboor, eattr in edges.items(): #Loop on connections (edges) separated to make sure that Wall2Cell[i] is complete
#        temp+=1
    

#At this point, nThickWallPolygonX equals 2 (two thick wall nodes for each wall node)
#We still have to add two more thick junction nodes for each row in nThickWallPolygonX (each row is a "wall polygon")
if Paraview==1 or ParTrack==1 or Apo_Contagion>0 or Sym_Contagion>0:
    for node, edges in G.adjacency_iter():
        i=indice[node]
        if i<Nwalls: #wall that is not a junction (connected to a cell)
            for neighboor, eattr in edges.items(): #Loop on connections (edges) separated to make sure that Wall2Cell[i] is complete
                j=indice[neighboor]
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
            for neighboor, eattr in edges.items(): #Loop on connections (edges) separated to make sure that Wall2Cell[i] is complete
                j=indice[neighboor]
                if G.node[j]['type']=='apo': #then j is a junction node 
                    for cid in Wall2Cell[i][0:int(nWall2Cell[i])]: #Wall2Cell[i][k] are the cell node ID associated to Wall i
                        if cid not in Junction2Wall2Cell[j-Nwalls]: 
                            Junction2Wall2Cell[j-Nwalls][int(nJunction2Wall2Cell[j-Nwalls])]=cid #Writes the cells indirectly associated to each junction
                            nJunction2Wall2Cell[j-Nwalls]+=1

if Paraview==1 or ParTrack==1 or Apo_Contagion>0 or Sym_Contagion>0:
    Apo_j_Zombies0=[]
    Apo_j_cc=[]
    for j in range(Nwalls, NwallsJun):
        for cid in Junction2Wall2Cell[j-Nwalls]:
            if cid-NwallsJun in Apo_Zombie0:
                cc=Apo_cc[Apo_Zombie0.index(cid-NwallsJun)]
                if j not in Apo_j_Zombies0:
                    Apo_j_Zombies0.append(j)
                    Apo_j_cc.append(cc)

#######################################
##Variables for potential calculation##
#######################################

#Import Hydraulic data
#print('Importing hydraulic data')
kwrange=etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('kwrange/kw')
kw_barrier_range=etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('kw_barrier_range/kw_barrier')
kmb=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('km')[0].get("value"))
kAQPrange=etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('kAQPrange/kAQP')

ratio_cortex=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('ratio_cortex')[0].get("value"))
Kplrange=etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('Kplrange/Kpl')
Fplxheight=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('Fplxheight')[0].get("value")) #Product of plasmodesmatal frequency (per square cm) by cell axial length when the frequency was counted (cm) => units: per cm of cell membrane perimeter, to be multiplied by cell perimeter in cm to obtain a quantity of plasmodesmata
Fplxheight_epi_exo=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('Fplxheight_epi_exo')[0].get("value")) #Because they are "cell axial length independent", these values conserve the quantity of plasmodesmata when cells elongate
Fplxheight_outer_cortex=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('Fplxheight_outer_cortex')[0].get("value"))
Fplxheight_cortex_cortex=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('Fplxheight_cortex_cortex')[0].get("value"))
Fplxheight_cortex_endo=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('Fplxheight_cortex_endo')[0].get("value"))
Fplxheight_endo_endo=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('Fplxheight_endo_endo')[0].get("value"))
Fplxheight_endo_peri=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('Fplxheight_endo_peri')[0].get("value"))
Fplxheight_peri_peri=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('Fplxheight_peri_peri')[0].get("value"))
Fplxheight_peri_stele=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('Fplxheight_peri_stele')[0].get("value"))
Fplxheight_stele_stele=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('Fplxheight_stele_stele')[0].get("value"))
Fplxheight_stele_comp=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('Fplxheight_stele_comp')[0].get("value"))
Fplxheight_peri_comp=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('Fplxheight_peri_comp')[0].get("value"))
Fplxheight_comp_comp=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('Fplxheight_comp_comp')[0].get("value"))
Fplxheight_comp_sieve=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('Fplxheight_comp_sieve')[0].get("value"))
Fplxheight_peri_sieve=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('Fplxheight_peri_sieve')[0].get("value"))
Fplxheight_stele_sieve=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('Fplxheight_stele_sieve')[0].get("value"))
K_sieve=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('K_sieve')[0].get("value")) #Sieve tube hydraulic conductance
K_xyl=float(etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('K_xyl')[0].get("value")) #Xylem vessel axial hydraulic conductance
Xcontactrange=etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('Xcontactrange/Xcontact')
path_hydraulics=etree.parse(dir + Project + 'in/' + Hydr).getroot().xpath('path_hydraulics/Output')

#Import boundary conditions
Psi_soil_range=etree.parse(dir + Project + 'in/' + BC).getroot().xpath('Psi_soil_range/Psi_soil')
BC_xyl_range=etree.parse(dir + Project + 'in/' + BC).getroot().xpath('BC_xyl_range/BC_xyl')
BC_sieve_range=etree.parse(dir + Project + 'in/' + BC).getroot().xpath('BC_sieve_range/BC_sieve')
Psi_cell_range=etree.parse(dir + Project + 'in/' + BC).getroot().xpath('Psi_cell_range/Psi_cell')
Elong_cell_range=etree.parse(dir + Project + 'in/' + BC).getroot().xpath('Elong_cell_range/Elong_cell')
#Elong_cell_kappa=float(etree.parse(dir + Project + 'in/' + BC).getroot().xpath('Elong_cell')[0].get("kappa_dot")) #Rate change of curvature (radian/micron/d)
Water_fraction_apo=float(etree.parse(dir + Project + 'in/' + BC).getroot().xpath('Water_fractions')[0].get("Apoplast")) #Relative volumetric fraction of water in the apoplast (dimensionless)
Water_fraction_sym=float(etree.parse(dir + Project + 'in/' + BC).getroot().xpath('Water_fractions')[0].get("Symplast")) #Relative volumetric fraction of water in the symplast (dimensionless)
path_BC=etree.parse(dir + Project + 'in/' + BC).getroot().xpath('path_scenarios/Output')[0].get("path")
Nhydraulics=len(path_hydraulics) #Total number of hydraulic parameters sets
Nkw=len(kwrange)
NKpl=len(Kplrange)
NkAQP=len(kAQPrange)
NXcontact=len(Xcontactrange)
Nkw_barrier=len(kw_barrier_range)
Nscenarios=len(Psi_soil_range) #Including the "forced scenario"
Nmaturity=len(Maturityrange)
Psi_soil=zeros((2,Nscenarios))
Os_soil=zeros((6,Nscenarios))
Os_xyl=zeros((6,Nscenarios))
C_flag=False #Do we calculate solute stationnary fluxes?
for count in range(1,Nscenarios):
    Psi_soil[0][count]=float(Psi_soil_range[count].get("pressure_left")) #Soil pressure potential (hPa)
    Psi_soil[1][count]=float(Psi_soil_range[count].get("pressure_right"))
    Os_soil[0][count]=float(Psi_soil_range[count].get("osmotic_left"))
    Os_soil[1][count]=float(Psi_soil_range[count].get("osmotic_right"))
    Os_soil[2][count]=float(Psi_soil_range[count].get("osmotic_symmetry"))
    Os_soil[3][count]=float(Psi_soil_range[count].get("osmotic_shape")) #1 for linear, >1 for outer slope flat, <1 for inner slope flat
    Os_soil[4][count]=float(etree.parse(dir + Project + 'in/' + BC).getroot().xpath('Psi_soil_range/osmotic_diffusivity')[0].get("value"))
    #Os_soil[5][count]=float(etree.parse(dir + Project + 'in/' + BC).getroot().xpath('Psi_soil_range/osmotic_convection')[0].get("flag"))
    Os_xyl[0][count]=float(BC_xyl_range[count].get("osmotic_xyl"))
    Os_xyl[1][count]=float(BC_xyl_range[count].get("osmotic_endo"))
    Os_xyl[2][count]=float(BC_xyl_range[count].get("osmotic_symmetry"))
    Os_xyl[3][count]=float(BC_xyl_range[count].get("osmotic_shape")) #The symmetry is central. 1 for linear, >1 for outer slope flat, <1 for inner slope flat
    Os_xyl[4][count]=float(etree.parse(dir + Project + 'in/' + BC).getroot().xpath('BC_xyl_range/osmotic_diffusivity')[0].get("value"))
    #Os_xyl[5][count]=float(etree.parse(dir + Project + 'in/' + BC).getroot().xpath('BC_xyl_range/osmotic_convection')[0].get("flag"))
    if not Os_xyl[4][count]==0 and not Os_soil[4][count]==0:
        C_flag=True
        print('Calculation of analytical solution for radial solute transport in cell walls')

#Unit changes
sperd=24.0*3600.0 #(seconds per day)
cmperm=100.0 #(cm per metre)

#Start the loop of hydraulic properties
for h in range(Nhydraulics):
    #print('   ')
    #print('Hydraulic network #'+str(h))
    newpath=dir+Project+Output_path+Plant+'/'+path_BC+path_hydraulics[h].get("path")
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    
    #System solving
    Psi_xyl=empty((Nmaturity,Nscenarios))
    Psi_xyl[:]=NAN
    dPsi_xyl=empty((Nmaturity,Nscenarios))
    dPsi_xyl[:]=NAN
    iEquil_xyl=NAN #index of the equilibrium root xylem pressure scenario
    Flow_xyl=empty((Nxyl+1,Nscenarios))
    Flow_xyl[:]=NAN
    Psi_sieve=empty((Nmaturity,Nscenarios))
    Psi_sieve[:]=NAN
    dPsi_sieve=empty((Nmaturity,Nscenarios))
    dPsi_sieve[:]=NAN
    iEquil_sieve=NAN #index of the equilibrium root phloem pressure scenario
    Flow_sieve=empty((Nsieve+1,Nscenarios))
    Flow_sieve[:]=NAN
    Os_sieve=zeros((1,Nscenarios))
    Os_cortex=zeros((1,Nscenarios))
    Os_hetero=zeros((1,Nscenarios))
    s_factor=zeros((1,Nscenarios))
    s_hetero=zeros((1,Nscenarios))
    Elong_cell=zeros((1,Nscenarios))
    Elong_cell_side_diff=zeros((1,Nscenarios))
    UptakeLayer_plus=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
    UptakeLayer_minus=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
    Q_xyl_layer=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
    Q_sieve_layer=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
    Q_elong_layer=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
    STFmb=zeros((Nmb,Nmaturity))
    STFcell_plus=zeros((Ncells,Nmaturity))
    STFcell_minus=zeros((Ncells,Nmaturity))
    STFlayer_plus=zeros((int(r_discret[0]),Nmaturity))
    STFlayer_minus=zeros((int(r_discret[0]),Nmaturity))
    PsiCellLayer=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
    PsiWallLayer=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
    OsCellLayer=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
    nOsCellLayer=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
    OsWallLayer=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
    nOsWallLayer=zeros((int(r_discret[0]),Nmaturity,Nscenarios)) #Used for averaging OsWallLayer
    NWallLayer=zeros((int(r_discret[0]),Nmaturity,Nscenarios))
    #UptakeDistri_plus=zeros((40,3,8))#the size will be adjusted, but won't be more than 40. Dimension 1: radial position, 2: compartment, 3: scenario
    #UptakeDistri_minus=zeros((40,3,8))
    Q_tot=zeros((Nmaturity,Nscenarios)) #(cm^3/d) Total flow rate at root surface
    kr_tot=zeros((Nmaturity,1))
    Hydropatterning=empty((Nmaturity,Nscenarios))
    Hydropatterning[:]=NAN
    Hydrotropism=empty((Nmaturity,Nscenarios))
    Hydrotropism[:]=NAN
    
    iMaturity=-1 #Iteration index in the Barriers loop
    for Maturity in Maturityrange:
        Barrier=int(Maturity.get("Barrier")) #Apoplastic barriers (0: No apoplastic barrier, 1:Endodermis radial walls, 2:Endodermis with passage cells, 3: Endodermis full, 4: Endodermis full and exodermis radial walls)
        height=int(Maturity.get("height")) #Cell length in the axial direction (microns)
        
        #Index for barriers loop
        iMaturity+=1
        print('Maturity #'+str(iMaturity)+' with apoplastic barrier type #'+str(Barrier))
        
        #Scenarios concern boundary conditions only
        count=0
        #print('Scenario #'+str(count))
        
        #Soil, xylem, and phloem pressure potentials
        Psi_xyl[iMaturity][count]=float(BC_xyl_range[count].get("pressure")) #Xylem pressure potential (hPa)
        dPsi_xyl[iMaturity][count]=float(BC_xyl_range[count].get("deltaP")) #Xylem pressure potential change as compared to equilibrium pressure (hPa)
        Flow_xyl[0][count]=float(BC_xyl_range[count].get("flowrate")) #Xylem flow rate (cm^3/d)
        if not isnan(Flow_xyl[0][count]):
            if isnan(Psi_xyl[iMaturity][count]) and isnan(dPsi_xyl[iMaturity][count]):
                tot_flow=Flow_xyl[0][count]
                sum_area=0
                i=1
                for cid in listxyl:
                    area=cellarea[cid-NwallsJun]
                    Flow_xyl[i][count]=tot_flow*area
                    sum_area+=area
                    i+=1
                i=1
                for cid in listxyl:
                    Flow_xyl[i][count]/=sum_area #Total xylem flow rate partitioned proportionnally to xylem cross-section area
                    i+=1
                if Flow_xyl[0][count]==0.0:
                    iEquil_xyl=count
            else:
                print('Error: Cannot have both pressure and flow BC at xylem boundary')
        elif not isnan(dPsi_xyl[iMaturity][count]):
            if isnan(Psi_xyl[iMaturity][count]):
                if not isnan(iEquil_xyl):
                    Psi_xyl[iMaturity][count]=Psi_xyl[iMaturity][iEquil_xyl]+dPsi_xyl[iMaturity][count]
                else:
                    print('Error: Cannot have xylem pressure change relative to equilibrium without having a prior scenario with equilibrium xylem boundary condition')
            else:
                print('Error: Cannot have both pressure and pressure change relative to equilibrium as xylem boundary condition')
        
        Psi_sieve[iMaturity][count]=float(BC_sieve_range[count].get("pressure")) #Phloem sieve element pressure potential (hPa)
        dPsi_sieve[iMaturity][count]=float(BC_sieve_range[count].get("deltaP")) #Phloem pressure potential change as compared to equilibrium pressure (hPa)
        Flow_sieve[0][count]=float(BC_sieve_range[count].get("flowrate")) #Phloem flow rate (cm^3/d)
        if not isnan(Flow_sieve[0][count]):
            if isnan(Psi_sieve[iMaturity][count]) and isnan(dPsi_sieve[iMaturity][count]):
                tot_flow=Flow_sieve[0][count]
                sum_area=0
                i=1
                for cid in listprotosieve:
                    area=cellarea[cid-NwallsJun]
                    Flow_sieve[i][count]=tot_flow*area
                    sum_area+=area
                    i+=1
                i=1
                for cid in listprotosieve:
                    Flow_sieve[i][count]/=sum_area #Total phloem flow rate partitioned proportionnally to phloem cross-section area
                    i+=1
                if Flow_sieve[0][count]==0.0:
                    iEquil_sieve=count
            else:
                print('Error: Cannot have both pressure and flow BC at phloem boundary')
        elif not isnan(dPsi_sieve[iMaturity][count]):
            if isnan(Psi_sieve[iMaturity][count]):
                if not isnan(iEquil_sieve):
                    Psi_sieve[iMaturity][count]=Psi_sieve[iMaturity][iEquil_sieve]+dPsi_sieve[iMaturity][count]
                else:
                    print('Error: Cannot have phloem pressure change relative to equilibrium without having a prior scenario with equilibrium phloem boundary condition')
            else:
                print('Error: Cannot have both pressure and pressure change relative to equilibrium as phloem boundary condition')
        
        #Soil - root contact limit
        if NXcontact == Nhydraulics:
            Xcontact=float(Xcontactrange[h].get("value")) #(micrometers) X threshold coordinate of contact between soil and root (lower X not in contact with soil)
        elif NXcontact == 1:
            Xcontact=float(Xcontactrange[0].get("value"))
        else:
            Xcontact=float(Xcontactrange[int(h/(NkAQP*NKpl*Nkw*Nkw_barrier))].get("value")) #OK
        
        #Cell wall hydraulic conductivity
        if Nkw == Nhydraulics:
            kw = float(kwrange[h].get("value"))
        elif Nkw == 1:
            kw = float(kwrange[0].get("value"))
        else:
            kw = float(kwrange[int(h/(NkAQP*NKpl))%Nkw].get("value"))
        if Nkw_barrier == Nhydraulics:
            kw_barrier = float(kw_barrier_range[h].get("value"))
        elif Nkw_barrier == 1:
            kw_barrier = float(kw_barrier_range[0].get("value"))
        else:
            kw_barrier = float(kw_barrier_range[int(h/(NkAQP*NKpl*Nkw))%Nkw_barrier].get("value"))
        #kw_barrier = kw/10.0
        if Barrier==0: #No Casparian strip ###Yet to come: Punctured Casparian strip as in Steudle et al. (1993)
            kw_endo_endo=kw
            kw_puncture=kw
            kw_exo_exo=kw #(cm^2/hPa/d) hydraulic conductivity of the suberised walls between exodermis cells
            kw_cortex_cortex=kw
            kw_endo_peri=kw #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_endo_cortex=kw #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_passage=kw #(cm^2/hPa/d) hydraulic conductivity of passage cells tangential walls
        elif Barrier==1: #Endodermis radial walls
            kw_endo_endo=kw_barrier
            kw_exo_exo=kw #(cm^2/hPa/d) hydraulic conductivity of the suberised walls between exodermis cells
            kw_cortex_cortex=kw
            kw_endo_peri=kw #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_endo_cortex=kw #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_passage=kw #(cm^2/hPa/d) hydraulic conductivity of passage cells tangential walls
        elif Barrier==2: #Endodermis with passage cells
            kw_endo_endo=kw_barrier
            kw_exo_exo=kw #(cm^2/hPa/d) hydraulic conductivity of the suberised walls between exodermis cells
            kw_cortex_cortex=kw
            kw_endo_peri=kw_barrier #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_endo_cortex=kw_barrier #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_passage=kw #(cm^2/hPa/d) hydraulic conductivity of passage cells tangential walls
        elif Barrier==3: #Endodermis full
            kw_endo_endo=kw_barrier
            kw_exo_exo=kw #(cm^2/hPa/d) hydraulic conductivity of the suberised walls between exodermis cells
            kw_cortex_cortex=kw
            kw_endo_peri=kw_barrier #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_endo_cortex=kw_barrier #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_passage=kw_barrier #(cm^2/hPa/d) hydraulic conductivity of passage cells tangential walls
        elif Barrier==4: #Endodermis full and exodermis radial walls
            kw_endo_endo=kw_barrier
            kw_exo_exo=kw_barrier #(cm^2/hPa/d) hydraulic conductivity of the suberised walls between exodermis cells
            kw_cortex_cortex=kw
            kw_endo_peri=kw_barrier #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_endo_cortex=kw_barrier #(cm^2/hPa/d) hydraulic conductivity of the walls between endodermis and pericycle cells
            kw_passage=kw_barrier #(cm^2/hPa/d) hydraulic conductivity of passage cells tangential walls
        
        #Plasmodesmatal hydraulic conductance
        if NKpl == Nhydraulics:
            iPD=h
        elif NKpl == 1:
            iPD=0
        else:
            iPD=int(h/NkAQP)%NKpl
        Kpl = float(Kplrange[iPD].get("value"))
        
        #Contribution of aquaporins to membrane hydraulic conductivity
        if NkAQP == Nhydraulics:
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
        if ratio_cortex==1: #Uniform AQP activity in all cortex membranes
            a_cortex=0.0  #(1/hPa/d)
            b_cortex=kaqp_cortex #(cm/hPa/d)
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
        
        matrix_W = np.zeros(((G.number_of_nodes()),G.number_of_nodes())) #Initializes the Doussan matrix
        if Apo_Contagion==2 and Sym_Contagion==2:
            matrix_C = np.zeros(((G.number_of_nodes()),G.number_of_nodes())) #Initializes the matrix of convection diffusion
            rhs_C = np.zeros((G.number_of_nodes(),1)) #Initializing the right-hand side matrix of solute apoplastic concentrations
            for i in range(Nwalls):
                if i in Apo_w_Zombies0:
                    matrix_C[i][i]=1.0
                    rhs_C[i][0]=Apo_w_cc[Apo_w_Zombies0.index(i)] #1.0 #Concentration in source wall i defined in Geom
                else: #Decomposition rate (mol decomp/mol-day * cm^3)
                    matrix_C[i][i]-=Degrad1*1.0E-12*(lat_dists[i][0]*thickness*lengths[i]+height*thickness*lengths[i]/2-square(thickness)*lengths[i])
            for j in range(Nwalls,NwallsJun):
                if j in Apo_j_Zombies0:
                    matrix_C[j][j]=1.0
                    rhs_C[j][0]=Apo_j_cc[Apo_j_Zombies0.index(j)] #1.0 #Concentration in source junction j defined in Geom
                else: #Decomposition rate (mol decomp/mol-day * cm^3)
                    matrix_C[j][j]-=Degrad1*1.0E-12*height*thickness*lengths[j]/2
            for cellnumber1 in range(Ncells):
                if cellnumber1 in Sym_Zombie0:
                    matrix_C[NwallsJun+cellnumber1][NwallsJun+cellnumber1]=1.0
                    rhs_C[NwallsJun+cellnumber1][0]=Sym_cc[Sym_Zombie0.index(cellnumber1)] #1.0 #Concentration in source protoplasts defined in Geom
                else: #Decomposition rate (mol decomp/mol-day * cm^3)
                    matrix_C[NwallsJun+cellnumber1][NwallsJun+cellnumber1]-=Degrad1*1.0E-12*cellarea[cellnumber1]*height
        elif Apo_Contagion==2:
            matrix_ApoC = np.zeros(((NwallsJun),NwallsJun)) #Initializes the matrix of convection
            rhs_ApoC = np.zeros((NwallsJun,1)) #Initializing the right-hand side matrix of solute apoplastic concentrations
            for i in range(Nwalls):
                if i in Apo_w_Zombies0:
                    matrix_ApoC[i][i]=1.0
                    rhs_ApoC[i][0]=Apo_w_cc[Apo_w_Zombies0.index(i)] #1 #Concentration in source wall i equals 1 by default
                else: #Decomposition rate (mol decomp/mol-day * cm^3)
                    matrix_ApoC[i][i]-=Degrad1*1.0E-12*(lat_dists[i][0]*thickness*lengths[i]+height*thickness*lengths[i]/2-square(thickness)*lengths[i])
            for j in range(Nwalls,NwallsJun):
                if j in Apo_j_Zombies0:
                    matrix_ApoC[j][j]=1.0
                    rhs_ApoC[j][0]=Apo_j_cc[Apo_j_Zombies0.index(j)] #1 #Concentration in source junction j equals 1 by default
                else: #Decomposition rate (mol decomp/mol-day * cm^3)
                    matrix_ApoC[j][j]-=Degrad1*1.0E-12*height*thickness*lengths[j]/2
        elif Sym_Contagion==2:
            matrix_SymC = np.zeros(((Ncells),Ncells)) #Initializes the matrix of convection
            rhs_SymC = np.zeros((Ncells,1)) #Initializing the right-hand side matrix of solute symplastic concentrations
            for cellnumber1 in range(Ncells):
                if cellnumber1 in Sym_Zombie0:
                    matrix_SymC[cellnumber1][cellnumber1]=1.0
                    rhs_SymC[cellnumber1][0]=Sym_cc[Sym_Zombie0.index(cellnumber1)] #1 #Concentration in source protoplasts equals 1 by default
                else: #Decomposition rate (mol decomp/mol-day * cm^3)
                    matrix_SymC[cellnumber1][cellnumber1]-=Degrad1*1.0E-12*cellarea[cellnumber1]*height
        
        Kmb=zeros((Nmb,1)) #Stores membranes conductances for the second K loop
        jmb=0 #Index of membrane in Kmb
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
                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                    if eattr['path'] == 'membrane': #Wall connection
                        if any(passage_cell_ID==array((indice[neighboor])-NwallsJun)):
                            count_passage+=1
                        if any(InterCid==array((indice[neighboor])-NwallsJun)):
                            count_interC+=1
                            if count_interC==2 and i not in list_ghostwalls:
                                list_ghostwalls.append(i)
                        if G.node[neighboor]['cgroup']==3:#Endodermis
                            count_endo+=1
                        elif G.node[neighboor]['cgroup']==13 or G.node[neighboor]['cgroup']==19 or G.node[neighboor]['cgroup']==20:#Xylem cell or vessel
                            count_xyl+=1
                            if (count_xyl==2 and Xylem_pieces) and i not in list_ghostwalls:
                                list_ghostwalls.append(i)
                        elif G.node[neighboor]['cgroup']>4:#Pericycle or stele but not xylem
                            count_stele_overall+=1
                        elif G.node[neighboor]['cgroup']==4:#Cortex
                            count_cortex+=1
                        elif G.node[neighboor]['cgroup']==1:#Exodermis
                            count_exo+=1
                        elif G.node[neighboor]['cgroup']==2:#Epidermis
                            count_epi+=1
            
            for neighboor, eattr in edges.items(): #Loop on connections (edges)
                j = (indice[neighboor]) #neighbouring node number
                if j > i: #Only treating the information one way to save time
                    path = eattr['path'] #eattr is the edge attribute (i.e. connection type)
                    if path == 'wall': #Wall connection
                        #K = eattr['kw']*1.0E-04*((eattr['lat_dist']+height)*eattr['thickness']-square(eattr['thickness']))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)
                        temp=1.0E-04*((eattr['lat_dist']+height)*thickness-square(thickness))/eattr['length'] #Wall section to length ratio (cm)
                        if (count_interC>=2 and Barrier>0) or (count_xyl==2 and Xylem_pieces): #"Fake wall" splitting an intercellular space or a xylem cell in two
                            K = 1.0E-16 #Non conductive
                            if j not in list_ghostjunctions:
                                fakeJ=True
                                for ind in range(int(nJunction2Wall[j-Nwalls])):
                                    if Junction2Wall[j-Nwalls][ind] not in list_ghostwalls:
                                        fakeJ=False #If any of the surrounding walls is real, the junction is real
                                if fakeJ:
                                    list_ghostjunctions.append(j)
                                    nGhostJunction2Wall+=int(nJunction2Wall[j-Nwalls])+2 #The first and second thick junction nodes each appear twice in the text file for Paraview
                        elif count_cortex>=2: #wall between two cortical cells
                            K = kw_cortex_cortex*temp #Junction-Wall conductance (cm^3/hPa/d)
                        elif count_endo>=2: #wall between two endodermis cells
                            K = kw_endo_endo*temp #Junction-Wall conductance (cm^3/hPa/d)  #(height*eattr['thickness'])/eattr['length']#
                        elif count_stele_overall>0 and count_endo>0: #wall between endodermis and pericycle
                            if count_passage>0:
                                K = kw_passage*temp #(height*eattr['thickness'])/eattr['length']#
                            else:
                                K = kw_endo_peri*temp #Junction-Wall conductance (cm^3/hPa/d) #(height*eattr['thickness'])/eattr['length']#
                        elif count_stele_overall==0 and count_endo==1: #wall between endodermis and cortex
                            if count_passage>0:
                                K = kw_passage*temp  #(height*eattr['thickness'])/eattr['length']#
                            else:
                                K = kw_endo_cortex*temp #Junction-Wall conductance (cm^3/hPa/d)  #(height*eattr['thickness'])/eattr['length']#
                        elif count_exo>=2: #wall between two exodermis cells
                            K = kw_exo_exo*temp #Junction-Wall conductance (cm^3/hPa/d)  #(height*eattr['thickness'])/eattr['length']#
                        else: #other walls
                            K = kw*temp #Junction-Wall conductance (cm^3/hPa/d)  #(height*eattr['thickness'])/eattr['length']#
                        ########Solute fluxes (diffusion across walls and junctions)
                        if Apo_Contagion==2:
                            temp_factor=1.0 #Factor for reduced diffusion across impermeable walls
                            if (count_interC>=2 and Barrier>0) or (count_xyl==2 and Xylem_pieces): #"fake wall" splitting an intercellular space or a xylem cell in two
                                temp_factor=1.0E-16 #Correction
                            elif count_endo>=2:
                                temp_factor=kw_endo_endo/kw
                            elif count_stele_overall>0 and count_endo>0: #wall between endodermis and pericycle
                                if count_passage>0:
                                    temp_factor=kw_passage/kw #(height*eattr['thickness'])/eattr['length']#
                                else:
                                    temp_factor=kw_endo_peri/kw #Junction-Wall conductance (cm^3/hPa/d) #(height*eattr['thickness'])/eattr['length']#
                            elif count_stele_overall==0 and count_endo==1: #wall between endodermis and cortex
                                if count_passage>0:
                                    temp_factor=kw_passage/kw  #(height*eattr['thickness'])/eattr['length']#
                                else:
                                    temp_factor=kw_endo_cortex/kw #Junction-Wall conductance (cm^3/hPa/d)  #(height*eattr['thickness'])/eattr['length']#
                            elif count_exo>=2: #wall between two exodermis cells
                                temp_factor=kw_exo_exo/kw #Junction-Wall conductance (cm^3/hPa/d)  #(height*eattr['thickness'])/eattr['length']#
                            DF=temp*temp_factor*Diff_PW1 #"Diffusive flux" (cm^3/d) temp is the section to length ratio of the wall to junction path
                            if Sym_Contagion==2: #Sym & Apo contagion
                                if i not in Apo_w_Zombies0:
                                    matrix_C[i][i] -= DF
                                    matrix_C[i][j] += DF #Convection will be dealt with further down
                                if j not in Apo_j_Zombies0:
                                    matrix_C[j][j] -= DF #temp_factor is the factor for reduced diffusion across impermeable walls
                                    matrix_C[j][i] += DF
                            else: #Only Apo contagion
                                if i not in Apo_w_Zombies0:
                                    matrix_ApoC[i][i] -= DF
                                    matrix_ApoC[i][j] += DF
                                if j not in Apo_j_Zombies0:
                                    matrix_ApoC[j][j] -= DF #Convection will be dealt with further down
                                    matrix_ApoC[j][i] += DF
                    elif path == "membrane": #Membrane connection
                        #K = (eattr['kmb']+eattr['kaqp'])*1.0E-08*(height+eattr['dist'])*eattr['length']
                        if Apo_Contagion==2 and Sym_Contagion==2:
                            for carrier in Active_transport_range:
                                if int(carrier.get("tissue"))==G.node[j]['cgroup']:
                                    #Condition is that the protoplast (j) is an actual protoplast with membranes
                                    if j-NwallsJun not in InterCid and not (Barrier>0 and (G.node[j]['cgroup']==13 or G.node[j]['cgroup']==19 or G.node[j]['cgroup']==20)):
                                        temp=float(carrier.get("constant"))*(height+eattr['dist'])*eattr['length'] #Linear transport constant (Vmax/KM) [liter/day^-1/micron^-2] * membrane surface [micron²]
                                        if int(carrier.get("direction"))==1: #Influx transporter
                                            if j-NwallsJun not in Sym_Zombie0: #Concentration not affected if set as boundary condition
                                                matrix_C[j][i] += temp #Increase of concentration in protoplast (j) depends on concentration in cell wall (i)
                                            if i not in Apo_w_Zombies0: #Concentration not affected if set as boundary condition
                                                matrix_C[i][i] -= temp #Decrease of concentration in apoplast (i) depends on concentration in apoplast (i)
                                        elif int(carrier.get("direction"))==int(-1): #Efflux transporter
                                            if j-NwallsJun not in Sym_Zombie0: #Concentration not affected if set as boundary condition
                                                matrix_C[j][j] -= temp #Increase of concentration in protoplast (j) depends on concentration in protoplast (j)
                                            if i not in Apo_w_Zombies0: #Concentration not affected if set as boundary condition
                                                matrix_C[i][j] += temp #Decrease of concentration in apoplast (i) depends on concentration in protoplast (j)
                                        else:
                                            error('Error, carrier direction is either 1 (influx) or -1 (efflux), please correct in *_Hormones_Carriers_*.xml')
                        if G.node[j]['cgroup']==1: #Exodermis
                            kaqp=kaqp_exo
                        elif G.node[j]['cgroup']==2: #Epidermis
                            kaqp=kaqp_epi
                        elif G.node[j]['cgroup']==3: #Endodermis
                            kaqp=kaqp_endo
                        elif G.node[j]['cgroup']==13 or G.node[j]['cgroup']==19 or G.node[j]['cgroup']==20: #xylem cell or vessel
                            if Barrier>0: #Xylem vessel
                                kaqp=kaqp_stele*10000 #No membrane resistance because no membrane
                                if Apo_Contagion==2 and Sym_Contagion==2:
                                    #Diffusion between mature xylem vessels and their walls
                                    temp=1.0E-04*(lengths[i]*height)/thickness #Section to length ratio (cm) for the xylem wall
                                    if i not in Apo_w_Zombies0:
                                        matrix_C[i][i] -= temp*Diff_PW1
                                        matrix_C[i][j] += temp*Diff_PW1
                                    if j-NwallsJun not in Sym_Zombie0: #Mature xylem vessels are referred to as cells, so they are on the Sym side even though they are part of the apoplast
                                        matrix_C[j][j] -= temp*Diff_PW1
                                        matrix_C[j][i] += temp*Diff_PW1
                            else:
                                kaqp=kaqp_stele
                        elif G.node[j]['cgroup']>4: #Stele and pericycle but not xylem
                            kaqp=kaqp_stele
                        elif (j-NwallsJun in InterCid) and Barrier>0: #the neighbour is an intercellular space "cell". Between j and i connected by a membrane, only j can be cell because j>i
                            kaqp=kInterC
                            #No carrier
                        elif G.node[j]['cgroup']==4: #Cortex
                            kaqp=float(a_cortex*dist_grav[i]*1.0E-04+b_cortex) #AQP activity (cm/hPa/d)
                            if kaqp < 0:
                                error('Error, negative kaqp in cortical cell, adjust Paqp_cortex')
                        #Calculating each conductance
                        if count_endo>=2: #wall between two endodermis cells, in this case the suberized wall can limit the transfer of water between cell and wall
                            if kw_endo_endo==0.00:
                                K=0.00
                            else:
                                K = 1/(1/(kw_endo_endo/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(height+eattr['dist'])*eattr['length']
                        elif count_exo>=2: #wall between two exodermis cells, in this case the suberized wall can limit the transfer of water between cell and wall
                            if kw_exo_exo==0.00:
                                K=0.00
                            else:
                                K = 1/(1/(kw_exo_exo/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(height+eattr['dist'])*eattr['length']
                        elif count_stele_overall>0 and count_endo>0: #wall between endodermis and pericycle, in this case the suberized wall can limit the transfer of water between cell and wall
                            if count_passage>0:
                                K = 1/(1/(kw_passage/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(height+eattr['dist'])*eattr['length']
                            else:
                                if kw_endo_peri==0.00:
                                    K=0.00
                                else:
                                    K = 1/(1/(kw_endo_peri/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(height+eattr['dist'])*eattr['length']
                        elif count_stele_overall==0 and count_endo==1: #wall between cortex and endodermis, in this case the suberized wall can limit the transfer of water between cell and wall
                            if kaqp==0.0:
                                K=1.00E-16
                            else:
                                if count_passage>0:
                                    K = 1/(1/(kw_passage/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(height+eattr['dist'])*eattr['length']
                                else:
                                    if kw_endo_cortex==0.00:
                                        K=0.00
                                    else:
                                        K = 1/(1/(kw_endo_cortex/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(height+eattr['dist'])*eattr['length']
                        else:
                            if kaqp==0.0:
                                K=1.00E-16
                            else:
                                K = 1/(1/(kw/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(height+eattr['dist'])*eattr['length']
                        Kmb[jmb]=K
                        #if jmb<=10:
                        #    print(jmb,'K init',K,'wid',i,'cid',j-NwallsJun)
                        jmb+=1
                    elif path == "plasmodesmata": #Plasmodesmata connection
                        cgroupi=G.node[i]['cgroup']
                        cgroupj=G.node[j]['cgroup']
                        if cgroupi==19 or cgroupi==20:  #Xylem in new Cellset version
                            cgroupi=13
                        elif cgroupi==21: #Xylem Pole Pericyle in new Cellset version
                            cgroupi=16
                        elif cgroupi==23: #Phloem in new Cellset version
                            cgroupi==11
                        elif cgroupi==26: #Companion Cell in new Cellset version
                            cgroupi==12
                        if cgroupj==19 or cgroupj==20:  #Xylem in new Cellset version
                            cgroupj=13
                        elif cgroupj==21: #Xylem Pole Pericyle in new Cellset version
                            cgroupj=16
                        elif cgroupj==23: #Phloem in new Cellset version
                            cgroupj==11
                        elif cgroupj==26: #Companion Cell in new Cellset version
                            cgroupj==12
                        temp_factor=1.0 #Quantity of plasmodesmata (adjusted by relative aperture)
                        if ((j-NwallsJun in InterCid) or (i-NwallsJun in InterCid)) and Barrier>0: #one of the connected cells is an intercellular space "cell".
                            temp_factor=0.0
                        elif cgroupj==13 and cgroupi==13: #Fake wall splitting a xylem cell or vessel, high conductance in order to ensure homogeneous pressure within the splitted cell
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
                                if i-NwallsJun not in Sym_Zombie0:
                                    matrix_C[i][i] -= DF
                                    matrix_C[i][j] += DF #Convection will be dealt with further down
                                if j-NwallsJun not in Sym_Zombie0:
                                    matrix_C[j][j] -= DF
                                    matrix_C[j][i] += DF
                            else: #Only Sym contagion
                                if i-NwallsJun not in Sym_Zombie0:
                                    matrix_SymC[i-NwallsJun][i-NwallsJun] -= DF
                                    matrix_SymC[i-NwallsJun][j-NwallsJun] += DF
                                if j-NwallsJun not in Sym_Zombie0:
                                    matrix_SymC[j-NwallsJun][j-NwallsJun] -= DF #Convection will be dealt with further down
                                    matrix_SymC[j-NwallsJun][i-NwallsJun] += DF
                    matrix_W[i][i] -= K #Filling the Doussan matrix (symmetric)
                    matrix_W[i][j] += K
                    matrix_W[j][i] += K
                    matrix_W[j][j] -= K
        
        #Adding matrix components at soil-wall and wall-xylem connections & rhs terms
        rhs = np.zeros((G.number_of_nodes(),1))
        rhs_s = np.zeros((G.number_of_nodes(),1)) #Initializing the right-hand side matrix of soil pressure potentials
        rhs_x = np.zeros((G.number_of_nodes(),1)) #Initializing the right-hand side matrix of xylem pressure potentials
        rhs_p = np.zeros((G.number_of_nodes(),1)) #Initializing the right-hand side matrix of hydrostatic potentials for phloem BC
        
        #Adding matrix components at soil-wall connections
        for wid in Borderwall:
            if (position[wid][0]>=Xcontact) or (Wall2Cell[wid][0]-NwallsJun in Contact): #Wall (not including junctions) connected to soil
                temp=1.0E-04*(lengths[wid]/2*height)/(thickness/2)
                K=kw*temp #Half the wall length is used here as the other half is attributed to the junction (Only for connection to soil)
                matrix_W[wid][wid] -= K #Doussan matrix
                rhs_s[wid][0] = -K    #Right-hand side vector, could become Psi_soil[idwall], which could be a function of the horizontal position
                #if C_flag:
                #    #Diffusion
                #    matrix_C[wid][wid] -= temp*Diff1
                #    rhs_C[wid][0] -= temp*Diff1*Os_soil[0][count]
                    
        #Adding matrix components at soil-junction connections
        for jid in Borderjunction:
            if (position[jid][0]>=Xcontact) or (Junction2Wall2Cell[jid-Nwalls][0]-NwallsJun in Contact) or (Junction2Wall2Cell[jid-Nwalls][1]-NwallsJun in Contact) or (Junction2Wall2Cell[jid-Nwalls][2]-NwallsJun in Contact): #Junction connected to soil
                temp=1.0E-04*(lengths[jid]*height)/(thickness/2)
                K=kw*temp
                matrix_W[jid][jid] -= K #Doussan matrix
                rhs_s[jid][0] = -K    #Right-hand side vector, could become Psi_soil[idwall], which could be a function of the horizontal position
                #if C_flag:
                #    matrix_C[jid][jid] -= temp*Diff1 #Diffusion BC at soil junction
                #    rhs_C[jid][0] -= temp*Diff1*Os_soil[0][count]
        
        #Creating connections to xylem & phloem BC elements for kr calculation (either xylem or phloem flow occurs depending on whether the segment is in the differentiation or elongation zone)
        if Barrier>0:
            if not isnan(Psi_xyl[iMaturity][count]): #Pressure xylem BC
                for cid in listxyl:
                    rhs_x[cid][0] = -K_xyl  #Axial conductance of xylem vessels
                    matrix_W[cid][cid] -= K_xyl
                    #if C_flag:
                    #    temp=10E-04*((cellperimeter[cid-NwallsJun]/2)**2)/pi/height #Cell approximative cross-section area (cm^2) per length (cm)
                    #    matrix_C[cid][cid] -= temp*Diff1*100 #Diffusion BC in xylem open vessels assumed 100 times easier than in walls
                    #    rhs_C[cid][0] -= temp*Diff1*100
                rhs = rhs_s*Psi_soil[0][count] + rhs_x*Psi_xyl[iMaturity][count] #multiplication of rhs components delayed till this point so that rhs_s & rhs_x can be re-used to calculate Q
            elif not isnan(Flow_xyl[0][count]):
                i=1
                for cid in listxyl:
                    rhs_x[cid][0] = Flow_xyl[i][count]
                    i+=1
                #    if C_flag:
                #        temp=10E-04*((cellperimeter[cid-NwallsJun]/2)**2)/pi/height #Cell approximative cross-section area (cm^2) per length (cm)
                #        matrix_C[cid][cid] -= temp*Diff1*100 #Diffusion BC in xylem open vessels assumed 100 times easier than in walls
                #        rhs_C[cid][0] -= temp*Diff1*100
                rhs = rhs_s*Psi_soil[0][count] + rhs_x #multiplication of rhs components delayed till this point so that rhs_s & rhs_x can be re-used to calculate Q
            else:
                rhs = rhs_s*Psi_soil[0][count]
        elif Barrier==0:
            if not isnan(Psi_sieve[iMaturity][count]):
                for cid in listprotosieve:
                    rhs_p[cid][0] = -K_sieve  #Axial conductance of phloem sieve tube
                    matrix_W[cid][cid] -= K_sieve
                rhs = rhs_s*Psi_soil[0][count] + rhs_p*Psi_sieve[iMaturity][count] #multiplication of rhs components delayed till this point so that rhs_s & rhs_x can be re-used to calculate Q
            elif not isnan(Flow_sieve[0][count]):
                i=1
                for cid in listprotosieve:
                    rhs_p[cid][0] = Flow_sieve[i][count]
                    i+=1
                rhs = rhs_s*Psi_soil[0][count] + rhs_p #multiplication of rhs components delayed till this point so that rhs_s & rhs_x can be re-used to calculate Q
            else:
                rhs = rhs_s*Psi_soil[0][count]
        
        
        ##################################################
        ##Solve Doussan equation, results in soln matrix##
        ##################################################
        
        soln = np.linalg.solve(matrix_W,rhs) #Solving the equation to get potentials inside the network
        
        #Verification that computation was correct
        verif1=np.allclose(np.dot(matrix_W,soln),rhs)
        
        #print("Correct computation on PSI ?", verif1)
        
        #Removing xylem and phloem BC terms
        if Barrier>0:
            if not isnan(Psi_xyl[iMaturity][count]): #Pressure xylem BC
                for cid in listxyl:
                    matrix_W[cid][cid] += K_xyl
        elif Barrier==0:
            if not isnan(Psi_sieve[iMaturity][count]):
                for cid in listprotosieve:
                    matrix_W[cid][cid] += K_sieve
        
        #Flow rates at interfaces
        Q_soil=[]
        for ind in Borderwall:
            Q_soil.append(rhs_s[ind]*(soln[ind]-Psi_soil[0][count])) #(cm^3/d) Positive for water flowing into the root
        for ind in Borderjunction:
            Q_soil.append(rhs_s[ind]*(soln[ind]-Psi_soil[0][count])) #(cm^3/d) Positive for water flowing into the root
        Q_xyl=[]
        Q_sieve=[]
        if Barrier>0:
            if not isnan(Psi_xyl[iMaturity][count]):
                for cid in listxyl:
                    Q=rhs_x[cid]*(soln[cid]-Psi_xyl[iMaturity][count])
                    Q_xyl.append(Q) #(cm^3/d) Negative for water flowing into xylem tubes
                    rank=int(Cell_rank[cid-NwallsJun])
                    row=int(rank2row[rank])
                    Q_xyl_layer[row][iMaturity][count] += Q
            elif not isnan(Flow_xyl[0][count]):
                for cid in listxyl:
                    Q=-rhs_x[cid]
                    Q_xyl.append(Q) #(cm^3/d) Negative for water flowing into xylem tubes
                    rank=int(Cell_rank[cid-NwallsJun])
                    row=int(rank2row[rank])
                    Q_xyl_layer[row][iMaturity][count] += Q
        elif Barrier==0:
            if not isnan(Psi_sieve[iMaturity][count]):
                for cid in listprotosieve:
                    Q=rhs_p[cid]*(soln[cid]-Psi_sieve[iMaturity][count])
                    Q_sieve.append(Q) #(cm^3/d) Negative for water flowing into phloem tubes
                    rank=int(Cell_rank[cid-NwallsJun])
                    row=int(rank2row[rank])
                    Q_sieve_layer[row][iMaturity][count] += Q
            elif not isnan(Flow_sieve[0][count]):
                for cid in listprotosieve:
                    Q=-rhs_p[cid]
                    Q_sieve.append(Q) #(cm^3/d) Negative for water flowing into xylem tubes
                    rank=int(Cell_rank[cid-NwallsJun])
                    row=int(rank2row[rank])
                    Q_sieve_layer[row][iMaturity][count] += Q
            
        Q_tot[iMaturity][0]=sum(Q_soil) #Total flow rate at root surface
        if Barrier>0:
            if not isnan(Psi_xyl[iMaturity][0]):
                kr_tot[iMaturity][0]=Q_tot[iMaturity][0]/(Psi_soil[0][0]-Psi_xyl[iMaturity][0])/perimeter/height/1.0E-04
            else:
                print('Error: Scenario 0 should have xylem pressure boundary conditions, except for the elongation zone')
        elif Barrier==0:
            if not isnan(Psi_sieve[iMaturity][0]):
                kr_tot[iMaturity][0]=Q_tot[iMaturity][0]/(Psi_soil[0][0]-Psi_sieve[iMaturity][0])/perimeter/height/1.0E-04
            else:
                print('Error: Scenario 0 should have phloem pressure boundary conditions in the elongation zone')
        #print("Flow rates per unit root length: soil ",(Q_tot[iMaturity][0]/height/1.0E-04),"cm^2/d, xylem ",(sum(Q_xyl)/height/1.0E-04),"cm^2/d, phloem ",(sum(Q_sieve)/height/1.0E-04),"cm^2/d")
        #print("Mass balance error:",(Q_tot[iMaturity][0]+sum(Q_xyl)+sum(Q_sieve))/height/1.0E-04,"cm^2/d")
        print("Radial conductivity:",kr_tot[iMaturity][0],"cm/hPa/d")#, Barrier:",Barrier,", height: ",height," microns")
        
        if Barrier>0 and isnan(Psi_xyl[iMaturity][0]):
            Psi_xyl[iMaturity][0]=0.0
            for cid in listxyl:
                Psi_xyl[iMaturity][0]+=soln[cid]/Nxyl #Average of xylem water pressures
        elif Barrier==0 and isnan(Psi_sieve[iMaturity][0]):
            Psi_sieve[iMaturity][0]=0.0
            for cid in listprotosieve:
                Psi_sieve[iMaturity][0]+=soln[cid]/Nprotosieve #Average of protophloem water pressures
        
        #Calculation of standard transmembrane fractions
        jmb=0 #Index for membrane conductance vector
        for node, edges in G.adjacency_iter() : #adjacency_iter returns an iterator of (node, adjacency dict) tuples for all nodes. This is the fastest way to look at every edge. For directed graphs, only outgoing adjacencies are included.
            i = indice[node] #Node ID number
            if i<Nwalls: #wall ID
                psi = soln[i]    #Node water potential
                #print('i',i,'psi',psi)
                psi_o_cell = inf #Opposite cell water potential
                #Here we count surrounding cell types in order to identify in which row of the endodermis or exodermis we are.
                count_endo=0 #total number of endodermis cells around the wall
                count_stele_overall=0 #total number of stelar cells around the wall
                count_exo=0 #total number of exodermis cells around the wall
                count_epi=0 #total number of epidermis cells around the wall
                #count_stele=0 #total number of epidermis cells around the wall
                count_cortex=0 #total number of epidermis cells around the wall
                count_passage=0 #total number of passage cells around the wall
                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                    if eattr['path'] == 'membrane': #Wall connection
                        if any(passage_cell_ID==array((indice[neighboor])-NwallsJun)):
                            count_passage+=1
                        if G.node[neighboor]['cgroup']==3:#Endodermis
                            count_endo+=1
                        elif G.node[neighboor]['cgroup']>4:#Pericycle or stele
                            count_stele_overall+=1
                        elif G.node[neighboor]['cgroup']==1:#Exodermis
                            count_exo+=1
                        elif G.node[neighboor]['cgroup']==2:#Epidermis
                            count_epi+=1
                        elif G.node[neighboor]['cgroup']==4:#Cortex
                            count_cortex+=1
                    # if G.node[neighboor]['cgroup']==5:#Stele
                    #     count_stele+=1
                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                    j = indice[neighboor] #Neighbouring node ID number
                    path = eattr['path'] #eattr is the edge attribute (i.e. connection type)
                    if path == "membrane": #Membrane connection
                        psin = soln[j] #Neighbouring node water potential
                        #print('j',j,'psin',psin)
                        K=Kmb[jmb]
                        #if jmb<=10:
                        #    print(jmb,'K STF',K,'wid',i,'cid',j-NwallsJun)
                        jmb+=1
                        #Flow densities calculation
                        #Macroscopic distributed parameter for transmembrane flow
                        #Discretization based on cell layers and apoplasmic barriers
                        rank = int(Cell_rank[j-NwallsJun])
                        row = int(rank2row[rank])
                        if rank == 1 and count_epi > 0: #Outer exodermis
                            row += 1
                        if rank == 3 and count_cortex > 0: #Outer endodermis
                            if any(passage_cell_ID==array(j-NwallsJun)) and Barrier==2:
                                row += 2
                            else:
                                row += 3
                        elif rank == 3 and count_stele_overall > 0: #Inner endodermis
                            if any(passage_cell_ID==array(j-NwallsJun)) and Barrier==2:
                                row += 1
                        Flow = K * (psi - psin) #Note that this is only valid because we are in the scenario 0 with no osmotic potentials
                        #print('Flow',Flow,'dP',soln[node] - soln[neighboor],'Pi',soln[node],'Pj',soln[neighboor])
                        if ((j-NwallsJun not in InterCid) and (j not in listxyl)) or Barrier==0: #Not part of STF if crosses an intercellular space "membrane" or mature xylem "membrane" (that is no membrane though still labelled like one)
                            if Flow > 0 :
                                UptakeLayer_plus[row][iMaturity][count] += Flow #grouping membrane flow rates in cell layers
                            else:
                                UptakeLayer_minus[row][iMaturity][count] += Flow
                            if Flow/Q_tot[iMaturity][0] > 0 :
                                STFlayer_plus[row][iMaturity] += Flow/Q_tot[iMaturity][0] #Cell standard transmembrane fraction (positive)
                                STFcell_plus[j-NwallsJun][iMaturity] += Flow/Q_tot[iMaturity][0] #Cell standard transmembrane fraction (positive)
                                #STFmb[jmb-1][iMaturity] = Flow/Q_tot[iMaturity][0]
                            else:
                                STFlayer_minus[row][iMaturity] += Flow/Q_tot[iMaturity][0] #Cell standard transmembrane fraction (negative)
                                STFcell_minus[j-NwallsJun][iMaturity] += Flow/Q_tot[iMaturity][0] #Cell standard transmembrane fraction (negative)
                                #STFmb[jmb-1][iMaturity] = Flow/Q_tot[iMaturity][0]
                            STFmb[jmb-1][iMaturity] = Flow/Q_tot[iMaturity][0]
        
        for count in range(1,Nscenarios):
            
            #Initializing the connectivity matrix including boundary conditions
            rhs = np.zeros((G.number_of_nodes(),1))
            rhs_x = np.zeros((G.number_of_nodes(),1)) #Initializing the right-hand side matrix of xylem pressure potentials
            rhs_p = np.zeros((G.number_of_nodes(),1)) #Initializing the right-hand side matrix of hydrostatic potentials for phloem BC
            rhs_e = np.zeros((G.number_of_nodes(),1)) #Initializing the right-hand side matrix of cell elongation
            rhs_o = np.zeros((G.number_of_nodes(),1)) #Initializing the right-hand side matrix of osmotic potentials
            Os_cells = np.zeros((Ncells,1)) #Initializing the cell osmotic potential vector
            Os_walls = np.zeros((Nwalls,1)) #Initializing the wall osmotic potential vector
            s_membranes = np.zeros((Nmb,1)) #Initializing the membrane reflection coefficient vector
            Os_membranes = np.zeros((Nmb,2)) #Initializing the osmotic potential storage side by side of membranes (0 for the wall, 1 for the protoplast)
            #rhs_s invariable between diferent scenarios but can vary for different hydraulic properties
            
            #Apoplastic & symplastic convective direction matrices initialization
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
            
            if C_flag:
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
                        for neighboor, eattr in edges.items(): #Loop on connections (edges)
                            if eattr['path'] == 'membrane': #Wall connection
                                if any(passage_cell_ID==array((indice[neighboor])-NwallsJun)):
                                    count_passage+=1
                                if G.node[neighboor]['cgroup']==3:#Endodermis
                                    count_endo+=1
                                elif G.node[neighboor]['cgroup']>4:#Pericycle or stele
                                    count_stele_overall+=1
                                elif G.node[neighboor]['cgroup']==4:#Cortex
                                    count_cortex+=1
                                elif G.node[neighboor]['cgroup']==1:#Exodermis
                                    count_exo+=1
                                elif G.node[neighboor]['cgroup']==2:#Epidermis
                                    count_epi+=1
                    for neighboor, eattr in edges.items(): #Loop on connections (edges)
                        j = (indice[neighboor]) #neighbouring node number
                        if j > i: #Only treating the information one way to save time
                            path = eattr['path'] #eattr is the edge attribute (i.e. connection type)
                            if path == "membrane": #Membrane connection
                                #Cell and wall osmotic potentials (cell types: 1=Exodermis;2=epidermis;3=endodermis;4=cortex;5=stele;16=pericycle)
                                rank=int(Cell_rank[int(j-NwallsJun)])
                                row=int(rank2row[rank])
                                if rank==1:#Exodermis
                                    Os_membranes[jmb][1]=Os_exo
                                    if count_epi==1: #wall between exodermis and epidermis
                                        s_membranes[jmb]=s_exo_epi
                                    elif count_epi==0: #wall between exodermis and cortex or between two exodermal cells
                                        s_membranes[jmb]=s_exo_cortex
                                elif rank==2:#Epidermis
                                    Os_membranes[jmb][1]=Os_epi
                                    s_membranes[jmb]=s_epi
                                elif rank==3:#Endodermis
                                    Os_membranes[jmb][1]=Os_endo
                                    if count_stele_overall==0: #wall between endodermis and cortex or between two endodermal cells
                                        s_membranes[jmb]=s_endo_cortex
                                    elif count_stele_overall>0 and count_endo>0: #wall between endodermis and pericycle
                                        s_membranes[jmb]=s_endo_peri
                                elif rank>=40 and rank<50:#Cortex
                                    if j-NwallsJun in InterCid:
                                        Os_membranes[jmb][1]=0
                                        s_membranes[jmb]=0
                                    else:
                                        if row==row_outercortex-7:
                                            Os_membranes[jmb][1]=Os_c8
                                        elif row==row_outercortex-6:
                                            Os_membranes[jmb][1]=Os_c7
                                        elif row==row_outercortex-5:
                                            Os_membranes[jmb][1]=Os_c6
                                        elif row==row_outercortex-4:
                                            Os_membranes[jmb][1]=Os_c5
                                        elif row==row_outercortex-3:
                                            Os_membranes[jmb][1]=Os_c4
                                        elif row==row_outercortex-2:
                                            Os_membranes[jmb][1]=Os_c3
                                        elif row==row_outercortex-1:
                                            Os_membranes[jmb][1]=Os_c2
                                        elif row==row_outercortex:
                                            Os_membranes[jmb][1]=Os_c1
                                        s_membranes[jmb]=s_cortex
                                elif G.node[j]['cgroup']==5:#Stelar parenchyma
                                    Os_membranes[jmb][1]=Os_stele
                                    s_membranes[jmb]=s_stele
                                elif rank==16:#Pericycle
                                    Os_membranes[jmb][1]=Os_peri
                                    s_membranes[jmb]=s_peri
                                elif G.node[j]['cgroup']==11 or G.node[j]['cgroup']==23:#Phloem sieve tube cell
                                    if not isnan(Os_sieve[0][count]):
                                        if Barrier>0 or j in listprotosieve:
                                            Os_membranes[jmb][1]=float(Os_sieve[0][count])
                                        else:
                                            Os_membranes[jmb][1]=Os_stele
                                    else:
                                        Os_membranes[jmb][1]=Os_stele
                                    s_membranes[jmb]=s_sieve
                                elif G.node[j]['cgroup']==12 or G.node[j]['cgroup']==26:#Companion cell
                                    if not isnan(Os_sieve[0][count]):
                                        Os_membranes[jmb][1]=Os_comp
                                    else:
                                        Os_membranes[jmb][1]=Os_stele
                                    s_membranes[jmb]=s_comp
                                elif G.node[j]['cgroup']==13 or G.node[j]['cgroup']==19 or G.node[j]['cgroup']==20:#Xylem cell or vessel
                                    if Barrier==0:
                                        Os_membranes[jmb][1]=Os_stele
                                        s_membranes[jmb]=s_stele
                                    else:
                                        Os_membranes[jmb][1]=0.0
                                        s_membranes[jmb]=0.0
                                jmb+=1
            
            #Soil and xylem water potentials
            #Psi_soil[0][count]=float(Psi_soil_range[count].get("pressure_left")) #Soil pressure potential (hPa)
            Psi_xyl[iMaturity][count]=float(BC_xyl_range[count].get("pressure")) #Xylem pressure potential (hPa)
            dPsi_xyl[iMaturity][count]=float(BC_xyl_range[count].get("deltaP")) #Xylem pressure potential change as compared to equilibrium pressure (hPa)
            Flow_xyl[0][count]=float(BC_xyl_range[count].get("flowrate")) #Xylem flow rate (cm^3/d)
            if not isnan(Flow_xyl[0][count]):
                if isnan(Psi_xyl[iMaturity][count]) and isnan(dPsi_xyl[iMaturity][count]):
                    tot_flow=Flow_xyl[0][count]
                    sum_area=0.0
                    i=1
                    for cid in listxyl:
                        area=cellarea[cid-NwallsJun]
                        Flow_xyl[i][count]=tot_flow*area
                        sum_area+=area
                        i+=1
                    i=1
                    for cid in listxyl:
                        Flow_xyl[i][count]/=sum_area #Total xylem flow rate partitioned proportionnally to xylem cross-section area
                        i+=1
                    if Flow_xyl[0][count]==0.0:
                        iEquil_xyl=count
                    if C_flag:
                        #Estimate the radial distribution of solutes later on from "u"
                        #First estimate water radial velocity in the apoplast
                        u=zeros((2,1))
                        u[0][0]=tot_flow/(height*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[0][0] #Cortex (cm/d)
                        u[1][0]=tot_flow/(height*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[1][0] #Stele (cm/d)
                else:
                    print('Error: Cannot have both pressure and flow BC at xylem boundary')
            elif not isnan(dPsi_xyl[iMaturity][count]):
                if isnan(Psi_xyl[iMaturity][count]):
                    Psi_xyl[iMaturity][count]=Psi_xyl[iMaturity][iEquil_xyl]+dPsi_xyl[iMaturity][count]
                else:
                    print('Error: Cannot have both pressure and pressure change relative to equilibrium as xylem boundary condition')
            if not isnan(Psi_xyl[iMaturity][count]):
                if C_flag:
                    #Estimate the radial distribution of solutes
                    #First estimate total flow rate (cm^3/d) from BC & kr
                    tot_flow1=0.0
                    u=zeros((2,1))
                    iter=0
                    tot_flow2=kr_tot[iMaturity][0]*perimeter*height*1.0E-04*(Psi_soil[0][count]+Os_soil[0][count]-Psi_xyl[iMaturity][count]-Os_xyl[0][count]) 
                    print('flow_rate =',tot_flow2,' iter =',iter)
                    #Convergence loop of water radial velocity and solute apoplastic convection-diffusion
                    while abs(tot_flow1-tot_flow2)/abs(tot_flow2)>0.001 and iter<30:
                        iter+=1
                        if iter==1:
                            tot_flow1=tot_flow2
                        elif iter>1 and sign(tot_flow1/tot_flow2)==1:
                            tot_flow1=(tot_flow1+tot_flow2)/2
                        else:
                            tot_flow1=tot_flow1/2
                        #Then estimate water radial velocity in the apoplast
                        u[0][0]=tot_flow1/(height*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[0][0] #Cortex apoplastic water velocity (cm/d) positive inwards
                        u[1][0]=tot_flow1/(height*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[1][0] #Stele apoplastic water velocity (cm/d) positive inwards
                        #Then estimate the radial solute distribution from an analytical solution (C(x)=C0+C0*(exp(u*x/D)-1)/(u/D*exp(u*x/D)-exp(u*L/D)+1)
                        Os_apo_cortex_eq=0.0
                        Os_apo_stele_eq=0.0
                        Os_sym_cortex_eq=0.0
                        Os_sym_stele_eq=0.0
                        #temp1=0.0
                        #temp2=0.0
                        jmb=0 #Index for membrane vector
                        for node, edges in G.adjacency_iter() : #adjacency_iter returns an iterator of (node, adjacency dict) tuples for all nodes. This is the fastest way to look at every edge. For directed graphs, only outgoing adjacencies are included.
                            i = indice[node] #Node ID number
                            if i<Nwalls: #wall ID
                                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                                    if eattr['path'] == 'membrane': #Wall connection
                                        if r_rel[i]>=0: #cortical side
                                            Os_apo=Os_soil[0][count]*exp(u[0][0]*abs(r_rel[i])*L_diff[0]/Os_soil[4][count])
                                            Os_apo_cortex_eq+=STFmb[jmb][iMaturity]*(Os_apo*s_membranes[jmb])
                                            Os_sym_cortex_eq+=STFmb[jmb][iMaturity]*(Os_membranes[jmb][1]*s_membranes[jmb])
                                            #temp1+=STFmb[jmb][iMaturity]
                                        else: #Stelar side
                                            Os_apo=Os_xyl[0][count]*exp(-u[1][0]*abs(r_rel[i])*L_diff[1]/Os_xyl[4][count])
                                            Os_apo_stele_eq-=STFmb[jmb][iMaturity]*(Os_apo*s_membranes[jmb])
                                            Os_sym_stele_eq-=STFmb[jmb][iMaturity]*(Os_membranes[jmb][1]*s_membranes[jmb])
                                            #temp2+=STFmb[jmb][iMaturity]
                                        Os_membranes[jmb][0]=Os_apo
                                        jmb+=1
                        tot_flow2=kr_tot[iMaturity][0]*perimeter*height*1.0E-04*(Psi_soil[0][count]+Os_apo_cortex_eq-Os_sym_cortex_eq-Psi_xyl[iMaturity][count]-Os_apo_stele_eq+Os_sym_stele_eq)
                        print('flow_rate =',tot_flow2,' iter =',iter)
                    u[0][0]=tot_flow2/(height*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[0][0] #Cortex (cm/d)
                    u[1][0]=tot_flow2/(height*1.0E-04)/(thickness*1.0E-04)/cell_per_layer[1][0] #Stele (cm/d)
                    ##Then estimate osmotic potentials in radial walls later on: C(x)=C0+C0*(exp(u*x/D)-1)/(u/D*exp(u*x/D)-exp(u*L/D)+1)
            
            #Elongation BC
            if Barrier==0: #No elongation from the Casparian strip on
                for wid in range(Nwalls):
                    rhs_e[wid][0]=lengths[wid]*thickness/2*1.0E-08*(Elong_cell[0][count]+(x_rel[wid]-0.5)*Elong_cell_side_diff[0][count])*Water_fraction_apo #cm^3/d Cell wall horizontal surface assumed to be rectangular (junctions are pointwise elements)
                for cid in range(Ncells):
                    if cellarea[cid]>cellperimeter[cid]*thickness/2:
                        rhs_e[NwallsJun+cid][0]=(cellarea[cid]-cellperimeter[cid]*thickness/2)*1.0E-8*(Elong_cell[0][count]+(x_rel[NwallsJun+cid]-0.5)*Elong_cell_side_diff[0][count])*Water_fraction_sym #cm^3/d Wall thickness removed from cell horizontal area to obtain protoplast horizontal area
                    else:
                        rhs_e[NwallsJun+cid][0]=0 #The cell elongation virtually does not imply water influx, though its walls do (typically intercellular spaces
                        #print('Cell too small to have '+str(thickness/2)+' micron thick walls')
                        #print('Cell ID & rank: '+str(cid)+' '+str(Cell_rank[cid]),'Total horizontal area: '+str(cellarea[cid])+' microns^2','Wall horizontal area: '+str(cellperimeter[cid]*thickness/2)+' microns^2')
            
            Psi_sieve[iMaturity][count]=float(BC_sieve_range[count].get("pressure")) #Phloem sieve element pressure potential (hPa)
            dPsi_sieve[iMaturity][count]=float(BC_sieve_range[count].get("deltaP")) #Phloem pressure potential change as compared to equilibrium pressure (hPa)
            Flow_sieve[0][count]=float(BC_sieve_range[count].get("flowrate")) #Phloem flow rate (cm^3/d)
            if not isnan(Flow_sieve[0][count]):
                if isnan(Psi_sieve[iMaturity][count]) and isnan(dPsi_sieve[iMaturity][count]):
                    if Barrier==0:
                        if Flow_sieve[0][count]==0:
                            tot_flow=-float(sum(rhs_e)) #"Equilibrium condition" with phloem water fully used for elongation
                        else:
                            tot_flow=Flow_sieve[0][count]
                        sum_area=0
                        i=1
                        for cid in listprotosieve:
                            area=cellarea[cid-NwallsJun]
                            Flow_sieve[i][count]=tot_flow*area
                            sum_area+=area
                            i+=1
                        i=1
                        for cid in listprotosieve:
                            Flow_sieve[i][count]/=sum_area #Total phloem flow rate partitioned proportionnally to phloem cross-section area
                            i+=1
                    elif Barrier>0:
                        tot_flow=Flow_sieve[0][count]
                        sum_area=0
                        i=1
                        for cid in listsieve:
                            area=cellarea[cid-NwallsJun]
                            Flow_sieve[i][count]=tot_flow*area
                            sum_area+=area
                            i+=1
                        i=1
                        for cid in listsieve:
                            Flow_sieve[i][count]/=sum_area #Total phloem flow rate partitioned proportionnally to phloem cross-section area
                            i+=1
                    if Flow_sieve[0][count]==0.0:
                        iEquil_sieve=count
                else:
                    print('Error: Cannot have both pressure and flow BC at phloem boundary')
            elif not isnan(dPsi_sieve[iMaturity][count]):
                if isnan(Psi_sieve[iMaturity][count]):
                    if not isnan(iEquil_sieve):
                        Psi_sieve[iMaturity][count]=Psi_sieve[iMaturity][iEquil_sieve]+dPsi_sieve[iMaturity][count]
                    else:
                        print('Error: Cannot have phloem pressure change relative to equilibrium without having a prior scenario with equilibrium phloem boundary condition')
                else:
                    print('Error: Cannot have both pressure and pressure change relative to equilibrium as phloem boundary condition')
            
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
                    if Os_soil[2][count] == 2: #Central symmetrical gradient for apoplastic osmotic potential
                        if Os_soil[4][count] == 0: #Not the analytical solution
                            Os_soil_local=float(Os_soil[0][count]+(Os_soil[1][count]-Os_soil[0][count])*abs(r_rel[i])**Os_soil[3][count])
                        else:
                            if r_rel[i]>=0: #cortical side
                                Os_soil_local=Os_soil[0][count]*exp(u[0][0]*abs(r_rel[i])*L_diff[0]/Os_soil[4][count])
                    elif Os_soil[2][count] == 1: #Left-right gradient for apoplastic osmotic potential
                        Os_soil_local=float(Os_soil[0][count]*(1-x_rel[i])+Os_soil[1][count]*x_rel[i])
                    if Os_xyl[2][count] == 2:
                        if Os_xyl[4][count] == 0: #Not the analytical solution
                            Os_xyl_local=float(Os_xyl[1][count]+(Os_xyl[0][count]-Os_xyl[1][count])*(1-abs(r_rel[i]))**Os_xyl[3][count])
                        else:
                            if r_rel[i]<0: #cortical side
                                Os_xyl_local=Os_xyl[0][count]*exp(-u[1][0]*abs(r_rel[i])*L_diff[1]/Os_xyl[4][count])
                    elif Os_xyl[2][count] == 1:
                        Os_xyl_local=float((Os_xyl[0][count]+Os_xyl[1][count])/2)
                    for neighboor, eattr in edges.items(): #Loop on connections (edges)
                        if eattr['path'] == 'membrane': #Wall connection
                            if any(passage_cell_ID==array((indice[neighboor])-NwallsJun)):
                                count_passage+=1
                            if G.node[neighboor]['cgroup']==3:#Endodermis
                                count_endo+=1
                            elif G.node[neighboor]['cgroup']>4:#Pericycle or stele
                                count_stele_overall+=1
                            elif G.node[neighboor]['cgroup']==4:#Cortex
                                count_cortex+=1
                            elif G.node[neighboor]['cgroup']==1:#Exodermis
                                count_exo+=1
                            elif G.node[neighboor]['cgroup']==2:#Epidermis
                                count_epi+=1
                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                    j = (indice[neighboor]) #neighbouring node number
                    if j > i: #Only treating the information one way to save time
                        path = eattr['path'] #eattr is the edge attribute (i.e. connection type)
                        if path == "membrane": #Membrane connection
                            #Cell and wall osmotic potentials (cell types: 1=Exodermis;2=epidermis;3=endodermis;4=cortex;5=stele;16=pericycle)
                            rank=int(Cell_rank[int(j-NwallsJun)])
                            row=int(rank2row[rank])
                            if rank==1:#Exodermis
                                Os_cells[j-NwallsJun]=Os_exo
                                Os_membranes[jmb][1]=Os_exo
                                OsCellLayer[row][iMaturity][count]+=Os_exo
                                nOsCellLayer[row][iMaturity][count]+=1
                                OsCellLayer[row+1][iMaturity][count]+=Os_exo
                                nOsCellLayer[row+1][iMaturity][count]+=1
                                Os_walls[i]=Os_soil_local
                                if count_epi==1: #wall between exodermis and epidermis
                                    s_membranes[jmb]=s_exo_epi
                                    OsWallLayer[row+1][iMaturity][count]+=Os_soil_local
                                    nOsWallLayer[row+1][iMaturity][count]+=1
                                elif count_epi==0: #wall between exodermis and cortex or between two exodermal cells
                                    s_membranes[jmb]=s_exo_cortex
                                    OsWallLayer[row][iMaturity][count]+=Os_soil_local
                                    nOsWallLayer[row][iMaturity][count]+=1
                            elif rank==2:#Epidermis
                                Os_cells[j-NwallsJun]=Os_epi
                                Os_membranes[jmb][1]=Os_epi
                                OsCellLayer[row][iMaturity][count]+=Os_epi
                                nOsCellLayer[row][iMaturity][count]+=1
                                Os_walls[i]=Os_soil_local
                                s_membranes[jmb]=s_epi
                                OsWallLayer[row][iMaturity][count]+=Os_soil_local
                                nOsWallLayer[row][iMaturity][count]+=1
                            elif rank==3:#Endodermis
                                Os_cells[j-NwallsJun]=Os_endo
                                Os_membranes[jmb][1]=Os_endo
                                OsCellLayer[row][iMaturity][count]+=Os_endo
                                nOsCellLayer[row][iMaturity][count]+=1
                                OsCellLayer[row+3][iMaturity][count]+=Os_endo
                                nOsCellLayer[row+3][iMaturity][count]+=1
                                if count_stele_overall==0 and count_cortex>0: #wall between endodermis and cortex or between two endodermal cells
                                    Os_walls[i]=Os_soil_local
                                    s_membranes[jmb]=s_endo_cortex
                                    #Not including the osmotic potential of walls that are located at the same place as the casparian strip
                                    OsWallLayer[row+3][iMaturity][count]+=Os_soil_local
                                    nOsWallLayer[row+3][iMaturity][count]+=1
                                elif count_stele_overall>0 and count_endo>0: #wall between endodermis and pericycle
                                    if Barrier==0: #No apoplastic barrier
                                        Os_walls[i]=Os_soil_local
                                        OsWallLayer[row][iMaturity][count]+=Os_soil_local
                                        nOsWallLayer[row][iMaturity][count]+=1
                                    else:
                                        Os_walls[i]=Os_xyl_local #float(Os_xyl[0][count])
                                        OsWallLayer[row][iMaturity][count]+=Os_xyl_local #float(Os_xyl[0][count])
                                        nOsWallLayer[row][iMaturity][count]+=1
                                    s_membranes[jmb]=s_endo_peri
                                else: #Wall between endodermal cells
                                    if Barrier==0: #No apoplastic barrier
                                        Os_walls[i]=Os_soil_local
                                        OsWallLayer[row][iMaturity][count]+=Os_soil_local
                                        nOsWallLayer[row][iMaturity][count]+=1
                                    else:
                                        Os_walls[i]=Os_xyl_local #float(Os_xyl[0][count])
                                        OsWallLayer[row][iMaturity][count]+=Os_xyl_local #float(Os_xyl[0][count])
                                        nOsWallLayer[row][iMaturity][count]+=1
                                    s_membranes[jmb]=s_endo_peri
                            elif rank>=40 and rank<50:#Cortex
                                if j-NwallsJun in InterCid: 
                                    Os_cells[j-NwallsJun]=Os_soil_local
                                    Os_walls[i]=Os_soil_local
                                    Os_membranes[jmb][1]=Os_soil_local
                                    Os_membranes[jmb][0]=Os_soil_local
                                    s_membranes[jmb]=0
                                    OsWallLayer[row][iMaturity][count]+=Os_soil_local
                                    nOsWallLayer[row][iMaturity][count]+=1
                                else:
                                    if row==row_outercortex-7:
                                        Os_cells[j-NwallsJun]=Os_c8
                                        Os_membranes[jmb][1]=Os_c8
                                        OsCellLayer[row][iMaturity][count]+=Os_c8
                                        nOsCellLayer[row][iMaturity][count]+=1
                                    elif row==row_outercortex-6:
                                        Os_cells[j-NwallsJun]=Os_c7
                                        Os_membranes[jmb][1]=Os_c7
                                        OsCellLayer[row][iMaturity][count]+=Os_c7
                                        nOsCellLayer[row][iMaturity][count]+=1
                                    elif row==row_outercortex-5:
                                        Os_cells[j-NwallsJun]=Os_c6
                                        Os_membranes[jmb][1]=Os_c6
                                        OsCellLayer[row][iMaturity][count]+=Os_c6
                                        nOsCellLayer[row][iMaturity][count]+=1
                                    elif row==row_outercortex-4:
                                        Os_cells[j-NwallsJun]=Os_c5
                                        Os_membranes[jmb][1]=Os_c5
                                        OsCellLayer[row][iMaturity][count]+=Os_c5
                                        nOsCellLayer[row][iMaturity][count]+=1
                                    elif row==row_outercortex-3:
                                        Os_cells[j-NwallsJun]=Os_c4
                                        Os_membranes[jmb][1]=Os_c4
                                        OsCellLayer[row][iMaturity][count]+=Os_c4
                                        nOsCellLayer[row][iMaturity][count]+=1
                                    elif row==row_outercortex-2:
                                        Os_cells[j-NwallsJun]=Os_c3
                                        Os_membranes[jmb][1]=Os_c3
                                        OsCellLayer[row][iMaturity][count]+=Os_c3
                                        nOsCellLayer[row][iMaturity][count]+=1
                                    elif row==row_outercortex-1:
                                        Os_cells[j-NwallsJun]=Os_c2
                                        Os_membranes[jmb][1]=Os_c2
                                        OsCellLayer[row][iMaturity][count]+=Os_c2
                                        nOsCellLayer[row][iMaturity][count]+=1
                                    elif row==row_outercortex:
                                        Os_cells[j-NwallsJun]=Os_c1
                                        Os_membranes[jmb][1]=Os_c1
                                        OsCellLayer[row][iMaturity][count]+=Os_c1
                                        nOsCellLayer[row][iMaturity][count]+=1
                                    Os_walls[i]=Os_soil_local
                                    s_membranes[jmb]=s_cortex
                                    OsWallLayer[row][iMaturity][count]+=Os_soil_local
                                    nOsWallLayer[row][iMaturity][count]+=1
                            elif G.node[j]['cgroup']==5:#Stelar parenchyma
                                Os_cells[j-NwallsJun]=Os_stele
                                Os_membranes[jmb][1]=Os_stele
                                OsCellLayer[row][iMaturity][count]+=Os_stele
                                nOsCellLayer[row][iMaturity][count]+=1
                                if Barrier==0: #No apoplastic barrier
                                    Os_walls[i]=Os_soil_local
                                    OsWallLayer[row][iMaturity][count]+=Os_soil_local
                                    nOsWallLayer[row][iMaturity][count]+=1
                                else:
                                    Os_walls[i]=Os_xyl_local #float(Os_xyl[0][count])
                                    OsWallLayer[row][iMaturity][count]+=Os_xyl_local #float(Os_xyl[0][count])
                                    nOsWallLayer[row][iMaturity][count]+=1
                                s_membranes[jmb]=s_stele
                            elif rank==16:#Pericycle
                                Os_cells[j-NwallsJun]=Os_peri
                                Os_membranes[jmb][1]=Os_peri
                                OsCellLayer[row][iMaturity][count]+=Os_peri
                                nOsCellLayer[row][iMaturity][count]+=1
                                if Barrier==0: #No apoplastic barrier
                                    Os_walls[i]=Os_soil_local
                                    OsWallLayer[row][iMaturity][count]+=Os_soil_local
                                    nOsWallLayer[row][iMaturity][count]+=1
                                else:
                                    Os_walls[i]=Os_xyl_local #float(Os_xyl[0][count])
                                    OsWallLayer[row][iMaturity][count]+=Os_xyl_local #float(Os_xyl[0][count])
                                    nOsWallLayer[row][iMaturity][count]+=1
                                s_membranes[jmb]=s_peri
                            elif G.node[j]['cgroup']==11 or G.node[j]['cgroup']==23:#Phloem sieve tube cell
                                if not isnan(Os_sieve[0][count]):
                                    if Barrier>0 or j in listprotosieve:
                                        Os_cells[j-NwallsJun]=float(Os_sieve[0][count])
                                        Os_membranes[jmb][1]=float(Os_sieve[0][count])
                                        OsCellLayer[row][iMaturity][count]+=float(Os_sieve[0][count])
                                        nOsCellLayer[row][iMaturity][count]+=1
                                    else:
                                        Os_cells[j-NwallsJun]=Os_stele
                                        Os_membranes[jmb][1]=Os_stele
                                        OsCellLayer[row][iMaturity][count]+=Os_stele
                                        nOsCellLayer[row][iMaturity][count]+=1
                                else:
                                    Os_cells[j-NwallsJun]=Os_stele
                                    Os_membranes[jmb][1]=Os_stele
                                    OsCellLayer[row][iMaturity][count]+=Os_stele
                                    nOsCellLayer[row][iMaturity][count]+=1
                                if Barrier==0: #No apoplastic barrier
                                    Os_walls[i]=Os_soil_local
                                    OsWallLayer[row][iMaturity][count]+=Os_soil_local
                                    nOsWallLayer[row][iMaturity][count]+=1
                                else:
                                    Os_walls[i]=Os_xyl_local #float(Os_xyl[0][count])
                                    OsWallLayer[row][iMaturity][count]+=Os_xyl_local #float(Os_xyl[0][count])
                                    nOsWallLayer[row][iMaturity][count]+=1
                                s_membranes[jmb]=s_sieve
                            elif G.node[j]['cgroup']==12 or G.node[j]['cgroup']==26:#Companion cell
                                if not isnan(Os_sieve[0][count]):
                                    Os_cells[j-NwallsJun]=Os_comp
                                    Os_membranes[jmb][1]=Os_comp
                                    OsCellLayer[row][iMaturity][count]+=Os_comp
                                    nOsCellLayer[row][iMaturity][count]+=1
                                else:
                                    Os_cells[j-NwallsJun]=Os_stele
                                    Os_membranes[jmb][1]=Os_stele
                                    OsCellLayer[row][iMaturity][count]+=Os_stele
                                    nOsCellLayer[row][iMaturity][count]+=1
                                if Barrier==0: #No apoplastic barrier
                                    Os_walls[i]=Os_soil_local
                                    OsWallLayer[row][iMaturity][count]+=Os_soil_local
                                    nOsWallLayer[row][iMaturity][count]+=1
                                else:
                                    Os_walls[i]=Os_xyl_local
                                    OsWallLayer[row][iMaturity][count]+=Os_xyl_local
                                    nOsWallLayer[row][iMaturity][count]+=1
                                s_membranes[jmb]=s_comp
                            elif G.node[j]['cgroup']==13 or G.node[j]['cgroup']==19 or G.node[j]['cgroup']==20:#Xylem cell or vessel
                                if Barrier==0:
                                    Os_cells[j-NwallsJun]=Os_stele
                                    Os_membranes[jmb][1]=Os_stele
                                    OsCellLayer[row][iMaturity][count]+=Os_stele
                                    nOsCellLayer[row][iMaturity][count]+=1
                                    Os_walls[i]=Os_soil_local
                                    s_membranes[jmb]=s_stele
                                    OsWallLayer[row][iMaturity][count]+=Os_soil_local
                                    nOsWallLayer[row][iMaturity][count]+=1
                                else:
                                    Os_cells[j-NwallsJun]=Os_xyl_local
                                    Os_membranes[jmb][0]=Os_xyl_local
                                    Os_membranes[jmb][1]=Os_xyl_local
                                    Os_membranes[jmb][1]=Os_xyl_local
                                    Os_walls[i]=Os_xyl_local
                                    s_membranes[jmb]=0.0
                                    OsWallLayer[row][iMaturity][count]+=Os_xyl_local #float(Os_xyl[0][count])
                                    nOsWallLayer[row][iMaturity][count]+=1
                            K=Kmb[jmb]
                            rhs_o[i]+= K*s_membranes[jmb]*(Os_walls[i] - Os_cells[j-NwallsJun]) #Wall node
                            rhs_o[j]+= K*s_membranes[jmb]*(Os_cells[j-NwallsJun] - Os_walls[i]) #Cell node 
                            jmb+=1
            for row in range(int(r_discret[0])):
                if nOsWallLayer[row][iMaturity][count]>0:
                    OsWallLayer[row][iMaturity][count]=OsWallLayer[row][iMaturity][count]/nOsWallLayer[row][iMaturity][count]
                if nOsCellLayer[row][iMaturity][count]>0:
                    OsCellLayer[row][iMaturity][count]=OsCellLayer[row][iMaturity][count]/nOsCellLayer[row][iMaturity][count]
            
            #Xylem BC
            if Barrier>0: #No mature xylem before the Casparian strip stage
                if not isnan(Psi_xyl[iMaturity][count]): #Pressure xylem BC
                    for cid in listxyl:
                        rhs_x[cid][0] = -K_xyl  #Axial conductance of xylem vessels
                        matrix_W[cid][cid] -= K_xyl
                elif not isnan(Flow_xyl[0][count]): #Flow xylem BC
                    i=1
                    for cid in listxyl:
                        rhs_x[cid][0] = Flow_xyl[i][count] #(cm^3/d)
                        i+=1
            
            #Phloem BC
            if Barrier==0: #Protophloem only
                if not isnan(Psi_sieve[iMaturity][count]):
                    for cid in listprotosieve:
                        rhs_p[cid][0] = -K_sieve  #Axial conductance of phloem sieve tube
                        matrix_W[cid][cid] -= K_sieve
                elif not isnan(Flow_sieve[0][count]):
                    i=1
                    for cid in listprotosieve:
                        rhs_p[cid][0] = Flow_sieve[i][count] #(cm^3/d)
                        i+=1
            elif Barrier>0: #Includes mature phloem
                if not isnan(Psi_sieve[iMaturity][count]): #Then there is a phloem BC in scenarios (assuming that we did not pick scenarios with and others without)
                    for cid in listsieve: #both proto and metaphloem
                        rhs_p[cid][0] = -K_sieve  #Axial conductance of xylem vessels
                        matrix_W[cid][cid] -= K_sieve
                elif not isnan(Flow_sieve[0][count]):
                    i=1
                    for cid in listsieve:
                        rhs_p[cid][0] = Flow_sieve[i][count] #(cm^3/d)
                        i+=1
            
            
            #Adding up all BC
            #Elongation BC
            rhs += rhs_e
            
            #Osmotic BC
            rhs += rhs_o
            
            #Soil BC
            rhs += np.multiply(rhs_s,Psi_soil[0][count]*(1-x_rel)+Psi_soil[1][count]*x_rel)
            
            #Xylem BC
            if not isnan(Psi_xyl[iMaturity][count]): #Pressure xylem BC
                rhs += rhs_x*Psi_xyl[iMaturity][count]  #multiplication of rhs components delayed till this point so that rhs_s & rhs_x can be re-used
            elif not isnan(Flow_xyl[0][count]): #Flow xylem BC
                rhs += rhs_x
            
            #Phloem BC
            if not isnan(Flow_sieve[0][count]):
                rhs += rhs_p
            elif not isnan(Psi_sieve[iMaturity][count]):
                rhs += rhs_p*Psi_sieve[iMaturity][count]
            
            ##################################################
            ##Solve Doussan equation, results in soln matrix##
            ##################################################
            
            soln = np.linalg.solve(matrix_W,rhs) #Solving the equation to get potentials inside the network
            
            #Verification that computation was correct
            verif1=np.allclose(np.dot(matrix_W,soln),rhs)
            
            #print("Correct computation on PSI ?", verif1)
            
            #Removing Xylem and phloem BC terms in "matrix" in case they would change in the next scenario
            if Barrier>0:
                if not isnan(Psi_xyl[iMaturity][count]): #Pressure xylem BC
                    for cid in listxyl:
                        matrix_W[cid][cid] += K_xyl
            if Barrier==0: #Protophloem only
                if not isnan(Psi_sieve[iMaturity][count]):
                    for cid in listprotosieve:
                        matrix_W[cid][cid] += K_sieve
            elif Barrier>0: #Includes mature phloem
                if not isnan(Psi_sieve[iMaturity][count]): #Then there is a phloem BC in scenarios (assuming that we did not pick scenarios with and others without)
                    for cid in listsieve: #both proto and metaphloem
                        matrix_W[cid][cid] += K_sieve
            
            #Flow rates at interfaces
            Q_soil=[]
            for ind in Borderwall:
                Q=rhs_s[ind]*(soln[ind]-(Psi_soil[0][count]*(1-x_rel[ind])+Psi_soil[1][count]*x_rel[ind]))
                Q_soil.append(Q) #(cm^3/d) Positive for water flowing into the root, rhs_s is minus the conductance at the soil root interface
                if Apo_Contagion==2:
                    if Sym_Contagion==2:
                        if ind not in Apo_w_Zombies0:
                            if Q<0.0:
                                matrix_C[ind][ind]+=Q
                    else:
                        if ind not in Apo_w_Zombies0:
                            if Q<0.0:
                                matrix_ApoC[ind][ind]+=Q
                #if C_flag and Os_soil[5][count]==1:
                #if Apo_Contagion==2:
                #    #if not Q==0:
                #    #    list_walls_apo_conv.append(ind)
                #    if Q>0:
                #        rhs_ApoC[ind][0] -= Q*Os_soil[0][count] #Flow rate to be multiplied by concentration BC
                #    else:
                #        matrix_ApoC[ind][ind] += Q
                        
            for ind in Borderjunction:
                Q=rhs_s[ind]*(soln[ind]-(Psi_soil[0][count]*(1-x_rel[ind])+Psi_soil[1][count]*x_rel[ind]))
                Q_soil.append(Q) #(cm^3/d) Positive for water flowing into the root
                if Apo_Contagion==2:
                    if Sym_Contagion==2:
                        if ind not in Apo_j_Zombies0:
                            if Q<0.0:
                                matrix_C[ind][ind]+=Q
                    else:
                        if ind not in Apo_j_Zombies0:
                            if Q<0.0:
                                matrix_ApoC[ind][ind]+=Q
                #if C_flag and Os_soil[5][count]==1:
                #    if not Q==0:
                #        list_walls_apo_conv.append(ind)
                #    if Q>0:
                #        rhs_C[ind][0] -= Q*Os_soil[0][count] #Flow rate times concentration BC
                #    else:
                #        matrix_C[ind][ind] += Q
            
            Q_xyl=[]
            if Barrier>0:
                if not isnan(Psi_xyl[iMaturity][count]): #Xylem pressure BC
                    for cid in listxyl:
                        Q=rhs_x[cid]*(soln[cid]-Psi_xyl[iMaturity][count])
                        Q_xyl.append(Q) #(cm^3/d) Negative for water flowing into xylem tubes
                        rank=int(Cell_rank[cid-NwallsJun])
                        row=int(rank2row[rank])
                        Q_xyl_layer[row][iMaturity][count] += Q
                        #if C_flag:
                        #    if Q>0: #Water leaving the cross-section
                        #        matrix_C[cid][cid] -= Q
                        #    else: #Water entering the cross-section through xylem
                        #        rhs_C[cid][0] += Q #Flow rate times concentration BC
                        #    rhs_C[cid][0] *= Os_xyl[1][count]
                elif not isnan(Flow_xyl[0][count]): #Xylem flow BC
                    for cid in listxyl:
                        Q=-rhs_x[cid]
                        Q_xyl.append(Q) #(cm^3/d) Negative for water flowing into xylem tubes
                        rank=int(Cell_rank[cid-NwallsJun])
                        row=int(rank2row[rank])
                        Q_xyl_layer[row][iMaturity][count] += Q
                        #if C_flag:
                        #    if Q>0: #Water leaving the cross-section
                        #        matrix_C[cid][cid] -= Q
                        #    else: #Water entering the cross-section through xylem
                        #        rhs_C[cid][0] += Q #Flow rate times concentration BC
                        #    rhs_C[cid][0] *= Os_xyl[1][count]
            
            Q_sieve=[]
            if Barrier==0:
                if not isnan(Psi_sieve[iMaturity][count]): #Phloem pressure BC
                    for cid in listprotosieve: #Q will be 0 for metaphloem if Barrier==0 because rhs_p=0 for these cells
                        Q=rhs_p[cid]*(soln[cid]-Psi_sieve[iMaturity][count])
                        Q_sieve.append(Q) #(cm^3/d) Positive for water flowing from sieve tubes
                        rank=int(Cell_rank[cid-NwallsJun])
                        row=int(rank2row[rank])
                        Q_sieve_layer[row][iMaturity][count] += Q
                elif not isnan(Flow_sieve[0][count]): #Phloem flow BC
                    for cid in listprotosieve:
                        Q=-rhs_p[cid]
                        Q_sieve.append(Q) #(cm^3/d) Negative for water flowing into xylem tubes
                        rank=int(Cell_rank[cid-NwallsJun])
                        row=int(rank2row[rank])
                        Q_sieve_layer[row][iMaturity][count] += Q
            elif Barrier>0:
                if not isnan(Psi_sieve[iMaturity][count]): #Phloem pressure BC
                    for cid in listsieve: #Q will be 0 for metaphloem if Barrier==0 because rhs_p=0 for these cells
                        Q=rhs_p[cid]*(soln[cid]-Psi_sieve[iMaturity][count])
                        Q_sieve.append(Q) #(cm^3/d) Positive for water flowing from sieve tubes
                        rank=int(Cell_rank[cid-NwallsJun])
                        row=int(rank2row[rank])
                        Q_sieve_layer[row][iMaturity][count] += Q
                elif not isnan(Flow_sieve[0][count]): #Phloem flow BC
                    for cid in listsieve:
                        Q=-rhs_p[cid]
                        Q_sieve.append(Q) #(cm^3/d) Negative for water flowing into xylem tubes
                        rank=int(Cell_rank[cid-NwallsJun])
                        row=int(rank2row[rank])
                        Q_sieve_layer[row][iMaturity][count] += Q
            Q_elong=-rhs_e #(cm^3/d) The elongation flux virtually disappears from the cross-section => negative
            for cid in range(Ncells):
                rank=int(Cell_rank[cid])
                row=int(rank2row[rank])
                Q_elong_layer[row][iMaturity][count] += Q_elong[NwallsJun+cid]
            Q_tot[iMaturity][count]=sum(Q_soil) #(cm^3/d) Total flow rate at root surface
            for ind in range(NwallsJun,len(G.node)): #NwallsJun is the index of the first cell
                cellnumber1=ind-NwallsJun
                rank = int(Cell_rank[cellnumber1])
                row = int(rank2row[rank])
                if rank == 1: #Exodermis
                    PsiCellLayer[row][iMaturity][count] += soln[ind]*(STFcell_plus[cellnumber1][iMaturity]+abs(STFcell_minus[cellnumber1][iMaturity]))/(STFlayer_plus[row][iMaturity]+abs(STFlayer_minus[row][iMaturity])+STFlayer_plus[row+1][iMaturity]+abs(STFlayer_minus[row+1][iMaturity])) #(hPa)
                    PsiCellLayer[row+1][iMaturity][count] += soln[ind]*(STFcell_plus[cellnumber1][iMaturity]+abs(STFcell_minus[cellnumber1][iMaturity]))/(STFlayer_plus[row][iMaturity]+abs(STFlayer_minus[row][iMaturity])+STFlayer_plus[row+1][iMaturity]+abs(STFlayer_minus[row+1][iMaturity])) #(hPa)
                elif rank == 3: #Endodermis
                    if any(passage_cell_ID==array(cellnumber1)) and Barrier==2: #Passage cell
                        PsiCellLayer[row+1][iMaturity][count] += soln[ind]*(STFcell_plus[cellnumber1][iMaturity]+abs(STFcell_minus[cellnumber1][iMaturity]))/(STFlayer_plus[row+1][iMaturity]+abs(STFlayer_minus[row+1][iMaturity])+STFlayer_plus[row+2][iMaturity]+abs(STFlayer_minus[row+2][iMaturity])) #(hPa)
                        PsiCellLayer[row+2][iMaturity][count] += soln[ind]*(STFcell_plus[cellnumber1][iMaturity]+abs(STFcell_minus[cellnumber1][iMaturity]))/(STFlayer_plus[row+1][iMaturity]+abs(STFlayer_minus[row+1][iMaturity])+STFlayer_plus[row+2][iMaturity]+abs(STFlayer_minus[row+2][iMaturity])) #(hPa)
                    else:
                        PsiCellLayer[row][iMaturity][count] += soln[ind]*(STFcell_plus[cellnumber1][iMaturity]+abs(STFcell_minus[cellnumber1][iMaturity]))/(STFlayer_plus[row][iMaturity]+abs(STFlayer_minus[row][iMaturity])+STFlayer_plus[row+3][iMaturity]+abs(STFlayer_minus[row+3][iMaturity])) #(hPa)
                        PsiCellLayer[row+3][iMaturity][count] += soln[ind]*(STFcell_plus[cellnumber1][iMaturity]+abs(STFcell_minus[cellnumber1][iMaturity]))/(STFlayer_plus[row][iMaturity]+abs(STFlayer_minus[row][iMaturity])+STFlayer_plus[row+3][iMaturity]+abs(STFlayer_minus[row+3][iMaturity])) #(hPa)
                        if not Barrier==2:
                            PsiCellLayer[row+1][iMaturity][count] = nan
                            PsiCellLayer[row+2][iMaturity][count] = nan
                elif (ind not in listxyl) or Barrier==0: #Not for mature xylem
                    PsiCellLayer[row][iMaturity][count] += soln[ind]*(STFcell_plus[cellnumber1][iMaturity]+abs(STFcell_minus[cellnumber1][iMaturity]))/(STFlayer_plus[row][iMaturity]+abs(STFlayer_minus[row][iMaturity])) #(hPa)
            
            if Barrier>0 and isnan(Psi_xyl[iMaturity][count]):
                Psi_xyl[iMaturity][count]=0.0
                for cid in listxyl:
                    Psi_xyl[iMaturity][count]+=soln[cid]/Nxyl
            if Barrier>0:
                if isnan(Psi_sieve[iMaturity][count]):
                    Psi_sieve[iMaturity][count]=0.0
                    for cid in listsieve:
                        Psi_sieve[iMaturity][count]+=soln[cid]/Nsieve #Average of phloem water pressures
            elif Barrier==0:
                if isnan(Psi_sieve[iMaturity][count]):
                    Psi_sieve[iMaturity][count]=0.0
                    for cid in listprotosieve:
                        Psi_sieve[iMaturity][count]+=soln[cid]/Nprotosieve #Average of protophloem water pressures
            
            print("Uptake rate per unit root length: soil ",(sum(Q_soil)/height/1.0E-04),"cm^2/d, xylem ",(sum(Q_xyl)/height/1.0E-04),"cm^2/d, phloem ",(sum(Q_sieve)/height/1.0E-04),"cm^2/d, elongation ",(sum(Q_elong)/height/1.0E-04),"cm^2/d")
            if not isnan(sum(Q_sieve)):
                print("Mass balance error:",(sum(Q_soil)+sum(Q_xyl)+sum(Q_sieve)+sum(Q_elong))/height/1.0E-04,"cm^2/d")
            else:
                print("Mass balance error:",(sum(Q_soil)+sum(Q_xyl)+sum(Q_elong))/height/1.0E-04,"cm^2/d")
            
            #################################################################
            ##Calul of Fluxes between nodes and Creating the edge_flux_list##
            #################################################################
            
            #Creating a list for the fluxes
            #edge_flux_list=[]
            
            #Filling the fluxes list
            MembraneFlowDensity=[]
            WallFlowDensity=[]
            WallFlowDensity_cos=[]
            PlasmodesmFlowDensity=[]
            Fjw_list=[]
            Fcw_list=[]
            Fcc_list=[]
            jmb=0 #Index for membrane conductance vector
            for node, edges in G.adjacency_iter() : #adjacency_iter returns an iterator of (node, adjacency dict) tuples for all nodes. This is the fastest way to look at every edge. For directed graphs, only outgoing adjacencies are included.
                i = indice[node] #Node ID number
                psi = soln[i] #Node water potential
                psi_o_cell = inf #Opposite cell water potential
                ind_o_cell = inf #Opposite cell index
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
                noPD=False #Initializes the flag for wall connected to an intercellular space -> does not have plasmodesmata
                if i<Nwalls: #wall ID
                    for neighboor, eattr in edges.items(): #Loop on connections (edges)
                        if eattr['path'] == 'membrane': #Wall connection
                            if any(passage_cell_ID==array((indice[neighboor])-NwallsJun)):
                                count_passage+=1
                            if any(InterCid==array((indice[neighboor])-NwallsJun)):
                                count_interC+=1
                            if G.node[neighboor]['cgroup']==3:#Endodermis
                                count_endo+=1
                            elif G.node[neighboor]['cgroup']==13 or G.node[neighboor]['cgroup']==19 or G.node[neighboor]['cgroup']==20:#Xylem cell or vessel
                                count_xyl+=1
                            elif G.node[neighboor]['cgroup']==16 or G.node[neighboor]['cgroup']==21:#Pericycle or stele
                                count_peri+=1
                                if neighboor in PPP:
                                    count_PPP+=1
                            elif G.node[neighboor]['cgroup']==1:#Exodermis
                                count_exo+=1
                            elif G.node[neighboor]['cgroup']==2:#Epidermis
                                count_epi+=1
                            elif G.node[neighboor]['cgroup']==4:#Cortex
                                count_cortex+=1
                            elif G.node[neighboor]['cgroup']==5:#Stelar parenchyma
                                count_stele+=1
                            elif G.node[neighboor]['cgroup']==11 or G.node[neighboor]['cgroup']==23:#Phloem sieve tube
                                count_sieve+=1
                            elif G.node[neighboor]['cgroup']==12 or G.node[neighboor]['cgroup']==26:#Companion cell
                                count_comp+=1
                            if G.node[neighboor]['cgroup']>4:#Stele overall
                                count_stele_overall+=1
                ijunction=0
                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                    j = indice[neighboor] #Neighbouring node ID number
                    #if j > i: #Only treating the information one way to save time
                    psin = soln[j] #Neighbouring node water potential
                    path = eattr['path'] #eattr is the edge attribute (i.e. connection type)
                    if i<Nwalls:
                        if Paraview==1 or ParTrack==1 or Apo_Contagion>0 or Sym_Contagion>0:
                            if path == "wall":
                                #K = eattr['kw']*1.0E-04*((eattr['lat_dist']+height)*eattr['thickness']-square(eattr['thickness']))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)
                                if (count_interC>=2 and Barrier>0) or (count_xyl==2 and Xylem_pieces): #"Fake wall" splitting an intercellular space or a xylem cell in two
                                    K = 1.0E-16 #Non conductive
                                elif count_cortex>=2: #wall between two cortical cells
                                    K = kw_cortex_cortex*1.0E-04*((eattr['lat_dist']+height)*thickness-square(thickness))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)
                                elif count_endo>=2: #wall between two endodermis cells
                                    K = kw_endo_endo*1.0E-04*((eattr['lat_dist']+height)*thickness-square(thickness))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)
                                elif count_stele_overall>0 and count_endo>0: #wall between endodermis and pericycle
                                    if count_passage>0:
                                        K = kw_passage*1.0E-04*((eattr['lat_dist']+height)*thickness-square(thickness))/eattr['length']
                                    else:
                                        K = kw_endo_peri*1.0E-04*((eattr['lat_dist']+height)*thickness-square(thickness))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)
                                elif count_stele_overall==0 and count_endo==1: #wall between endodermis and cortex
                                    if count_passage>0:
                                        K = kw_passage*1.0E-04*((eattr['lat_dist']+height)*thickness-square(thickness))/eattr['length']
                                    else:
                                        K = kw_endo_cortex*1.0E-04*((eattr['lat_dist']+height)*thickness-square(thickness))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)
                                elif count_exo>=2: #wall between two exodermis cells
                                    K = kw_exo_exo*1.0E-04*((eattr['lat_dist']+height)*thickness-square(thickness))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)
                                else: #other walls
                                    K = kw*1.0E-04*((eattr['lat_dist']+height)*thickness-square(thickness))/eattr['length'] #Junction-Wall conductance (cm^3/hPa/d)
                                Fjw = K * (psin - psi) * sign(j-i) #(cm^3/d) Water flow rate positive from junction to wall
                                Fjw_list.append((i,j,Fjw))
                                #The ordering in WallFlowDensity will correspond to the one of ThickWallsX, saved for display only
                                WallFlowDensity.append((i,j, Fjw / (((eattr['lat_dist']+height)*thickness-square(thickness))*1.0E-08))) # (cm/d) Positive towards lower node ID 
                                cos_angle=(position[i][0]-position[j][0])/(hypot(position[j][0]-position[i][0],position[j][1]-position[i][1])) #Vectors junction1-wall
                                WallFlowDensity_cos.append((i,j, cos_angle * Fjw / (((eattr['lat_dist']+height)*thickness-square(thickness))*1.0E-08))) # (cm/d) Positive towards lower node ID 
                                #if C_flag and Os_soil[5][count]*Os_xyl[5][count]==1:
                                if Apo_Contagion==2:
                                    if Sym_Contagion==2: # Apo & Sym contagion
                                        if Fjw>0: #Flow from junction to wall
                                            if i not in Apo_w_Zombies0:
                                                matrix_C[i][j] += Fjw
                                            if j not in Apo_j_Zombies0:
                                                matrix_C[j][j] -= Fjw
                                        else: #Flow from wall to junction
                                            if i not in Apo_w_Zombies0:
                                                matrix_C[i][i] += Fjw
                                            if j not in Apo_j_Zombies0:
                                                matrix_C[j][i] -= Fjw
                                    else: #Only Apo contagion
                                        if Fjw>0: #Flow from junction to wall
                                            if i not in Apo_w_Zombies0:
                                                matrix_ApoC[i][j] += Fjw
                                            if j not in Apo_j_Zombies0:
                                                matrix_ApoC[j][j] -= Fjw
                                        else: #Flow from wall to junction
                                            if i not in Apo_w_Zombies0:
                                                matrix_ApoC[i][i] += Fjw
                                            if j not in Apo_j_Zombies0:
                                                matrix_ApoC[j][i] -= Fjw
                                
                                if Apo_Contagion==1:
                                    if Fjw>0:
                                        Apo_connec_flow[j][nApo_connec_flow[j]]=i
                                        nApo_connec_flow[j]+=1
                                    elif Fjw<0:
                                        Apo_connec_flow[i][nApo_connec_flow[i]]=j
                                        nApo_connec_flow[i]+=1
                            elif path == "membrane": #Membrane connection
                                #K = (eattr['kmb']+eattr['kaqp'])*1.0E-08*(height+eattr['dist'])*eattr['length']
                                if G.node[j]['cgroup']==1: #Exodermis
                                    kaqp=kaqp_exo
                                elif G.node[j]['cgroup']==2: #Epidermis
                                    kaqp=kaqp_epi
                                elif G.node[j]['cgroup']==3: #Endodermis
                                    kaqp=kaqp_endo
                                elif G.node[j]['cgroup']==13 or G.node[j]['cgroup']==19 or G.node[j]['cgroup']==20: #xylem cell or vessel
                                    if Barrier>0: #Xylem vessel
                                        kaqp=kaqp_stele*10000 #No membrane resistance because no membrane
                                        noPD=True
                                    elif Barrier==0: #Xylem cell
                                        kaqp=kaqp_stele
                                        if (count_xyl==2 and Xylem_pieces):
                                            noPD=True
                                elif G.node[j]['cgroup']>4: #Stele and pericycle
                                    kaqp=kaqp_stele
                                elif (j-NwallsJun in InterCid) and Barrier>0: #the neighbour is an intercellular space "cell"
                                    kaqp=kInterC
                                    noPD=True
                                elif G.node[j]['cgroup']==4: #Cortex
                                    kaqp=float(a_cortex*dist_grav[wid]*1.0E-04+b_cortex) #AQP activity (cm/hPa/d)
                                    if kaqp < 0:
                                        error('Error, negative kaqp in cortical cell, adjust Paqp_cortex')
                                #Calculating conductances
                                if count_endo>=2: #wall between two endodermis cells, in this case the suberized wall can limit the transfer of water between cell and wall
                                    if kw_endo_endo==0.00:
                                        K=0.00
                                    else:
                                        K = 1/(1/(kw_endo_endo/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(height+eattr['dist'])*eattr['length'] #(cm^3/hPa/d)
                                elif count_exo>=2: #wall between two exodermis cells, in this case the suberized wall can limit the transfer of water between cell and wall
                                    if kw_exo_exo==0.00:
                                        K=0.00
                                    else:
                                        K = 1/(1/(kw_exo_exo/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(height+eattr['dist'])*eattr['length'] #(cm^3/hPa/d)
                                elif count_stele_overall>0 and count_endo>0: #wall between endodermis and pericycle, in this case the suberized wall can limit the transfer of water between cell and wall
                                    if count_passage>0:
                                        K = 1/(1/(kw_passage/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(height+eattr['dist'])*eattr['length']
                                    else:
                                        if kw_endo_peri==0.00:
                                            K=0.00
                                        else:
                                            K = 1/(1/(kw_endo_peri/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(height+eattr['dist'])*eattr['length']
                                elif count_stele_overall==0 and count_endo==1: #wall between endodermis and cortex, in this case the suberized wall can limit the transfer of water between cell and wall
                                    if kaqp==0.0:
                                        K=1.00E-16
                                    else:
                                        if count_passage>0:
                                            K = 1/(1/(kw_passage/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(height+eattr['dist'])*eattr['length']
                                        else:
                                            if kw_endo_cortex==0.00:
                                                K=0.00
                                            else:
                                                K = 1/(1/(kw_endo_cortex/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(height+eattr['dist'])*eattr['length']
                                else:
                                    if kaqp==0.0:
                                        K=1.00E-16
                                    else:
                                        K = 1/(1/(kw/(thickness/2*1.0E-04))+1/(kmb+kaqp))*1.0E-08*(height+eattr['dist'])*eattr['length'] #(cm^3/hPa/d)
                                Fcw = K * (psi - psin + s_membranes[jmb]*(Os_walls[i] - Os_cells[j-NwallsJun])) #(cm^3/d) Water flow rate positive from wall to protoplast
                                Fcw_list.append((i,j,-Fcw,s_membranes[jmb])) #Water flow rate positive from protoplast to wall
                                #Flow densities calculation
                                #The ordering in MembraneFlowDensity will correspond to the one of ThickWalls, saved for display only 
                                MembraneFlowDensity.append(Fcw / (1.0E-08*(height+eattr['dist'])*eattr['length']))
                                ####Solute convection across membranes####
                                if Apo_Contagion==2 and Sym_Contagion==2:
                                    if Fcw>0: #Flow from wall to protoplast
                                        if i not in Apo_w_Zombies0:
                                            if D2O1==1:#Solute that moves across membranes like water 
                                                matrix_C[i][i] -= Fcw
                                            else: #Solute that moves across membranes independently of water (the membrane is possibly not one) 
                                                matrix_C[i][i] -= Fcw*(1-s_membranes[jmb])
                                        if j-NwallsJun not in Sym_Zombie0:
                                            if D2O1==1:#Solute that moves across membranes like water 
                                                matrix_C[j][i] += Fcw
                                            else: #Solute that moves across membranes independently of water (the membrane is possibly not one) 
                                                matrix_C[j][i] += Fcw*(1-s_membranes[jmb])
                                    else: #Flow from protoplast to wall
                                        if j-NwallsJun not in Sym_Zombie0:
                                            if D2O1==1:#Solute that moves across membranes like water 
                                                matrix_C[j][j] += Fcw
                                            else: #Solute that moves across membranes independently of water (the membrane is possibly not one) 
                                                matrix_C[j][j] += Fcw*(1-s_membranes[jmb])
                                        if i not in Apo_w_Zombies0:
                                            if D2O1==1:#Solute that moves across membranes like water 
                                                matrix_C[i][j] -= Fcw
                                            else: #Solute that moves across membranes independently of water (the membrane is possibly not one) 
                                                matrix_C[i][j] -= Fcw*(1-s_membranes[jmb])
                                
                                #Macroscopic distributed parameter for transmembrane flow
                                #Discretization based on cell layers and apoplasmic barriers
                                rank = int(Cell_rank[j-NwallsJun])
                                row = int(rank2row[rank])
                                if rank == 1 and count_epi > 0: #Outer exodermis
                                    row += 1
                                if rank == 3 and count_cortex > 0: #Outer endodermis
                                    if any(passage_cell_ID==array(j-NwallsJun)) and Barrier==2:
                                        row += 2
                                    else:
                                        row += 3
                                elif rank == 3 and count_stele_overall > 0: #Inner endodermis
                                    if any(passage_cell_ID==array(j-NwallsJun)) and Barrier==2:
                                        row += 1
                                        #print('PsiWallPassage:',psi)
                                Flow = K * (psi - psin + s_membranes[jmb]*(Os_walls[i] - Os_cells[j-NwallsJun]))
                                jmb+=1
                                if ((j-NwallsJun not in InterCid) and (j not in listxyl)) or Barrier==0: #No aerenchyma in the elongation zone
                                    if Flow > 0 :
                                        UptakeLayer_plus[row][iMaturity][count] += Flow #grouping membrane flow rates in cell layers
                                    else:
                                        UptakeLayer_minus[row][iMaturity][count] += Flow
                                
                                if K>1.0e-18: #Not an impermeable wall
                                    PsiWallLayer[row][iMaturity][count] += psi
                                    NWallLayer[row][iMaturity][count] += 1
                                
                                if psi_o_cell == inf:
                                    psi_o_cell=psin
                                    ind_o_cell=j
                                else:
                                    if noPD: #No plasmodesmata because the wall i is connected to an intercellular space or xylem vessel
                                        temp=0 #The ordering in PlasmodesmFlowDensity will correspond to the one of ThickWalls except for boderline walls, saved for display only                        
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
                                    PlasmodesmFlowDensity.append(temp/(1.0E-04*height))
                                    PlasmodesmFlowDensity.append(-temp/(1.0E-04*height))
                                    Fcc=temp*1.0E-04*eattr['length']*sign(j-ind_o_cell)
                                    if ind_o_cell<j:
                                        Fcc_list.append((ind_o_cell,j,Fcc)) #(cm^3/d) Water flow rate positive from high index to low index cell
                                    else:
                                        Fcc_list.append((j,ind_o_cell,Fcc))
                                    #if C_flag:
                                    if Sym_Contagion==2: #Convection across plasmodesmata
                                        if Apo_Contagion==2: #Apo & Sym Contagion
                                            if Fcc>0: #Flow from high index to low index cell
                                                if ind_o_cell<j: #From j to ind_o_cell
                                                    if j-NwallsJun not in Sym_Zombie0:
                                                        matrix_C[j][j] -= Fcc
                                                    if ind_o_cell-NwallsJun not in Sym_Zombie0:
                                                        matrix_C[ind_o_cell][j] += Fcc
                                                else: #From ind_o_cell to j
                                                    if ind_o_cell-NwallsJun not in Sym_Zombie0:
                                                        matrix_C[ind_o_cell][ind_o_cell] -= Fcc
                                                    if j-NwallsJun not in Sym_Zombie0:
                                                        matrix_C[j][ind_o_cell] += Fcc
                                            else: #Flow from low index to high index cell
                                                if ind_o_cell<j: #From ind_o_cell to j
                                                    if ind_o_cell-NwallsJun not in Sym_Zombie0:
                                                        matrix_C[ind_o_cell][ind_o_cell] += Fcc
                                                    if j-NwallsJun not in Sym_Zombie0:
                                                        matrix_C[j][ind_o_cell] -= Fcc
                                                else: #From j to ind_o_cell
                                                    if j-NwallsJun not in Sym_Zombie0:
                                                        matrix_C[j][j] += Fcc
                                                    if ind_o_cell-NwallsJun not in Sym_Zombie0:
                                                        matrix_C[ind_o_cell][j] -= Fcc
                                        else: #Only Sym contagion
                                            if Fcc>0: #Flow from high index to low index cell
                                                if ind_o_cell<j: #From j to ind_o_cell
                                                    if j-NwallsJun not in Sym_Zombie0:
                                                        matrix_SymC[j-NwallsJun][j-NwallsJun] -= Fcc
                                                    if ind_o_cell-NwallsJun not in Sym_Zombie0:
                                                        matrix_SymC[ind_o_cell-NwallsJun][j-NwallsJun] += Fcc
                                                else: #From ind_o_cell to j
                                                    if ind_o_cell-NwallsJun not in Sym_Zombie0:
                                                        matrix_SymC[ind_o_cell-NwallsJun][ind_o_cell-NwallsJun] -= Fcc
                                                    if j-NwallsJun not in Sym_Zombie0:
                                                        matrix_SymC[j-NwallsJun][ind_o_cell-NwallsJun] += Fcc
                                            else: #Flow from low index to high index cell
                                                if ind_o_cell<j: #From ind_o_cell to j
                                                    if ind_o_cell-NwallsJun not in Sym_Zombie0:
                                                        matrix_SymC[ind_o_cell-NwallsJun][ind_o_cell-NwallsJun] += Fcc
                                                    if j-NwallsJun not in Sym_Zombie0:
                                                        matrix_SymC[j-NwallsJun][ind_o_cell-NwallsJun] -= Fcc
                                                else: #From j to ind_o_cell
                                                    if j-NwallsJun not in Sym_Zombie0:
                                                        matrix_SymC[j-NwallsJun][j-NwallsJun] += Fcc
                                                    if ind_o_cell-NwallsJun not in Sym_Zombie0:
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
                                K=Kmb[jmb]
                                #Flow densities calculation
                                #Macroscopic distributed parameter for transmembrane flow
                                #Discretization based on cell layers and apoplasmic barriers
                                rank = int(Cell_rank[j-NwallsJun])
                                row = int(rank2row[rank])
                                if rank == 1 and count_epi > 0: #Outer exodermis
                                    row += 1
                                if rank == 3 and count_cortex > 0: #Outer endodermis
                                    if any(passage_cell_ID==array(j-NwallsJun)) and Barrier==2:
                                        row += 2
                                    else:
                                        row += 3
                                elif rank == 3 and count_stele_overall > 0: #Inner endodermis
                                    if any(passage_cell_ID==array(j-NwallsJun)) and Barrier==2:
                                        row += 1
                                Flow = K * (psi - psin + s_membranes[jmb]*(Os_walls[i] - Os_cells[j-NwallsJun]))
                                jmb+=1
                                if ((j-NwallsJun not in InterCid) and (j not in listxyl)) or Barrier==0:
                                    if Flow > 0 :
                                        UptakeLayer_plus[row][iMaturity][count] += Flow #grouping membrane flow rates in cell layers
                                    else:
                                        UptakeLayer_minus[row][iMaturity][count] += Flow
                                
                                if K>1.0e-12: #Not an impermeable wall
                                    PsiWallLayer[row][iMaturity][count] += psi
                                    NWallLayer[row][iMaturity][count] += 1
            
            #if C_flag: #Calculates stationary solute concentration
            if Apo_Contagion==2 or Sym_Contagion==2: #Sym & Apo contagion
                if Apo_Contagion==2 and Sym_Contagion==2: #Sym & Apo contagion
                    #Solving apoplastic & symplastic concentrations
                    soln_C = np.linalg.solve(matrix_C,rhs_C) #Solving the equation to get apoplastic relative concentrations
                elif Apo_Contagion==2:
                    #Solving apoplastic concentrations
                    soln_ApoC = np.linalg.solve(matrix_ApoC,rhs_ApoC) #Solving the equation to get apoplastic & symplastic relative concentrations
                else: # Only Symplastic contagion
                    #Solving apoplastic concentrations
                    soln_SymC = np.linalg.solve(matrix_SymC,rhs_SymC) #Solving the equation to get symplastic relative concentrations
                
                ##Including BC diffusion terms
                #for wid in listxylwalls:
                #    temp=1.0E-04*(lengths[wid]*height)/thickness #Section to length ratio (cm) for the xylem wall
                #    if not temp==0:
                #        list_walls_apo_diff.append(wid)
                #    matrix_C[wid][wid] -= temp*Diff1 #Adding BC diffusion term
                #    rhs_C[wid][0] -= temp*Diff1*Os_xyl[0][count] #new #Adding BC diffusion term
                #for wid in Borderwall:
                #    if (position[wid][0]>=Xcontact) or (Wall2Cell[wid][0]-NwallsJun in Contact): #Wall (not including junctions) connected to soil
                #        temp=1.0E-04*(lengths[wid]/2*height)/(thickness/2)
                #        if not temp==0:
                #            list_walls_apo_diff.append(wid)
                #        matrix_C[wid][wid] -= temp*Diff1 #Adding diffusion BC at soil junction
                #        rhs_C[wid][0] -= temp*Diff1*Os_soil[0][count] #Adding BC diffusion term
                #for jid in Borderjunction:
                #    if (position[jid][0]>=Xcontact) or (Junction2Wall2Cell[jid-Nwalls][0]-NwallsJun in Contact) or (Junction2Wall2Cell[jid-Nwalls][1]-NwallsJun in Contact) or (Junction2Wall2Cell[jid-Nwalls][2]-NwallsJun in Contact): #Junction connected to soil
                #        temp=1.0E-04*(lengths[jid]*height)/(thickness/2)
                #        if not temp==0:
                #            list_walls_apo_diff.append(jid)
                #        matrix_C[jid][jid] -= temp*Diff1 #Adding diffusion BC at soil junction
                #        rhs_C[jid][0] -= temp*Diff1*Os_soil[0][count] #Adding BC diffusion term
                #
                #Nwalls_apo_diff=np.zeros((NwallsJun,2))
                #Nwalls_apo_conv=np.zeros((NwallsJun,2))
                #Nwalls_TM_conv=np.zeros((NwallsJun,2))
                #for wid in list_walls_apo_diff:
                #    Nwalls_apo_diff[wid][1]+=1
                #    Nwalls_apo_diff[wid][0]=wid
                #for wid in list_walls_apo_conv:
                #    Nwalls_apo_conv[wid][1]+=1
                #    Nwalls_apo_conv[wid][0]=wid
                #for wid in list_walls_TM_conv:
                #    Nwalls_TM_conv[wid][1]+=1
                #    Nwalls_TM_conv[wid][0]=wid
                #
                
                #Solving apoplastic concentrations
                #soln_C = np.linalg.solve(matrix_C,rhs_C) #Solving the equation to get potentials inside the network
                
            #Resets matrix_C and rhs_C to geometrical factor values
            if Apo_Contagion==2:
                if Sym_Contagion==2: # Apo & Sym contagion
                    for i,j,Fjw in Fjw_list:
                        if Fjw>0: #Flow from junction to wall
                            if i not in Apo_w_Zombies0:
                                matrix_C[i][j] -= Fjw #Removing convective term
                            if j not in Apo_j_Zombies0:
                                matrix_C[j][j] += Fjw #Removing convective term
                        else: #Flow from wall to junction
                            if i not in Apo_w_Zombies0:
                                matrix_C[i][i] -= Fjw #Removing convective term
                            if j not in Apo_j_Zombies0:
                                matrix_C[j][i] += Fjw #Removing convective term
                else: #Only Apo contagion
                    for i,j,Fjw in Fjw_list:
                        if Fjw>0: #Flow from junction to wall
                            if i not in Apo_w_Zombies0:
                                matrix_ApoC[i][j] -= Fjw #Removing convective term
                            if j not in Apo_j_Zombies0:
                                matrix_ApoC[j][j] += Fjw #Removing convective term
                        else: #Flow from wall to junction
                            if i not in Apo_w_Zombies0:
                                matrix_ApoC[i][i] -= Fjw #Removing convective term
                            if j not in Apo_j_Zombies0:
                                matrix_ApoC[j][i] += Fjw #Removing convective term
            
            if Sym_Contagion==2: #Convection across plasmodesmata
                if Apo_Contagion==2: #Apo & Sym Contagion
                    for i,j,Fcc in Fcc_list:
                        if Fcc>0: #Flow from j to i
                            if j-NwallsJun not in Sym_Zombie0:
                                matrix_C[j][j] += Fcc #Removing convective term
                            if i-NwallsJun not in Sym_Zombie0:
                                matrix_C[i][j] -= Fcc #Removing convective term
                        else: #Flow from i to j
                            if i-NwallsJun not in Sym_Zombie0:
                                matrix_C[i][i] -= Fcc #Removing convective term
                            if j-NwallsJun not in Sym_Zombie0:
                                matrix_C[j][i] += Fcc #Removing convective term
                else: #Only Sym contagion
                    for i,j,Fcc in Fcc_list:
                        if Fcc>0: #Flow from j to i
                            if j-NwallsJun not in Sym_Zombie0:
                                matrix_SymC[j-NwallsJun][j-NwallsJun] += Fcc #Removing convective term
                            if ind_o_cell-NwallsJun not in Sym_Zombie0:
                                matrix_SymC[i-NwallsJun][j-NwallsJun] -= Fcc #Removing convective term
                        else: #Flow from i to j
                            if i-NwallsJun not in Sym_Zombie0:
                                matrix_SymC[i-NwallsJun][i-NwallsJun] -= Fcc #Removing convective term
                            if j-NwallsJun not in Sym_Zombie0:
                                matrix_SymC[j-NwallsJun][i-NwallsJun] += Fcc #Removing convective term
            
            if Apo_Contagion==2 and Sym_Contagion==2:
                for i,j,Fcw,s in Fcw_list:
                    Fcw=-Fcw #Attention, -Fcw was saved
                    if Fcw>0: #Flow from wall to protoplast
                        if i not in Apo_w_Zombies0:
                            if D2O1==1:#Solute that moves across membranes like water 
                                matrix_C[i][i] += Fcw #Removing convective term
                            else: #Solute that moves across membranes independently of water (the membrane is possibly not one) 
                                matrix_C[i][i] += Fcw*(1-s) #Removing convective term
                        if j-NwallsJun not in Sym_Zombie0:
                            if D2O1==1:#Solute that moves across membranes like water 
                                matrix_C[j][i] -= Fcw #Removing convective term
                            else: #Solute that moves across membranes independently of water (the membrane is possibly not one) 
                                matrix_C[j][i] -= Fcw*(1-s) #Removing convective term
                    else: #Flow from protoplast to wall
                        if j-NwallsJun not in Sym_Zombie0:
                            if D2O1==1:#Solute that moves across membranes like water 
                                matrix_C[j][j] -= Fcw #Removing convective term
                            else: #Solute that moves across membranes independently of water (the membrane is possibly not one) 
                                matrix_C[j][j] -= Fcw*(1-s) #Removing convective term
                        if i not in Apo_w_Zombies0:
                            if D2O1==1:#Solute that moves across membranes like water 
                                matrix_C[i][j] += Fcw #Removing convective term
                            else: #Solute that moves across membranes independently of water (the membrane is possibly not one) 
                                matrix_C[i][j] += Fcw*(1-s) #Removing convective term
            
            if Apo_Contagion==2:
                if Sym_Contagion==2: # Apo & Sym contagion
                    i=0
                    for ind in Borderwall:
                        if ind not in Apo_w_Zombies0:
                            Q=Q_soil[i] #(cm^3/d) Positive for water flowing into the root, rhs_s is minus the conductance at the soil root interface
                            if Q<0.0:
                                matrix_C[ind][ind]-=Q #Removing convective term
                        i+=1
                    for ind in Borderjunction:
                        if ind not in Apo_j_Zombies0:
                            Q=Q_soil[i] #(cm^3/d) Positive for water flowing into the root, rhs_s is minus the conductance at the soil root interface
                            if Q<0.0:
                                matrix_C[ind][ind]-=Q #Removing convective term
                        i+=1
                else:
                    i=0
                    for ind in Borderwall:
                        if ind not in Apo_w_Zombies0:
                            Q=Q_soil[i] #(cm^3/d) Positive for water flowing into the root, rhs_s is minus the conductance at the soil root interface
                            if Q<0.0:
                                matrix_ApoC[ind][ind]-=Q #Removing convective term
                        i+=1
                    for ind in Borderjunction:
                        if ind not in Apo_j_Zombies0:
                            Q=Q_soil[i] #(cm^3/d) Positive for water flowing into the root, rhs_s is minus the conductance at the soil root interface
                            if Q<0.0:
                                matrix_ApoC[ind][ind]-=Q #Removing convective term
                        i+=1
            
                ##Removing diffusion terms linked to BC
                #for jid in Borderjunction:
                #    if (position[jid][0]>=Xcontact) or (Junction2Wall2Cell[jid-Nwalls][0]-NwallsJun in Contact) or (Junction2Wall2Cell[jid-Nwalls][1]-NwallsJun in Contact) or (Junction2Wall2Cell[jid-Nwalls][2]-NwallsJun in Contact): #Junction connected to soil
                #        temp=1.0E-04*(lengths[jid]*height)/(thickness/2)
                #        matrix_C[jid][jid] += temp*Diff1 #Removing diffusion BC at soil junction
                #        rhs_C[jid][0] += temp*Diff1*Os_soil[0][count] #Removing BC diffusion term
                #for wid in Borderwall:
                #    if (position[wid][0]>=Xcontact) or (Wall2Cell[wid][0]-NwallsJun in Contact): #Wall (not including junctions) connected to soil
                #        temp=1.0E-04*(lengths[wid]/2*height)/(thickness/2)
                #        matrix_C[wid][wid] += temp*Diff1 #Removing diffusion BC at soil junction
                #        rhs_C[wid][0] += temp*Diff1*Os_soil[0][count] #Removing BC diffusion term
                #for wid in listxylwalls:
                #    temp=1.0E-04*(lengths[wid]*height)/thickness #Section to length ratio (cm) for the xylem wall
                #    matrix_C[wid][wid] += temp*Diff1 #Removing BC diffusion term
                #    rhs_C[wid][0] += temp*Diff1*Os_xyl[0][count] #new #Removing BC diffusion term
                
            
            ####################################
            ## Creates .vtk file for Paraview ##
            ####################################
            
            if Sym_Contagion==1:
                Sym_Zombies=[]
                for source in Sym_source_range:
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
                
                
                text_file = open(newpath+"Sym_Contagion_bottomb"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "w")
                with open(newpath+"Sym_Contagion_bottomb"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "a") as myfile:
                    myfile.write("# vtk DataFile Version 4.0 \n")
                    myfile.write("Contaminated symplastic space geometry \n")
                    myfile.write("ASCII \n")
                    myfile.write(" \n")
                    myfile.write("DATASET UNSTRUCTURED_GRID \n")
                    myfile.write("POINTS "+str(len(ThickWalls))+" float \n")
                    for ThickWallNode in ThickWalls:
                        myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(height/200) + " \n")
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
                
            elif Sym_Contagion==2:
                text_file = open(newpath+"Sym_Contagion_bottomb"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "w")
                with open(newpath+"Sym_Contagion_bottomb"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "a") as myfile:
                    myfile.write("# vtk DataFile Version 4.0 \n")
                    myfile.write("Symplastic hormone concentration \n")
                    myfile.write("ASCII \n")
                    myfile.write(" \n")
                    myfile.write("DATASET UNSTRUCTURED_GRID \n")
                    myfile.write("POINTS "+str(len(ThickWalls))+" float \n")
                    for ThickWallNode in ThickWalls:
                        myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(height/200) + " \n")
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
                    if Apo_Contagion==2:
                        for ThickWallNode in ThickWalls:
                            cellnumber1=ThickWallNode[2]-NwallsJun
                            #print(cellnumber1, soln_C[int(cellnumber1)+NwallsJun])
                            myfile.write(str(float(soln_C[int(cellnumber1+NwallsJun)])) + " \n")
                    else:
                        for ThickWallNode in ThickWalls:
                            cellnumber1=ThickWallNode[2]-NwallsJun
                            #print(cellnumber1, soln_SymC[int(cellnumber1)])
                            myfile.write(str(float(soln_SymC[int(cellnumber1)])) + " \n") #
                myfile.close()
                text_file.close()
                
                #text_file = open(newpath+"Contagion"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "w")
                ##sath01=max(soln[NwallsJun:NwallsJun+Ncells-1])
                ##satl01=min(soln[NwallsJun:NwallsJun+Ncells-1])
                #with open(newpath+"Contagion"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "a") as myfile:
                #    myfile.write("# vtk DataFile Version 4.0 \n")     #("Purchase Amount: %s" % TotalAmount)
                #    myfile.write("Symplastic hormonal spread by convection \n")
                #    myfile.write("ASCII \n")
                #    myfile.write(" \n")
                #    myfile.write("DATASET UNSTRUCTURED_GRID \n")
                #    myfile.write("POINTS "+str(len(G.node))+" float \n")
                #    for node in G:
                #        myfile.write(str(float(position[node][0])) + " " + str(float(position[node][1])) + " " + str(0.0) + " \n")
                #    myfile.write(" \n")
                #    myfile.write("CELLS " + str(Ncells) + " " + str(Ncells*2) + " \n") #
                #    for node, edges in G.adjacency_iter():
                #        i=indice[node]
                #        if i>=NwallsJun: #Cell node
                #            myfile.write("1 " + str(i) + " \n")
                #    myfile.write(" \n")
                #    myfile.write("CELL_TYPES " + str(Ncells) + " \n") #
                #    for node, edges in G.adjacency_iter():
                #        i=indice[node]
                #        if i>=NwallsJun: #Cell node
                #            myfile.write("1 \n")
                #    myfile.write(" \n")
                #    myfile.write("POINT_DATA " + str(len(G.node)) + " \n")
                #    myfile.write("SCALARS Cell_pressure float \n")
                #    myfile.write("LOOKUP_TABLE default \n")
                #    for node in G:
                #        if node-NwallsJun in [Zombie0]: #Source cell
                #            myfile.write(str(float(0.0)) + " \n")
                #        elif node-NwallsJun in Zombies:
                #            myfile.write(str(float(1.0)) + " \n")
                #        else:
                #            myfile.write(str(float(-1.0)) + " \n")
                #myfile.close()
                #text_file.close()
            
            if Apo_Contagion==1:
                Apo_w_Zombies=Apo_w_Zombies0
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
                
                
                text_file = open(newpath+"Apo_Contagion_bottomb"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "w")
                with open(newpath+"Apo_Contagion_bottomb"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "a") as myfile:
                    myfile.write("# vtk DataFile Version 4.0 \n")
                    myfile.write("Contaminated Apoplastic space geometry \n")
                    myfile.write("ASCII \n")
                    myfile.write(" \n")
                    myfile.write("DATASET UNSTRUCTURED_GRID \n")
                    myfile.write("POINTS "+str(len(ThickWallsX))+" float \n")
                    for ThickWallNodeX in ThickWallsX:
                        myfile.write(str(ThickWallNodeX[1]) + " " + str(ThickWallNodeX[2]) + " 0.0 \n")
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
                
            elif Apo_Contagion==2:
                text_file = open(newpath+"Apo_Contagion_bottomb"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "w")
                with open(newpath+"Apo_Contagion_bottomb"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "a") as myfile:
                    myfile.write("# vtk DataFile Version 4.0 \n")
                    myfile.write("Apoplastic hormone concentration \n")
                    myfile.write("ASCII \n")
                    myfile.write(" \n")
                    myfile.write("DATASET UNSTRUCTURED_GRID \n")
                    myfile.write("POINTS "+str(len(ThickWallsX))+" float \n")
                    for ThickWallNodeX in ThickWallsX:
                        myfile.write(str(ThickWallNodeX[1]) + " " + str(ThickWallNodeX[2]) + " 0.0 \n")
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
                    myfile.write("SCALARS Hormone_Symplastic_Relative_Concentration_(-) float \n")
                    myfile.write("LOOKUP_TABLE default \n")
                    if Sym_Contagion==2:
                        Newsoln_C=zeros((len(ThickWallsX),1))
                        j=0
                        for PolygonX in Wall2NewWallX:
                            for id1 in range(int(nWall2NewWallX[j])):
                                Newsoln_C[int(PolygonX[id1])]=soln_C[j]
                            j+=1
                        for i in range(len(ThickWallsX)):
                            myfile.write(str(float(Newsoln_C[i])) + " \n")
                    else:
                        Newsoln_ApoC=zeros((len(ThickWallsX),1))
                        j=0
                        for PolygonX in Wall2NewWallX:
                            for id1 in range(int(nWall2NewWallX[j])):
                                Newsoln_ApoC[int(PolygonX[id1])]=soln_ApoC[j]
                            j+=1
                        for i in range(len(ThickWallsX)):
                            myfile.write(str(float(Newsoln_ApoC[i])) + " \n")
                myfile.close()
                text_file.close()
            
            
            if Paraview==1:
                if ParaviewWP==1: #2D visualization of walls pressure potentials
                    text_file = open(newpath+"Walls2Db"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "w")
                    #sath0=max(soln[0:NwallsJun-1])
                    #satl0=min(soln[0:NwallsJun-1])
                    with open(newpath+"Walls2Db"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "a") as myfile:
                        myfile.write("# vtk DataFile Version 4.0 \n")     #("Purchase Amount: %s" % TotalAmount)
                        myfile.write("Wall geometry 2D \n")
                        myfile.write("ASCII \n")
                        myfile.write(" \n")
                        myfile.write("DATASET UNSTRUCTURED_GRID \n")
                        myfile.write("POINTS "+str(len(G.node))+" float \n")
                        for node in G:
                            myfile.write(str(float(position[node][0])) + " " + str(float(position[node][1])) + " " + str(0.0) + " \n")
                        myfile.write(" \n")
                        myfile.write("CELLS " + str(Nwalls*2-len(list_ghostwalls)*2) + " " + str(Nwalls*6-len(list_ghostwalls)*6) + " \n") #len(G.node)
                        for node, edges in G.adjacency_iter():
                            i=indice[node]
                            if i not in list_ghostwalls:
                                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                                    j=indice[neighboor]
                                    if j>i and eattr['path']=='wall':
                                        #print(nx.get_node_attributes(edges,'path'))
                                        myfile.write(str(2) + " " + str(i) + " " + str(j) + " \n")
                        myfile.write(" \n")
                        myfile.write("CELL_TYPES " + str(Nwalls*2-len(list_ghostwalls)*2) + " \n") #The number of nodes corresponds to the number of wall to wall connections.... to be checked, might not be generality
                        for node, edges in G.adjacency_iter():
                            i=indice[node]
                            if i not in list_ghostwalls:
                                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                                    j=indice[neighboor]
                                    if j>i and eattr['path']=='wall':
                                        #print(nx.get_node_attributes(edges,'path'))
                                        myfile.write(str(3) + " \n") #Line cell type
                        myfile.write(" \n")
                        myfile.write("POINT_DATA " + str(len(G.node)) + " \n")
                        myfile.write("SCALARS Wall_pressure float \n")
                        myfile.write("LOOKUP_TABLE default \n")
                        for node in G:
                            myfile.write(str(float(soln[node])) + " \n") #Line cell type      min(sath0,max(satl0,   ))
                    myfile.close()
                    text_file.close()
                
                if ParaviewWP==1 and ParaviewCP: #2D visualization of walls & cells osmotic potentials
                    text_file = open(newpath+"WallsOsAndCellsOs2Db"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "w")
                    with open(newpath+"WallsOsAndCellsOs2Db"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "a") as myfile:
                        myfile.write("# vtk DataFile Version 4.0 \n")     #("Purchase Amount: %s" % TotalAmount)
                        myfile.write("Wall geometry 2D \n")
                        myfile.write("ASCII \n")
                        myfile.write(" \n")
                        myfile.write("DATASET UNSTRUCTURED_GRID \n")
                        myfile.write("POINTS "+str(len(G.node))+" float \n")
                        for node in G:
                            myfile.write(str(float(position[node][0])) + " " + str(float(position[node][1])) + " " + str(0.0) + " \n")
                        myfile.write(" \n")                                     
                        myfile.write("CELLS " + str(Nwalls*2-len(list_ghostwalls)*2+Ncells) + " " + str(Nwalls*6-len(list_ghostwalls)*6+Ncells*2) + " \n") #
                        for node, edges in G.adjacency_iter():
                            i=indice[node]
                            if i not in list_ghostwalls:
                                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                                    j=indice[neighboor]
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
                                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                                    j=indice[neighboor]
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
                                myfile.write(str(float(Os_walls[i])) + " \n")
                            elif i<NwallsJun: #Junction node
                                myfile.write(str(float(0.0)) + " \n")
                            else: #Cell node
                                myfile.write(str(float(Os_cells[i-NwallsJun])) + " \n")
                    myfile.close()
                    text_file.close()
                    

                
                if ParaviewWP==1 and ParaviewCP==1: #2D visualization of walls & cells water potentials
                    text_file = open(newpath+"WallsAndCells2Db"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "w")
                    with open(newpath+"WallsAndCells2Db"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "a") as myfile:
                        myfile.write("# vtk DataFile Version 4.0 \n")     #("Purchase Amount: %s" % TotalAmount)
                        myfile.write("Water potential distribution in cells and walls 2D \n")
                        myfile.write("ASCII \n")
                        myfile.write(" \n")
                        myfile.write("DATASET UNSTRUCTURED_GRID \n")
                        myfile.write("POINTS "+str(len(G.node))+" float \n")
                        for node in G:
                            myfile.write(str(float(position[node][0])) + " " + str(float(position[node][1])) + " " + str(0.0) + " \n")
                        myfile.write(" \n")
                        myfile.write("CELLS " + str(Nwalls*2-len(list_ghostwalls)*2+Ncells) + " " + str(Nwalls*6-len(list_ghostwalls)*6+Ncells*2) + " \n") #
                        for node, edges in G.adjacency_iter():
                            i=indice[node]
                            if i not in list_ghostwalls:
                                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                                    j=indice[neighboor]
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
                                for neighboor, eattr in edges.items(): #Loop on connections (edges)
                                    j=indice[neighboor]
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
                            myfile.write(str(float(soln[node])) + " \n") #Line cell type
                    myfile.close()
                    text_file.close()
                
                if ParaviewCP==1: #2D visualization of cells water potentials
                    text_file = open(newpath+"Cells2Db"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "w")
                    #sath01=max(soln[NwallsJun:NwallsJun+Ncells-1])
                    #satl01=min(soln[NwallsJun:NwallsJun+Ncells-1])
                    with open(newpath+"Cells2Db"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "a") as myfile:
                        myfile.write("# vtk DataFile Version 4.0 \n")     #("Purchase Amount: %s" % TotalAmount)
                        myfile.write("Pressure potential distribution in cells 2D \n")
                        myfile.write("ASCII \n")
                        myfile.write(" \n")
                        myfile.write("DATASET UNSTRUCTURED_GRID \n")
                        myfile.write("POINTS "+str(len(G.node))+" float \n")
                        for node in G:
                            myfile.write(str(float(position[node][0])) + " " + str(float(position[node][1])) + " " + str(0.0) + " \n")
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
                            myfile.write(str(float(soln[node])) + " \n") #Line cell type      min(sath01,max(satl01,   ))
                    myfile.close()
                    text_file.close()
                    
                
                if ParaviewMF==1: #3D visualization of membrane fluxes
                    text_file = open(newpath+"Membranes3Db"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "w")
                    #sath1=max(MembraneFlowDensity)*color_threshold
                    #satl1=min(MembraneFlowDensity)*color_threshold
                    #if satl1<-sath1: #min(MembraneFlowDensity)<0:
                    #    sath1=-satl1
                    #else:
                    #    satl1=-sath1
                    with open(newpath+"Membranes3Db"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "a") as myfile:
                        myfile.write("# vtk DataFile Version 4.0 \n")
                        myfile.write("Membranes geometry 3D \n")
                        myfile.write("ASCII \n")
                        myfile.write(" \n")
                        myfile.write("DATASET UNSTRUCTURED_GRID \n")
                        myfile.write("POINTS "+str(len(ThickWalls)*2)+" float \n")
                        for ThickWallNode in ThickWalls:
                            myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " 0.0 \n")
                        for ThickWallNode in ThickWalls:
                            myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(height) + " \n")
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
                                myfile.write(str(float(MembraneFlowDensity[int(ThickWallNode[0])])/sperd/cmperm) + " \n") #Flow rate from wall (non junction) to cell   min(sath1,max(satl1,  ))
                            else:
                                myfile.write(str(float((MembraneFlowDensity[int(ThickWallNode[5])]+MembraneFlowDensity[int(ThickWallNode[6])])/2)/sperd/cmperm) + " \n") #Flow rate from junction wall to cell is the average of the 2 neighbouring wall flow rates   min(sath1,max(satl1,  ))
                        for ThickWallNode in ThickWalls:
                            if ThickWallNode[0]<len(MembraneFlowDensity):
                                myfile.write(str(float(MembraneFlowDensity[int(ThickWallNode[0])])/sperd/cmperm) + " \n") #Flow rate from wall (non junction) to cell   min(sath1,max(satl1,  ))
                            else:
                                myfile.write(str(float((MembraneFlowDensity[int(ThickWallNode[5])]+MembraneFlowDensity[int(ThickWallNode[6])])/2)/sperd/cmperm) + " \n") #Flow rate from junction wall to cell is the average of the 2 neighbouring wall flow rates   min(sath1,max(satl1,  ))
                    myfile.close()
                    text_file.close()
                
                if ParaviewWF==1: #Wall flow density data
                    maxWallFlowDensity=0.0
                    for ir in range(int(len(WallFlowDensity))):
                        maxWallFlowDensity=max(maxWallFlowDensity,abs(WallFlowDensity[ir][2]))
                    sath2=maxWallFlowDensity*color_threshold #(1-(1-color_threshold)/2)
                    #satl2=0.0
                    text_file = open(newpath+"WallsThick3D_bottomb"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "w")
                    with open(newpath+"WallsThick3D_bottomb"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "a") as myfile:
                        myfile.write("# vtk DataFile Version 4.0 \n")
                        myfile.write("Wall geometry 3D including thickness bottom \n")
                        myfile.write("ASCII \n")
                        myfile.write(" \n")
                        myfile.write("DATASET UNSTRUCTURED_GRID \n")
                        myfile.write("POINTS "+str(len(ThickWallsX))+" float \n")
                        for ThickWallNodeX in ThickWallsX:
                            myfile.write(str(ThickWallNodeX[1]) + " " + str(ThickWallNodeX[2]) + " 0.0 \n")
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
                                if abs(float(WallFlowDensity[i][2]))>min(NewWallFlowDensity[int(PolygonX[id1])]):
                                    NewWallFlowDensity[int(PolygonX[id1])][0]=max(NewWallFlowDensity[int(PolygonX[id1])])
                                    NewWallFlowDensity[int(PolygonX[id1])][1]=abs(float(WallFlowDensity[i][2]))
                            i+=1
                        for i in range(len(ThickWallsX)):
                            myfile.write(str(float(mean(NewWallFlowDensity[i]))/sperd/cmperm) + " \n")  # min(sath2,  )
                    myfile.close()
                    text_file.close()
                    
                    text_file = open(newpath+"WallsThick3Dcos_bottomb"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "w")
                    with open(newpath+"WallsThick3Dcos_bottomb"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "a") as myfile:
                        myfile.write("# vtk DataFile Version 4.0 \n")
                        myfile.write("Wall geometry 3D including thickness bottom \n")
                        myfile.write("ASCII \n")
                        myfile.write(" \n")
                        myfile.write("DATASET UNSTRUCTURED_GRID \n")
                        myfile.write("POINTS "+str(len(ThickWallsX))+" float \n")
                        for ThickWallNodeX in ThickWallsX:
                            myfile.write(str(ThickWallNodeX[1]) + " " + str(ThickWallNodeX[2]) + " 0.0 \n")
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
                                if abs(float(WallFlowDensity_cos[i][2]))>min(abs(NewWallFlowDensity_cos[int(PolygonX[id1])])):
                                    #Horizontal component of the flux
                                    if abs(NewWallFlowDensity_cos[int(PolygonX[id1])][1])>abs(NewWallFlowDensity_cos[int(PolygonX[id1])][0]): #Keeping the most extreme value
                                        NewWallFlowDensity_cos[int(PolygonX[id1])][0]=NewWallFlowDensity_cos[int(PolygonX[id1])][1]
                                    NewWallFlowDensity_cos[int(PolygonX[id1])][1]=float(WallFlowDensity_cos[i][2])
                            i+=1
                        for i in range(len(ThickWallsX)):
                            myfile.write(str(float(mean(NewWallFlowDensity_cos[i]))/sperd/cmperm) + " \n")  # min(sath2,  )
                    myfile.close()
                    text_file.close()
                
                    if Barrier>0:
                        text_file = open(newpath+"InterC3D_bottomb"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "w")
                        with open(newpath+"InterC3D_bottomb"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "a") as myfile:
                            myfile.write("# vtk DataFile Version 4.0 \n")
                            myfile.write("Intercellular space geometry 3D \n")
                            myfile.write("ASCII \n")
                            myfile.write(" \n")
                            myfile.write("DATASET UNSTRUCTURED_GRID \n")
                            myfile.write("POINTS "+str(len(ThickWalls))+" float \n")
                            for ThickWallNode in ThickWalls:
                                myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(height/200) + " \n")
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
                                    InterCFlowDensity[cid]+=abs(MembraneFlowDensity[int(twpid)])/n #Mean absolute flow density calculation
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
                    text_file = open(newpath+"Plasmodesm3Db"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "w")
                    #sath3=max(PlasmodesmFlowDensity)*color_threshold
                    #satl3=min(PlasmodesmFlowDensity)*color_threshold
                    #if satl3<-sath3: #min(PlasmodesmFlowDensity)<0:
                    #    sath3=-satl3
                    #else:
                    #    satl3=-sath3
                    with open(newpath+"Plasmodesm3Db"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "a") as myfile:
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
                                        myfile.write(str(x) + " " + str(y) + " " + str(z) + " \n")
                            else:
                                break #interrupts the for loop in case we reached the new junction nodes
                        myfile.write(" \n")
                        myfile.write("CELLS " + str(len(PlasmodesmFlowDensity)) + " " + str(len(PlasmodesmFlowDensity)*13) + " \n") #The number of cells corresponds to the number of lines in ThickWalls
                        for i in range(len(PlasmodesmFlowDensity)):
                            if PlasmodesmFlowDensity[i]==0:
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
                                myfile.write(str(float(PlasmodesmFlowDensity[i])/sperd/cmperm) + " \n") #min(sath3,max(satl3, ))
                    myfile.close()
                    text_file.close()
                
                
                if ParaviewMF==1 and ParaviewPF==1: #Membranes and plasmodesms in the same file
                    text_file = open(newpath+"Membranes_n_plasmodesm3Db"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "w")
                    with open(newpath+"Membranes_n_plasmodesm3Db"+str(Barrier)+","+str(iMaturity)+"s"+str(count)+".pvtk", "a") as myfile:
                        myfile.write("# vtk DataFile Version 4.0 \n")
                        myfile.write("Membranes geometry and plasmodesm disks 3D \n")
                        myfile.write("ASCII \n")
                        myfile.write(" \n")
                        myfile.write("DATASET UNSTRUCTURED_GRID \n")
                        myfile.write("POINTS "+str(len(ThickWalls)*2+len(PlasmodesmFlowDensity)*12)+" float \n")
                        for ThickWallNode in ThickWalls:
                            myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " 0.0 \n")
                        for ThickWallNode in ThickWalls:
                            myfile.write(str(ThickWallNode[3]) + " " + str(ThickWallNode[4]) + " " + str(height) + " \n")
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
                                        myfile.write(str(x) + " " + str(y) + " " + str(z) + " \n")
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
                            if PlasmodesmFlowDensity[i]==0:
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
                                myfile.write(str(float(MembraneFlowDensity[int(ThickWallNode[0])])/sperd/cmperm) + " \n") #Flow rate from wall (non junction) to cell
                            else:
                                myfile.write(str(float((MembraneFlowDensity[int(ThickWallNode[5])]+MembraneFlowDensity[int(ThickWallNode[6])])/2)/sperd/cmperm) + " \n") #Flow rate from junction wall to cell is the average of the 2 neighbouring wall flow rates
                        for ThickWallNode in ThickWalls:
                            if ThickWallNode[0]<len(MembraneFlowDensity):
                                myfile.write(str(float(MembraneFlowDensity[int(ThickWallNode[0])])/sperd/cmperm) + " \n") #Flow rate from wall (non junction) to cell
                            else:
                                myfile.write(str(float((MembraneFlowDensity[int(ThickWallNode[5])]+MembraneFlowDensity[int(ThickWallNode[6])])/2)/sperd/cmperm) + " \n") #Flow rate from junction wall to cell is the average of the 2 neighbouring wall flow rates
                        for i in range(len(PlasmodesmFlowDensity)):
                            for j in range(12):
                                myfile.write(str(float(PlasmodesmFlowDensity[i])/sperd/cmperm) + " \n")
                    myfile.close()
                    text_file.close()
        
    
    #write down kr_tot and Uptake distributions in matrices
    iMaturity=-1
    for Maturity in Maturityrange:
        Barrier=int(Maturity.get("Barrier"))
        height=int(Maturity.get("height")) #(microns)
        iMaturity+=1
        text_file = open(newpath+"Macro_prop_"+str(Barrier)+","+str(iMaturity)+".txt", "w")
        with open(newpath+"Macro_prop_"+str(Barrier)+","+str(iMaturity)+".txt", "a") as myfile:
            myfile.write("Macroscopic root radial hydraulic properties, apoplastic barrier "+str(Barrier)+","+str(iMaturity)+" \n")
            myfile.write("\n")
            myfile.write(str(Nscenarios-1)+" scenarios \n")
            myfile.write("\n")
            myfile.write("Cross-section height: "+str(height*1.0E-04)+" cm \n")
            myfile.write("\n")
            myfile.write("Cross-section perimeter: "+str(perimeter[0])+" cm \n")
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
                myfile.write(str(STFlayer_plus[j][iMaturity]*100)+" \n")
            myfile.write("\n")
            myfile.write("Standard Transmembrane release Fractions (%): \n")
            for j in range(int(r_discret[0])):
                myfile.write(str(STFlayer_minus[j][iMaturity]*100)+" \n")
            for i in range(1,Nscenarios):
                myfile.write("\n")
                myfile.write("\n")
                myfile.write("Scenario "+str(i)+" \n")
                myfile.write("\n")
                myfile.write("h_x: "+str(Psi_xyl[iMaturity][i])+" hPa \n")
                myfile.write("\n")
                myfile.write("h_s: "+str(Psi_soil[0][i])+" to "+str(Psi_soil[1][i])+" hPa \n")
                myfile.write("\n")
                myfile.write("h_p: "+str(Psi_sieve[iMaturity][i])+" hPa \n")
                myfile.write("\n")
                myfile.write("O_x: "+str(Os_xyl[0][i])+" to "+str(Os_xyl[1][i])+" hPa \n")
                myfile.write("\n")
                myfile.write("O_s: "+str(Os_soil[0][i])+" to "+str(Os_soil[1][i])+" hPa \n")
                myfile.write("\n")
                myfile.write("O_p: "+str(Os_sieve[0][i])+" hPa \n")
                myfile.write("\n")
                myfile.write("Xcontact: "+str(Xcontact)+" microns \n")
                myfile.write("\n")
                if Barrier==0:
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
                myfile.write("q_tot: "+str(Q_tot[iMaturity][i]/height/1.0E-04)+" cm^2/d \n")
                myfile.write("\n")
                myfile.write("Stele, cortex, and epidermis uptake distribution cm^3/d: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(UptakeLayer_plus[j][iMaturity][i])+" \n")
                myfile.write("\n")
                myfile.write("Stele, cortex, and epidermis release distribution cm^3/d: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(UptakeLayer_minus[j][iMaturity][i])+" \n")
                myfile.write("\n")
                myfile.write("Xylem uptake distribution cm^3/d: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(Q_xyl_layer[j][iMaturity][i])+" \n")
                myfile.write("\n")
                myfile.write("Phloem uptake distribution cm^3/d: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(Q_sieve_layer[j][iMaturity][i])+" \n")
                myfile.write("\n")
                myfile.write("Elongation flow convergence distribution cm^3/d: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(Q_elong_layer[j][iMaturity][i])+" \n")
                myfile.write("\n")
                myfile.write("Cell layers pressure potentials: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(PsiCellLayer[j][iMaturity][i])+" \n")
                myfile.write("\n")
                myfile.write("Cell layers osmotic potentials: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(OsCellLayer[j][iMaturity][i])+" \n")
                myfile.write("\n")
                myfile.write("Wall layers pressure potentials: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(PsiWallLayer[j][iMaturity][i]/NWallLayer[j][iMaturity][i])+" \n")
                myfile.write("\n")
                myfile.write("Wall layers osmotic potentials: \n")
                for j in range(int(r_discret[0])):
                    myfile.write(str(OsWallLayer[j][iMaturity][i])+" \n")
        myfile.close()
        text_file.close()
    
    if Sym_Contagion == 1: #write down results of the hydropatterning study
        iMaturity=-1
        for Maturity in Maturityrange:
            Barrier=int(Maturity.get("Barrier"))
            height=int(Maturity.get("height")) #(microns)
            iMaturity+=1
            text_file = open(newpath+"Hydropatterning_"+str(Barrier)+","+str(iMaturity)+".txt", "w")
            with open(newpath+"Hydropatterning_"+str(Barrier)+","+str(iMaturity)+".txt", "a") as myfile:
                myfile.write("Is there symplastic mass flow from source to target cells? Apoplastic barrier "+str(Barrier)+","+str(iMaturity)+" \n")
                myfile.write("\n")
                myfile.write(str(Nscenarios-1)+" scenarios \n")
                myfile.write("\n")
                myfile.write("Template: "+path+" \n")
                myfile.write("\n")
                myfile.write("Source cell: "+str(Sym_Zombie0)+" \n")
                myfile.write("\n")
                myfile.write("Target cells: "+str(Sym_Target)+" \n")
                myfile.write("\n")
                myfile.write("Immune cells: "+str(Sym_Immune)+" \n")
                myfile.write("\n")
                myfile.write("Cross-section height: "+str(height*1.0E-04)+" cm \n")
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
                if Barrier==0:
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
                    myfile.write("h_x: "+str(Psi_xyl[iMaturity][i])+" hPa, h_s: "+str(Psi_soil[0][i])+" to "+str(Psi_soil[1][i])+" hPa, h_p: "+str(Psi_sieve[iMaturity][i])+" hPa \n")
                    myfile.write("\n")
                    myfile.write("O_x: "+str(Os_xyl[0][i])+" to "+str(Os_xyl[0][i])+" hPa, O_s: "+str(Os_soil[0][i])+" to "+str(Os_soil[1][i])+" hPa, O_p: "+str(Os_sieve[0][i])+" hPa \n")
                    myfile.write("\n")
                    myfile.write("Os_cortex: "+str(Os_cortex[0][count])+" hPa, Os_hetero: "+str(Os_hetero[0][count])+", s_hetero: "+str(s_hetero[0][count])+", s_factor: "+str(s_factor[0][count])+" \n")
                    myfile.write("\n")
                    myfile.write("q_tot: "+str(Q_tot[iMaturity][i]/height/1.0E-04)+" cm^2/d \n")
                    myfile.write("\n")
            myfile.close()
            text_file.close()
    
    if Apo_Contagion == 1: #write down results of the hydrotropism study
        iMaturity=-1
        for Maturity in Maturityrange:
            Barrier=int(Maturity.get("Barrier"))
            height=int(Maturity.get("height")) #(microns)
            iMaturity+=1
            text_file = open(newpath+"Hydrotropism_"+str(Barrier)+","+str(iMaturity)+".txt", "w")
            with open(newpath+"Hydrotropism_"+str(Barrier)+","+str(iMaturity)+".txt", "a") as myfile:
                myfile.write("Is there apoplastic mass flow from source to target cells? Apoplastic barrier "+str(Barrier)+","+str(iMaturity)+" \n")
                myfile.write("\n")
                myfile.write(str(Nscenarios-1)+" scenarios \n")
                myfile.write("\n")
                myfile.write("Template: "+path+" \n")
                myfile.write("\n")
                myfile.write("Source cell: "+str(Apo_Zombie0)+" \n")
                myfile.write("\n")
                myfile.write("Target cells: "+str(Apo_Target)+" \n")
                myfile.write("\n")
                myfile.write("Immune cells: "+str(Apo_Immune)+" \n")
                myfile.write("\n")
                myfile.write("Cross-section height: "+str(height*1.0E-04)+" cm \n")
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
                if Barrier==0:
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
                    myfile.write("h_x: "+str(Psi_xyl[iMaturity][i])+" hPa, h_s: "+str(Psi_soil[0][i])+" to "+str(Psi_soil[1][i])+" hPa, h_p: "+str(Psi_sieve[iMaturity][i])+" hPa \n")
                    myfile.write("\n")
                    myfile.write("O_x: "+str(Os_xyl[0][i])+" to "+str(Os_xyl[0][i])+" hPa, O_s: "+str(Os_soil[0][i])+" to "+str(Os_soil[1][i])+" hPa, O_p: "+str(Os_sieve[0][i])+" hPa \n")
                    myfile.write("\n")
                    myfile.write("Os_cortex: "+str(Os_cortex[0][count])+" hPa, Os_hetero: "+str(Os_hetero[0][count])+", s_hetero: "+str(s_hetero[0][count])+", s_factor: "+str(s_factor[0][count])+" \n")
                    myfile.write("\n")
                    myfile.write("q_tot: "+str(Q_tot[iMaturity][i]/height/1.0E-04)+" cm^2/d \n")
                    myfile.write("\n")
            myfile.close()
            text_file.close()
        
    
    #text_file = open(newpath+"Cortical_cell_perimeters.txt", "w")
    #with open(newpath+"Cortical_cell_perimeters.txt", "a") as myfile:
    #    for j in range(10):
    #        myfile.write(str(cortex_cellperimeters[j][:])+" \n")
    #myfile.close()
    #text_file.close()
    
    
    #text_file = open(newpath+"PD_flow_rates.txt", "w")
    #with open(newpath+"Cortical_cell_perimeters.txt", "a") as myfile:
    #    for j in range(10):
    #        myfile.write(str(cortex_cellperimeters[j][:])+" \n")
    #myfile.close()
    #text_file.close()
