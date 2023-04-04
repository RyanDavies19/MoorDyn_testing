import ctypes
import os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
import matplotlib as mpl
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import moorpy as mp
from moorpy.helpers import set_axes_equal, read_mooring_file
import moordyn


class load_infile():
    # Built from system class and contained functions in MoorPy dev branch

    def __init__(self, in_dirname = "", out_dirname = "", rootname = "", extension = ""):

        
        self.in_dirname = in_dirname
        self.out_dirname = out_dirname
        self.rootname = rootname
        self.extension = extension
        self.in_file = str(in_dirname+rootname+extension)
        self.qs = 0

        # lists to hold mooring system objects
        self.bodyList = []
        self.rodList = []  # note: Rods are currently only fully supported when plotting MoorDyn output, not in MoorPy modeling
        # <<< TODO: add support for Rods eventually, for compatability with MoorDyn systems
        self.pointList = []
        self.lineList = []
        self.lineTypes = {}
        self.rodTypes = {}

        # For RMS comparison
        self.version1 = ""
        self.version2 = ""
        
        self.MDoptions = {} # dictionary that can hold any MoorDyn options read in from an input file, so they can be saved in a new MD file if need be
        self.outputList = [] # Output channels for unload

        # read in data from an input file if a filename was provided
        if len(self.in_file) > 0:
            self.load()

            # For plotting later
            self.depth = float(self.MDoptions["WtrDpth"])

            # So that timeseries data is not loaded multiple times by plot_start_end and plot_animation
            self.md_branch_old = ""
            self.first_plot = True # TODO: This and timeseries_not_loaded are redundant 
            self.timeseries_not_loaded = True

        else:
            print("Warning: input file empty, water depth unknown for plotting ouputs")

    def load(self):
        '''Loads a MoorPy System from a MoorDyn-style input file

        Parameters
        ----------
        filename : string
            the file name of a MoorDyn-style input file.

        Raises
        ------
        ValueError
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        
        # create/empty the lists to start with

        RodDict   = {}  # create empty dictionary for rod types
        self.lineTypes = {}  # create empty dictionary for line types
        self.rodTypes = {}  # create empty dictionary for line types

        # ensure the mooring system's object lists are empty before adding to them
        self.bodyList = []
        self.rodList  = []
        self.pointList= []
        self.lineList = []

        # Values to calculate w value for line and rods, artifact of MoorPy, not necessary for use here
        g = 9.81
        rho = 1025.0

        # Artifact from MoorPy, for debugging printing
        self.display = 0    # a flag that controls how much printing occurs in methods within the System (Set manually. Values > 0 cause increasing output.)

        
        # the ground body (number 0, type 1[fixed]) never moves but is the parent of all anchored things
        self.groundBody = mp.Body(self, 0, 1, np.zeros(6))   # <<< implementation not complete <<<< be careful here if/when MoorPy is split up

        # number of fairleads, for determining vector size in running MoorDyn
        self.numfair = 0 
        
        # figure out if it's a YAML file or MoorDyn-style file based on the extension, then open and process
        print('attempting to read '+self.in_file)

        # assuming normal form
        f = open(self.in_file, 'r')

        # read in the data
        
        for line in f:          # loop through each line in the file

            # get line type property sets
            if line.count('---') > 0 and (line.upper().count('LINE DICTIONARY') > 0 or line.upper().count('LINE TYPES') > 0):
                line = next(f) # skip this header line, plus channel names and units lines
                line = next(f)
                line = next(f)
                while line.count('---') == 0:
                    entries = line.split()  # entries: TypeName   Diam    Mass/m     EA     BA/-zeta    EI         Cd     Ca     CdAx    CaAx
                    #self.addLineType(entries[0], float(entries[1]), float(entries[2]), float(entries[3])) 
                    
                    type_string = entries[0]
                    d    = float(entries[1])
                    mass = float(entries[2])
                    w = (mass - np.pi/4*d**2 *rho)*g                        
                    lineType = dict(name=type_string, d_vol=d, w=w, m=mass)  # make dictionary for this rod type
                    
                    # support linear (EA) or nonlinear (filename string) option for elasticity
                    #if there is a text file in the EA input 
                    if entries[3].find(".txt") != -1:
                        #Then we read in ten-strain file
                        ten_str_fname = entries[3]
                        ten_str = open(ten_str_fname[1:-1], 'r') 
                        
                        #Read line in ten-strain file until we hit '---' signifying the end of the file
                        for line in ten_str:
                                #skip first 3 lines (Header for input file)
                                line = next(ten_str)
                                line = next(ten_str)
                                line = next(ten_str)
                                #Preallocate Arrays
                                str_array = []
                                ten_array = []
                                #Loop through lines until you hit '---' signifying the end of the file 
                                while line.count('---') == 0:
                                    ten_str_entries = line.split() #split entries ten_str_entries: strain tension
                                    str_array.append(ten_str_entries[0]) #First one is strain
                                    ten_array.append(ten_str_entries[1]) #Second one is tension
                                    line = next(ten_str) #go to next line
                        lineType['Str'] = str_array #make new entry in the dictionary to carry tension and strain arrays
                        lineType['Ten'] = ten_array

                    else:

                        try:
                            lineType['EA'] = float(entries[3].split('|')[0])         # get EA, and only take first value if multiples are given
                        except:
                            lineType['EA'] = 1e9
                            print('EA entry not recognized - using placeholder value of 1000 MN')
                    if len(entries) >= 10: # read in other elasticity and hydro coefficients as well if enough columns are provided
                        lineType['BA'  ] = float(entries[4].split('|')[0])
                        lineType['EI'  ] = float(entries[5])
                        lineType['Cd'  ] = float(entries[6])
                        lineType['Ca'  ] = float(entries[7])
                        lineType['CdAx'] = float(entries[8])
                        lineType['CaAx'] = float(entries[9])
                        lineType['material'] = type_string
                    else:
                        print("WARNING: line elasticity and hydro coeffecients not read")
                    
                    if type_string in self.lineTypes:                         # if there is already a line type with this name
                        self.lineTypes[type_string].update(lineType)          # update the existing dictionary values rather than overwriting with a new dictionary
                    else:
                        self.lineTypes[type_string] = lineType
                    
                    line = next(f)
                    
                    
            # get rod type property sets
            if line.count('---') > 0 and (line.upper().count('ROD DICTIONARY') > 0 or line.upper().count('ROD TYPES') > 0):
                line = next(f) # skip this header line, plus channel names and units lines
                line = next(f)
                line = next(f)
                while line.count('---') == 0:
                    entries = line.split()  # entries: TypeName      Diam     Mass/m    Cd     Ca      CdEnd    CaEnd
                    #RodTypesName.append(entries[0]) # name string
                    #RodTypesD.append(   entries[1]) # diameter
                    #RodDict[entries[0]] = entries[1] # add dictionary entry with name and diameter
                    
                    type_string = entries[0]
                    d    = float(entries[1])
                    mass = float(entries[2])
                    w = (mass - np.pi/4*d**2 *rho)*g
                    
                    rodType = dict(name=type_string, d_vol=d, w=w, m=mass)  # make dictionary for this rod type
                    
                    if len(entries) >= 7: # read in hydro coefficients as well if enough columns are provided
                        rodType['Cd'   ] = float(entries[3])
                        rodType['Ca'   ] = float(entries[4])
                        rodType['CdEnd'] = float(entries[5])
                        rodType['CaEnd'] = float(entries[6])
                    else:
                        print("WARNING: hydro coefficients not read in")
                    
                    if type_string in self.rodTypes:                        # if there is already a rod type with this name
                        self.rodTypes[type_string].update(rodType)          # update the existing dictionary values rather than overwriting with a new dictionary
                    else:
                        self.rodTypes[type_string] = rodType
                    
                    line = next(f)
                    
                    
            # get properties of each Body
            if line.count('---') > 0 and (line.upper().count('BODIES') > 0 or line.upper().count('BODY LIST') > 0 or line.upper().count('BODY PROPERTIES') > 0):
                line = next(f) # skip this header line, plus channel names and units lines
                line = next(f)
                line = next(f)
                while line.count('---') == 0:
                    entries = line.split()  # entries: ID   Attachment  X0  Y0  Z0  r0  p0  y0    M  CG*  I*    V  CdA*  Ca*            
                    num = int(entries[0])
                    entry0 = entries[1].lower()                         
                    #num = np.int_("".join(c for c in entry0 if not c.isalpha()))  # remove alpha characters to identify Body #
                    
                    if ("fair" in entry0) or ("coupled" in entry0) or ("ves" in entry0):       # coupled case
                        bodyType = -1                        
                    elif ("con" in entry0) or ("free" in entry0):                              # free case
                        bodyType = 0
                    else:                                                                      # for now assuming unlabeled free case
                        bodyType = 0
                        # if we detected there were unrecognized chars here, could: raise ValueError(f"Body type not recognized for Body {num}")
                    #bodyType = -1   # manually setting the body type as -1 for FAST.Farm SM investigation
                    
                    r6  = np.array(entries[2:8], dtype=float)   # initial position and orientation [m, rad]
                    r6[3:] = r6[3:]*np.pi/180.0                 # convert from deg to rad
                    #rCG = np.array(entries[7:10], dtype=float)  # location of body CG in body reference frame [m]
                    m = np.float_(entries[8])                   # mass, centered at CG [kg]
                    v = np.float_(entries[11])                   # volume, assumed centered at reference point [m^3]
                    
                    # process CG
                    strings_rCG = entries[ 9].split("|")                   # split by braces, if any
                    if len(strings_rCG) == 1:                              # if only one entry, it is the z coordinate
                        rCG = np.array([0.0, 0.0, float(strings_rCG[0])])
                    elif len(strings_rCG) == 3:                            # all three coordinates provided
                        rCG = np.array(strings_rCG, dtype=float)
                    else:
                        raise Exception(f"Body {num} CG entry (col 10) must have 1 or 3 numbers.")
                        
                    # process mements of inertia
                    strings_I = entries[10].split("|")                     # split by braces, if any
                    if len(strings_I) == 1:                                # if only one entry, use it for all directions
                        Inert = np.array(3*strings_I, dtype=float)
                    elif len(strings_I) == 3:                              # all three coordinates provided
                        Inert = np.array(strings_I, dtype=float)
                    else:
                        raise Exception(f"Body {num} inertia entry (col 11) must have 1 or 3 numbers.")
                    
                    # process drag ceofficient by area product
                    strings_CdA = entries[12].split("|")                   # split by braces, if any
                    if len(strings_CdA) == 1:                              # if only one entry, use it for all directions
                        CdA = np.array(3*strings_CdA, dtype=float)
                    elif len(strings_CdA) == 3:                            # all three coordinates provided
                        CdA = np.array(strings_CdA, dtype=float)
                    else:
                        raise Exception(f"Body {num} CdA entry (col 13) must have 1 or 3 numbers.")
                    
                    # process added mass coefficient
                    strings_Ca = entries[13].split("|")                    # split by braces, if any				
                    if len(strings_Ca) == 1:                               # if only one entry, use it for all directions
                        Ca = np.array(strings_Ca, dtype=float)
                    elif len(strings_Ca) == 3:                             #all three coordinates provided
                        Ca = np.array(strings_Ca, dtype=float)
                    else:
                        raise Exception(f"Body {num} Ca entry (col 14) must have 1 or 3 numbers.")
                    
                    # add the body
                    self.bodyList.append(mp.Body(self, num, bodyType, r6, m=m, v=v, rCG=rCG, I=Inert, CdA=CdA, Ca=Ca) )
                                
                    line = next(f)
                    
                    
            # get properties of each rod
            if line.count('---') > 0 and (line.upper().count('RODS') > 0 or line.upper().count('ROD LIST') > 0 or line.upper().count('ROD PROPERTIES') > 0):
                line = next(f) # skip this header line, plus channel names and units lines
                line = next(f)
                line = next(f)
                while line.count('---') == 0:
                    entries = line.split()  # entries: RodID  RodType  Attachment  Xa   Ya   Za   Xb   Yb   Zb  NumSegs  Flags/Outputs
                    num = int(entries[0])
                    rodType = self.rodTypes[entries[1]]
                    attachment = entries[2].lower()
                    dia = rodType['d_vol']  # find diameter based on specified rod type string
                    rA = np.array(entries[3:6], dtype=float)
                    rB = np.array(entries[6:9], dtype=float)
                    nSegs = int(entries[9])
                    # >>> note: this is currently only set up for use with MoorDyn output data <<<
                    
                    if nSegs==0:       # this is the zero-length special case
                        lUnstr = 0
                        self.rodList.append(mp.Point(self, num, 0, rA) )
                    else:
                        lUnstr = np.linalg.norm(rB-rA)
                        self.rodList.append(mp.Line(self, num, lUnstr, rodType, nSegs=nSegs, isRod=1) )
                        
                        if ("body" in attachment) or ("turbine" in attachment):
                            # attach to body here
                            BodyID = int("".join(filter(str.isdigit, attachment)))
                            if len(self.bodyList) < BodyID:
                                self.bodyList.append(mp.Body(self, 1, 0, np.zeros(6)))
                                
                            self.bodyList[BodyID-1].attachRod(num, np.hstack([rA,rB]))
                            
                        else: # (in progress - unsure if htis works) <<<
                            self.rodList[-1].rA = rA #.setEndPosition(rA, 0)  # set initial end A position
                            self.rodList[-1].rB = rB #.setEndPosition(rB, 1)  # set initial end B position
                        
                    line = next(f)
                    
            
            # get properties of each Point
            if line.count('---') > 0 and (line.upper().count('POINTS') > 0 or line.upper().count('POINT LIST') > 0 or line.upper().count('POINT PROPERTIES') > 0 or line.upper().count('CONNECTION PROPERTIES') > 0 or line.upper().count('NODE PROPERTIES') > 0):
                line = next(f) # skip this header line, plus channel names and units lines
                line = next(f)
                line = next(f)
                while line.count('---') == 0:
                    entries = line.split()         # entries:  ID   Attachment  X       Y     Z      Mass   Volume  CdA    Ca
                    entry0 = entries[0].lower()          
                    entry1 = entries[1].lower() 
                    
                    
                    num = np.int_("".join(c for c in entry0 if not c.isalpha()))  # remove alpha characters to identify Point #
                    
                    
                    if ("anch" in entry1) or ("fix" in entry1):
                        pointType = 1
                        # attach to ground body for ease of identifying anchors
                        self.groundBody.attachPoint(num,entries[2:5]) 
                        
                    elif ("body" in entry1) or ("turbine" in entry1):
                        pointType = 1
                        # attach to body here
                        BodyID = int("".join(filter(str.isdigit, entry1)))
                        if len(self.bodyList) < BodyID:
                            self.bodyList.append(mp.Body(self, 1, 0, np.zeros(6)))
                        
                        rRel = np.array(entries[2:5], dtype=float)
                        self.bodyList[BodyID-1].attachPoint(num, rRel)
                        
                    elif ("fair" in entry1) or ("ves" in entry1) or ("couple" in entry1):
                        # for coupled point type, just set it up that same way in MoorPy (attachment to a body not needed, right?)
                        pointType = -1     
                        self.numfair += 1                       
                        '''
                        # attach to a generic platform body (and make it if it doesn't exist)
                        if len(self.bodyList) > 1:
                            raise ValueError("Generic Fairlead/Vessel-type points aren't supported when multiple bodies are defined.")
                        if len(self.bodyList) == 0:
                            #print("Adding a body to attach fairlead points to.")
                            self.bodyList.append( Body(self, 1, 0, np.zeros(6)))#, m=m, v=v, rCG=rCG) )
                        
                        rRel = np.array(entries[2:5], dtype=float)
                        self.bodyList[0].attachPoint(num, rRel)    
                        '''
                       
                            
                    elif ("con" in entry1) or ("free" in entry1):
                        pointType = 0
                    else:
                        print("Point type not recognized")
                    
                    if 'seabed' in entries[4]:
                        entries[4] = -self.depth
                    r = np.array(entries[2:5], dtype=float)
                    m  = float(entries[5])
                    v  = float(entries[6])
                    CdA= float(entries[7])
                    Ca = float(entries[8])
                    self.pointList.append(mp.Point(self, num, pointType, r, m=m, v=v, CdA=CdA, Ca=Ca) )
                    line = next(f)
                    
                    
            # get properties of each line
            if line.count('---') > 0 and (line.upper().count('LINES') > 0 or line.upper().count('LINE LIST') > 0 or line.upper().count('LINE PROPERTIES') > 0):
                line = next(f) # skip this header line, plus channel names and units lines
                line = next(f)
                line = next(f)
                while line.count('---') == 0:
                    entries = line.split()  # entries: ID  LineType  AttachA  AttachB  UnstrLen  NumSegs   Outputs
                                            
                    num    = np.int_(entries[0])
                    lUnstr = np.float_(entries[4])
                    lineType = self.lineTypes[entries[1]]
                    nSegs  = np.int_(entries[5])         
                    
                    #lineList.append( Line(dirName, num, lUnstr, dia, nSegs) )
                    self.lineList.append(mp.Line(self, num, lUnstr, lineType, nSegs=nSegs)) #attachments = [int(entries[4]), int(entries[5])]) )
                    
                    # attach end A
                    numA = int("".join(filter(str.isdigit, entries[2])))  # get number from the attachA string
                    if entries[2][0] in ['r','R']:    # if id starts with an "R" or "Rod"  
                        if numA <= len(self.rodList) and numA > 0:
                            if entries[2][-1] in ['a','A']:
                                self.rodList[numA-1].attachLine(num, 0)  # add line (end A, denoted by 0) to rod >>end A, denoted by 0<<
                            elif entries[2][-1] in ['b','B']: 
                                self.rodList[numA-1].attachLine(num, 0)  # add line (end A, denoted by 0) to rod >>end B, denoted by 1<<
                            else:
                                raise ValueError(f"Rod end (A or B) must be specified for line {num} end A attachment. Input was: {entries[2]}")
                        else:
                            raise ValueError(f"Rod ID ({numA}) out of bounds for line {num} end A attachment.") 
                    
                    else:     # if J starts with a "C" or "Con" or goes straight ot the number then it's attached to a Connection
                        if numA <= len(self.pointList) and numA > 0:  
                            self.pointList[numA-1].attachLine(num, 0)  # add line (end A, denoted by 0) to Point
                        else:
                            raise ValueError(f"Point ID ({numA}) out of bounds for line {num} end A attachment.") 

                    # attach end B
                    numB = int("".join(filter(str.isdigit, entries[3])))  # get number from the attachA string
                    if entries[3][0] in ['r','R']:    # if id starts with an "R" or "Rod"  
                        if numB <= len(self.rodList) and numB > 0:
                            if entries[3][-1] in ['a','A']:
                                self.rodList[numB-1].attachLine(num, 1)  # add line (end B, denoted by 1) to rod >>end A, denoted by 0<<
                            elif entries[3][-1] in ['b','B']: 
                                self.rodList[numB-1].attachLine(num, 1)  # add line (end B, denoted by 1) to rod >>end B, denoted by 1<<
                            else:
                                raise ValueError(f"Rod end (A or B) must be specified for line {num} end B attachment. Input was: {entries[2]}")
                        else:
                            raise ValueError(f"Rod ID ({numB}) out of bounds for line {num} end B attachment.") 
                    
                    else:     # if J starts with a "C" or "Con" or goes straight ot the number then it's attached to a Connection
                        if numB <= len(self.pointList) and numB > 0:  
                            self.pointList[numB-1].attachLine(num, 1)  # add line (end B, denoted by 1) to Point
                        else:
                            raise ValueError(f"Point ID ({numB}) out of bounds for line {num} end B attachment.") 

                    line = next(f)  # advance to the next line

            # get options entries
            if line.count('---') > 0 and "options" in line.lower():
                line = next(f) # skip this header line
                
                while line.count('---') == 0:
                    entries = line.split()       
                    entry0 = entries[0] #.lower() 
                    entry1 = entries[1] #.lower() 
                    
                    # also store a dict of all parameters that can be regurgitated during an unload
                    self.MDoptions[entry1] = entry0
                    
                    line = next(f)
            if  line.count('---') > 0 and ("outputs" or "output") in line.lower():
                line = next(f) # skip this header line
            
                while line.count('---') == 0: 
                    self.outputList.append(line)
                    line = next(f)

        f.close()  # close data file
        
        print(f"Mooring input file '{self.in_file}' loaded successfully.")

    def unload(self, fileName, MDversion_out=2, MDversion_in = 2, line_dL=0, rod_dL=0, flag='p', outputList = []):

        # Does not work for rods 

        '''Unloads a MoorPy system into a MoorDyn-style input file

        Parameters
        ----------
        fileName : string
            file name of output file to hold MoorPy System.
        line_dL : float, optional
            Optional specified for target segment length when discretizing Lines
        rod_dL : float, optional
            Optional specified for target segment length when discretizing Rods
        outputList : list of strings, optional
            Optional list of additional requested output channels

        Returns
        -------
        None.

        '''
        if MDversion_out==1:
            #For version MoorDyn v1

            #Collection of default values, each can be customized when the method is called
            
            # Set up the dictionary that will be used to write the OPTIONS section
            # MDoptionsDict = dict(dtM=0.001, kb=3.0e6, cb=3.0e5, TmaxIC=60)        # start by setting some key default values
            # # Other available options: Echo=False, dtIC=2, CdScaleIC=10, threshIC=0.01
            # MDoptionsDict = self.MDoptions                                 # update the dict with any settings saved from an input file

            # Some default settings to fill in if coefficients aren't set
            #lineTypeDefaults = dict(BA=-1.0, EI=0.0, Cd=1.2, Ca=1.0, CdAx=0.2, CaAx=0.0)
            # lineTypeDefaults = dict(BA=-1.0, cIntDamp=-0.8, EI=0.0, Can=1.0, Cat=1.0, Cdn=1.0, Cdt=0.5)
            # rodTypeDefaults  = dict(Cd=1.2, Ca=1.0, CdEnd=1.0, CaEnd=1.0)
            lineTypeDefaults = {}
            rodTypeDefaults = {}
            
            # bodyDefaults = dict(IX=0, IY=0, IZ=0, CdA_xyz=[0,0,0], Ca_xyz=[0,0,0])
            
            # Figure out mooring line attachments (Create a ix2 array of connection points from a list of m points)
            connection_points = np.empty([len(self.lineList),2])                   #First column is Anchor Node, second is Fairlead node
            for point_ind,point in enumerate(self.pointList,start = 1):                    #Loop through all the points
                for (line,line_pos) in zip(point.attached,point.attachedEndB):          #Loop through all the lines #s connected to this point
                    if line_pos == 0:                                                       #If the A side of this line is connected to the point
                        connection_points[line -1,0] = point_ind                                #Save as as an Anchor Node
                        #connection_points[line -1,0] = self.pointList.index(point) + 1
                    elif line_pos == 1:                                                     #If the B side of this line is connected to the point
                        connection_points[line -1,1] = point_ind                                #Save as a Fairlead node
                        #connection_points[line -1,1] = self.pointList.index(point) + 1
            
            #Outputs List
            Outputs = []
            # Outputs = [f"FairTen{i+1}" for i in range(len(self.lineList))]        # for now, have a fairlead tension output for each line
            #Outputs.append("Con2Fz","Con3Fz","Con6Fz","Con7Fz","Con10Fz","Con11Fz","L3N20T","L6N20T","L9N20T")
   
            
   
            print('attempting to write '+fileName +' for MoorDyn v'+str(MDversion_out))
            #Array to add strings to for each line of moordyn input file
            L = []                   
            
            
            # Generate text for the MoorDyn input file
            L.append('Output file generated by run_all_scripts')
                
            
            #L.append("---------------------- LINE TYPES -----------------------------------------------------")
            L.append("---------------------- LINE DICTIONARY -----------------------------------------------------")
            #L.append(f"{len(self.lineTypes)}    NTypes   - number of LineTypes")
            #L.append("LineType         Diam     MassDen   EA        cIntDamp     EI     Can    Cat    Cdn    Cdt")
            #L.append("   (-)           (m)      (kg/m)    (N)        (Pa-s)    (N-m^2)  (-)    (-)    (-)    (-)")
            L.append("LineType         Diam     MassDenInAir   EA        BA/-zeta     Can    Cat    Cdn    Cdt")
            L.append("   (-)           (m)        (kg/m)       (N)       (Pa-s/-)     (-)    (-)    (-)    (-)")

            if MDversion_in ==1:
                for key, lineType in self.lineTypes.items(): 
                    di = lineTypeDefaults.copy()  # start with a new dictionary of just the defaults
                    di.update(lineType)           # then copy in the lineType's existing values
                    print(di)
                    #L.append("{:<12} {:7.4f} {:8.2f}  {:7.3e} {:7.3e} {:7.3e}   {:<7.3f} {:<7.3f} {:<7.2f} {:<7.2f}".format(
                            #key, di['d_vol'], di['m'], di['EA'], di['cIntDamp'], di['EI'], di['Can'], di['Cat'], di['Cdn'], di['Cdt']))
                    L.append("{:<12} {:7.4f} {:8.2f}  {:7.3e} {:7.3e}       {:<7.3f} {:<7.3f} {:<7.2f} {:<7.2f}".format(
                            key, di['d_vol'], di['m'], di['EA'], di['BA'], di['Can'], di['Cat'], di['Cdn'], di['Cdt']))
            elif MDversion_in == 2:

                # Untested, likely not working

                for key, lineType in self.lineTypes.items(): 
                    di = lineTypeDefaults.copy()  # start with a new dictionary of just the defaults
                    di.update(lineType)           # then copy in the lineType's existing values
                    L.append("{:<12} {:7.4f} {:8.2f}  {:7.3e} {:7.3e}  {:<7.3f} {:<7.3f} {:<7.2f} {:<7.2f}".format(
                            key, di['d_vol'], di['m'], di['EA'], di['BA'], di['Ca'], di['CaAx'], di['Cd'], di['CdAx']))
                
            
            
            #L.append("---------------------- POINTS ---------------------------------------------------------")
            L.append("---------------------- NODE PROPERTIES ---------------------------------------------------------")
            #L.append(f"{len(self.pointList)}    NConnects   - number of connections including anchors and fairleads")
            L.append("Node    Type         X        Y        Z        M      V      FX     FY     FZ    CdA    CA ")
            L.append("(-)     (-)         (m)      (m)      (m)      (kg)   (m^3)  (kN)   (kN)   (kN)   (m^2)  (-)")
            #L.append("ID  Attachment     X       Y       Z          Mass   Volume  CdA    Ca")
            #L.append("(#)   (-)         (m)     (m)     (m)         (kg)   (m^3)  (m^2)   (-)")
            
            for point in self.pointList:
                point_pos = point.r             # get point position in global reference frame to start with
                if point.type == 1:             # point is fixed or attached (anch, body, fix)
                    point_type = 'Fixed'
                    
                    #Check if the point is attached to body
                    for body in self.bodyList:
                        for attached_Point in body.attachedP:
                            
                            if attached_Point == point.number:
                                #point_type = "Body" + str(body.number)
                                point_type = "Vessel"
                                point_pos = body.rPointRel[body.attachedP.index(attached_Point)]   # get point position in the body reference frame
                    
                elif point.type == 0:           # point is coupled externally (con, free)
                    point_type = 'Connect'
                        
                elif point.type == -1:          # point is free to move (fair, ves)
                    point_type = 'Vessel'
                
                L.append("{:<4d} {:9} {:8.2f} {:8.2f} {:8.2f} {:9.2f} {:6.2f} {:6.2f} {:6.2f} {:6.2f} {:6.2f} {:6.2f}".format(
                          point.number,point_type, point_pos[0],point_pos[1],point_pos[2], point.m, point.v, 0, 0, 0, point.CdA, point.Ca))
                
            
            #L.append("---------------------- LINES -----------------------------------------------------")
            L.append("---------------------- LINE PROPERTIES -----------------------------------------------------")
            #L.append(f"{len(self.lineList)}    NLines   - number of line objects")
            #L.append("Line      LineType   UnstrLen  NumSegs  AttachA  AttachB  Outputs")
            #L.append("(-)         (-)       (m)        (-)     (-)      (-)     (-)")
            #L.append("ID    LineType      AttachA  AttachB  UnstrLen  NumSegs  LineOutputs")
            #L.append("(#)    (name)        (#)      (#)       (m)       (-)     (-)")
            L.append("Line      LineType   UnstrLen  NumSegs  NodeAnch  NodeFair  Flags/Outputs")
            L.append("(-)         (-)       (m)        (-)      (-)       (-)         (-)")

            for i,line in enumerate(self.lineList):
                L.append("{:<4d} {:<15} {:8.3f} {:5d} {:7d} {:8d}      {}"
                        .format(line.number, line.type['name'], line.L, line.nNodes-1, int(connection_points[i,0]), int(connection_points[i,1]), flag))

            
            #L.append("---------------------- OPTIONS ----------------------------------------")
            L.append("---------------------- SOLVER OPTIONS ----------------------------------------")

            for key, val in self.MDoptions.items():
                L.append(f"{val:<15}  {key}")
            
            """
            #Solver Options
            L.append("{:<9.3f}dtM          - time step to use in mooring integration (s)".format(float(dtm)))
            L.append("{:<9.0e}kbot           - bottom stiffness (Pa/m)".format(kbot))
            L.append("{:<9.0e}cbot           - bottom damping (Pa-s/m)".format(cbot))
            L.append("{:<9.0f}dtIC      - time interval for analyzing convergence during IC gen (s)".format(int(dtIC)))
            L.append("{:<9.0f}TmaxIC      - max time for ic gen (s)".format(int(TmaxIC)))
            L.append("{:<9.0f}CdScaleIC      - factor by which to scale drag coefficients during dynamic relaxation (-)".format(int(CdScaleIC)))
            L.append("{:<9.2f}threshIC      - threshold for IC convergence (-)".format(threshIC))
            
            #Failure Header
            """
            
            L.append("--------------------------- OUTPUTS --------------------------------------------")
            
            Outputs = Outputs+outputList   # add any user-specified outputs passed to unload
            
            for Output in Outputs:
                L.append(Output)
            #L.append("END")
                
                
            L.append('--------------------- need this line ------------------')
            
            
            #Write the text file
            with open(fileName, 'w') as out:
                for x in range(len(L)):
                    out.write(L[x])
                    out.write('\n')
            
            print('Successfully written '+fileName +' input file using MoorDyn v1')
       
    def v1_build(self, outfile  = "Mooring/lines.txt", flag='p'):
        # read in data from an input file if a filename was provided
        
        if os.path.exists("Mooring/lines.txt"): #protect from overwriting
            os.system("mv Mooring/lines.txt Mooring/lines_archived.txt")
            print("Mooring/lines.txt moved to Mooring/lines_archived.txt to protect from overwriting.")
            print("Mooring/lines_archived.txt will be overwritten next run of run_all_scripts.py")

        self.unload(outfile, MDversion_out = 1, MDversion_in = 2, flag = flag, outputList = self.outputList)

    def plot(self, ax=None, bounds='default', rbound=0, color=None, **kwargs):
        '''Plots the mooring system objects in their current positions
        Parameters
        ----------
        bounds : string, optional
            signifier for the type of bounds desired in the plot. The default is "default".
        ax : axes, optional
            Plot on an existing set of axes
        color : string, optional
            Some way to control the color of the plot ... TBD <<<
        hidebox : bool, optional
            If true, hides the axes and box so just the plotted lines are visible.
        rbound : float, optional
            A bound to be placed on each axis of the plot. If 0, the bounds will be the max values on each axis. The default is 0.
        title : string, optional
            A title of the plot. The default is "".
        linelabels : bool, optional
            Adds line numbers to plot in text. Default is False.
        pointlabels: bool, optional
            Adds point numbers to plot in text. Default is False.
        endpoints: bool, optional
            Adds visible end points to lines. Default is False.
        bathymetry: bool, optional
            Creates a bathymetry map of the seabed based on an input file. Default is False.
            
        Returns
        -------
        fig : figure object
            To hold the axes of the plot
        ax: axis object
            To hold the points and drawing of the plot
            
        '''
        
        # kwargs that can be used for plot or plot2d
        title           = kwargs.get('title'          , ""        )     # optional title for the plot
        time            = kwargs.get("time"           , 0         )     # the time in seconds of when you want to plot
        linelabels      = kwargs.get('linelabels'     , False     )     # toggle to include line number labels in the plot
        pointlabels     = kwargs.get('pointlabels'    , False     )     # toggle to include point number labels in the plot
        draw_body       = kwargs.get("draw_body"      , True      )     # toggle to draw the Bodies or not
        draw_anchors    = kwargs.get('draw_anchors'   , False     )     # toggle to draw the anchors of the mooring system or not  
        bathymetry      = kwargs.get("bathymetry"     , False     )     # toggle (and string) to include bathymetry or not. Can do full map based on text file, or simple squares
        cmap_bath       = kwargs.get("cmap"           , 'ocean'   )     # matplotlib colormap specification
        alpha           = kwargs.get("opacity"        , 1.0       )     # the transparency of the bathymetry plot_surface
        rang            = kwargs.get('rang'           , 'hold'    )     # colorbar range: if range not used, set it as a placeholder, it will get adjusted later
        cbar_bath       = kwargs.get('cbar_bath'      , False     )     # toggle to include a colorbar for a plot or not
        colortension    = kwargs.get("colortension"   , False     )     # toggle to draw the mooring lines in colors based on node tensions
        cmap_tension    = kwargs.get('cmap_tension'   , 'rainbow' )     # the type of color spectrum desired for colortensions
        cbar_tension    = kwargs.get('cbar_tension'   , False     )     # toggle to include a colorbar of the tensions when colortension=True
        figsize         = kwargs.get('figsize'        , (6,4)     )     # the dimensions of the figure to be plotted
        # kwargs that are currently only used in plot
        hidebox         = kwargs.get('hidebox'        , False     )     # toggles whether to show the axes or not
        endpoints       = kwargs.get('endpoints'      , False     )     # toggle to include the line end points in the plot
        waterplane      = kwargs.get("waterplane"     , False     )     # option to plot water surface
        shadow          = kwargs.get("shadow"         , True      )     # toggle to draw the mooring line shadows or not
        cbar_bath_size  = kwargs.get('colorbar_size'  , 1.0       )     # the scale of the colorbar. Not the same as aspect. Aspect adjusts proportions
        # bound kwargs
        xbounds         = kwargs.get('xbounds'        , None      )     # the bounds of the x-axis. The midpoint of these bounds determines the origin point of orientation of the plot
        ybounds         = kwargs.get('ybounds'        , None      )     # the bounds of the y-axis. The midpoint of these bounds determines the origin point of orientation of the plot
        zbounds         = kwargs.get('zbounds'        , None      )     # the bounds of the z-axis. The midpoint of these bounds determines the origin point of orientation of the plot
        
        # sort out bounds
        xs = []
        ys = []
        zs = [0, -self.depth]
        
        for point in self.pointList:
            xs.append(point.r[0])
            ys.append(point.r[1])
            zs.append(point.r[2])
            

        # if axes not passed in, make a new figure
        if ax == None:    
            fig = plt.figure(figsize=figsize)
            ax = plt.axes(projection='3d')
        else:
            fig = ax.get_figure()

        if color == None:
            if time == 0:
                color = "r"
            else:
                color = "b"
        
        # set bounds
        if rbound==0:
            if len(xs) > 0:
                rbound = max([max(xs), max(ys), -min(xs), -min(ys)]) # this is the most extreme coordinate
            else:
                rbound = self.depth
                
        # set the DATA bounds on the axis
        if bounds=='default':
            ax.set_zlim([-self.depth, 0])
        elif bounds=='rbound':   
            ax.set_xlim([-rbound,rbound])
            ax.set_ylim([-rbound,rbound])
            ax.set_zlim([-rbound, rbound])
        elif bounds=='mooring':
            ax.set_xlim([-rbound,0])
            ax.set_ylim([-rbound/2,rbound/2])
            ax.set_zlim([-self.depth, 0])
        
        # set the AXIS bounds on the axis (changing these bounds can change the perspective of the matplotlib figure)
        if (np.array([xbounds, ybounds, zbounds]) != None).any():
            ax.autoscale(enable=False,axis='both')
        if xbounds != None:
            ax.set_xbound(xbounds[0], xbounds[1])
        if ybounds != None:
            ax.set_ybound(ybounds[0], ybounds[1])
        if zbounds != None:
            ax.set_zbound(zbounds[0], zbounds[1])
        
        # draw things
        if draw_body:
            for body in self.bodyList:
                body.draw(ax)
        
        for rod in self.rodList:
            if len(self.rodList)==0:    # usually, there are no rods in the rodList
                pass
            else:
                #if self.qs==0 and len(rod.Tdata) == 0:
                #    pass
                if isinstance(rod, mp.Line) and rod.show:
                    rod.drawLine(time, ax, color=color, shadow=shadow)
                #if isinstance(rod, Point):  # zero-length special case
                #    not plotting points for now
            
        
        if draw_anchors:
            for line in self.lineList:
                if line.zp[0,0]==-self.depth:
                    itime = int(time/line.dt)
                    r = [line.xp[itime,0], line.yp[itime,0], line.zp[itime,0]]
                    if color==None:
                        c='tab:blue'
                    else:
                        c=color
                    plt.plot(r[0], r[1], r[2], 'v', color=c, markersize=5)
        
        j = 0
        for line in self.lineList:
            if self.qs==0 and len(line.Tdata) == 0:
                pass
            else:
                j = j + 1
                if color==None and 'material' in line.type:
                    if 'chain' in line.type['material'] or 'Cadena80' in line.type['material']:
                        line.drawLine(time, ax, color=[.1, 0, 0], endpoints=endpoints, shadow=shadow, colortension=colortension, cmap_tension=cmap_tension)
                    elif 'rope' in line.type['material'] or 'polyester' in line.type['material'] or 'Dpoli169' in line.type['material']:
                        line.drawLine(time, ax, color=[.3,.5,.5], endpoints=endpoints, shadow=shadow, colortension=colortension, cmap_tension=cmap_tension)
                    elif 'nylon' in line.type['material']:
                        line.drawLine(time, ax, color=[.8,.8,.2], endpoints=endpoints, shadow=shadow, colortension=colortension, cmap_tension=cmap_tension)
                    else:
                        line.drawLine(time, ax, color=[0.5,0.5,0.5], endpoints=endpoints, shadow=shadow, colortension=colortension, cmap_tension=cmap_tension)
                else:
                    line.drawLine(time, ax, color=color, endpoints=endpoints, shadow=shadow, colortension=colortension, cmap_tension=cmap_tension)
                
                # Add line labels 
                if linelabels == True:
                    ax.text((line.rA[0]+line.rB[0])/2, (line.rA[1]+line.rB[1])/2, (line.rA[2]+line.rB[2])/2, j)
            
        if cbar_tension:
            maxten = max([max(line.getLineTens()) for line in self.lineList])   # find the max tension in the System
            minten = min([min(line.getLineTens()) for line in self.lineList])   # find the min tension in the System
            bounds = range(int(minten),int(maxten), int((maxten-minten)/256)) 
            norm = mpl.colors.BoundaryNorm(bounds, 256)     # set the bounds in a norm object, with 256 being the length of all colorbar strings
            fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap_tension), label='Tension (N)')  # add the colorbar
            fig.tight_layout()
            
        # Add point labels
        i = 0
        for point in self.pointList:
            points = []
            i = i + 1
            if pointlabels == True:
                ax.text(point.r[0], point.r[1], point.r[2], i, c = 'r')
            
            if bathymetry==True:     # if bathymetry is true, then make squares at each anchor point
                if point.attachedEndB[0] == 0 and point.r[2] < -400:
                    points.append([point.r[0]+250, point.r[1]+250, point.r[2]])
                    points.append([point.r[0]+250, point.r[1]-250, point.r[2]])
                    points.append([point.r[0]-250, point.r[1]-250, point.r[2]])
                    points.append([point.r[0]-250, point.r[1]+250, point.r[2]])
                    
                    Z = np.array(points)
                    verts = [[Z[0],Z[1],Z[2],Z[3]]]
                    ax.add_collection3d(Poly3DCollection(verts, facecolors='limegreen', linewidths=1, edgecolors='g', alpha=alpha))
            
        if isinstance(bathymetry, str):   # or, if it's a string, load in the bathymetry file

            # parse through the MoorDyn bathymetry file
            bathGrid_Xs, bathGrid_Ys, bathGrid = self.readBathymetryFile(bathymetry)
            if rang=='hold':
                rang = (np.min(-bathGrid), np.max(-bathGrid))
            '''
            # First method: plot nice 2D squares using Poly3DCollection
            nX = len(bathGrid_Xs)
            nY = len(bathGrid_Ys)
            # store a list of points in the grid
            Z = [[bathGrid_Xs[j],bathGrid_Ys[i],-bathGrid[i,j]] for i in range(nY) for j in range(nX)]
            # plot every square in the grid (e.g. 16 point grid yields 9 squares)
            verts = []
            for i in range(nY-1):
                for j in range(nX-1):
                    verts.append([Z[j+nX*i],Z[(j+1)+nX*i],Z[(j+1)+nX*(i+1)],Z[j+nX*(i+1)]])
                    ax.add_collection3d(Poly3DCollection(verts, facecolors='limegreen', linewidths=1, edgecolors='g', alpha=0.5))
                    verts = []
            '''
            # Second method: plot a 3D surface, plot_surface
            X, Y = np.meshgrid(bathGrid_Xs, bathGrid_Ys)
            
            bath = ax.plot_surface(X,Y,-bathGrid, cmap=cmap_bath, vmin=rang[0], vmax=rang[1], alpha=alpha)
            
            if cbar_bath_size!=1.0:    # make sure the colorbar is turned on just in case it isn't when the other colorbar inputs are used
                cbar_bath=True
            if cbar_bath:
                fig.colorbar(bath, shrink=cbar_bath_size, label='depth (m)')
        
        # draw water surface if requested
        if waterplane:
            waterXs = np.array([min(xs), max(xs)])
            waterYs = np.array([min(ys), max(ys)])
            waterX, waterY = np.meshgrid(waterXs, waterYs)
            ax.plot_surface(waterX, waterY, np.array([[-50,-50],[-50,-50]]), alpha=0.5)
    
        
        fig.suptitle(title)
        
        set_axes_equal(ax)
        
        ax.set_zticks([-self.depth, 0])  # set z ticks to just 0 and seabed
        
        if hidebox:
            ax.axis('off')
        
        return fig, ax  # return the figure and axis object in case it will be used later to update the plot

    def plot2d(self, Xuvec=[1,0,0], Yuvec=[0,0,1], ax=None, color=None, **kwargs):
        '''Makes a 2D plot of the mooring system objects in their current positions

        Parameters
        ----------
        Xuvec : list, optional
            plane at which the x-axis is desired. The default is [1,0,0].
        Yuvec : lsit, optional
            plane at which the y-axis is desired. The default is [0,0,1].            
        ax : axes, optional
            Plot on an existing set of axes
        color : string, optional
            Some way to control the color of the plot ... TBD <<<
        title : string, optional
            A title of the plot. The default is "".

        Returns
        -------
        fig : figure object
            To hold the axes of the plot
        ax: axis object
            To hold the points and drawing of the plot

        '''
        
        # kwargs that can be used for plot or plot2d
        title            = kwargs.get('title'           , ""        )   # optional title for the plot
        time             = kwargs.get("time"            , 0         )   # the time in seconds of when you want to plot
        linelabels       = kwargs.get('linelabels'      , False     )   # toggle to include line number labels in the plot
        pointlabels      = kwargs.get('pointlabels'     , False     )   # toggle to include point number labels in the plot
        draw_body        = kwargs.get("draw_body"       , False     )   # toggle to draw the Bodies or not
        draw_anchors     = kwargs.get('draw_anchors'    , False     )   # toggle to draw the anchors of the mooring system or not   
        bathymetry       = kwargs.get("bathymetry"      , False     )   # toggle (and string) to include bathymetry contours or not based on text file
        cmap_bath        = kwargs.get("cmap_bath"       , 'ocean'   )   # matplotlib colormap specification
        alpha            = kwargs.get("opacity"         , 1.0       )   # the transparency of the bathymetry plot_surface
        rang             = kwargs.get('rang'            , 'hold'    )   # colorbar range: if range not used, set it as a placeholder, it will get adjusted later
        cbar_bath        = kwargs.get('colorbar'        , False     )   # toggle to include a colorbar for a plot or not
        colortension     = kwargs.get("colortension"    , False     )   # toggle to draw the mooring lines in colors based on node tensions
        cmap_tension     = kwargs.get('cmap_tension'    , 'rainbow' )   # the type of color spectrum desired for colortensions
        cbar_tension     = kwargs.get('cbar_tension'    , False     )   # toggle to include a colorbar of the tensions when colortension=True
        figsize          = kwargs.get('figsize'         , (6,4)     )   # the dimensions of the figure to be plotted
        # kwargs that are currently only used in plot2d
        levels           = kwargs.get("levels"          , 7         )   # the number (or array) of levels in the contour plot
        cbar_bath_aspect = kwargs.get('cbar_bath_aspect', 20        )   # the proportion of the colorbar. Default is 20 height x 1 width
        cbar_bath_ticks  = kwargs.get('cbar_bath_ticks' , None      )   # the desired tick labels on the colorbar (can be an array)
        plotnodes        = kwargs.get('plotnodes'       , []        )   # the list of node numbers that are desired to be plotted
        plotnodesline    = kwargs.get('plotnodesline'   , []        )   # the list of line numbers that match up with the desired node to be plotted
        label            = kwargs.get('label'           , ""        )   # the label/marker name of a line in the System
        draw_fairlead    = kwargs.get('draw_fairlead'   , False     )   # toggle to draw large points for the fairleads
        
        
        
        # if axes not passed in, make a new figure
        if ax == None:
            fig, ax = plt.subplots(1,1, figsize=figsize)
        else:
            fig = ax.get_figure()

        if color == None:
            if time == 0:
                color = "r"
            else:
                color = "b"

        
        if draw_body:
            for body in self.bodyList:
                #body.draw(ax)
                r = body.r6[0:3]
                x = r[Xuvec.index(1)]
                y = r[Yuvec.index(1)]
                plt.plot(x, y, 'ko', markersize=5)
            
        for rod in self.rodList:
            if isinstance(rod, mp.Line):
                rod.drawLine2d(time, ax, color=color, Xuvec=Xuvec, Yuvec=Yuvec)
            
        if draw_fairlead:
            for line in self.lineList:
                if line.number==1:
                    itime = int(time/line.dt)
                    r = [line.xp[itime,-1], line.yp[itime,-1], line.zp[itime,-1]]
                    x = r[Xuvec.index(1)]
                    y = r[Yuvec.index(1)]
                    if color==None:
                        c='tab:blue'
                    else:
                        c=color
                    plt.plot(x, y, 'o', color=c, markersize=5)
                
                
        if draw_anchors:
            for line in self.lineList:
                if line.zp[0,0]==-self.depth:
                    itime = int(time/line.dt)
                    r = [line.xp[itime,0], line.yp[itime,0], line.zp[itime,0]]
                    x = r[Xuvec.index(1)]
                    y = r[Yuvec.index(1)]
                    if color==None:
                        c='tab:blue'
                    else:
                        c=color
                    plt.plot(x, y, 'v', color=c, markersize=5)
            
        
        j = 0
        for line in self.lineList:
            if line!=self.lineList[0]:
                label=""
            j = j + 1
            if color==None and 'material' in line.type:
                if 'chain' in line.type['material']:
                    line.drawLine2d(time, ax, color=[.1, 0, 0], Xuvec=Xuvec, Yuvec=Yuvec, colortension=colortension, cmap=cmap_tension, plotnodes=plotnodes, plotnodesline=plotnodesline, label=label, alpha=alpha)
                elif 'rope' in line.type['material'] or 'polyester' in line.type['material']:
                    line.drawLine2d(time, ax, color=[.3,.5,.5], Xuvec=Xuvec, Yuvec=Yuvec, colortension=colortension, cmap=cmap_tension, plotnodes=plotnodes, plotnodesline=plotnodesline, label=label, alpha=alpha)
                else:
                    line.drawLine2d(time, ax, color=[0.3,0.3,0.3], Xuvec=Xuvec, Yuvec=Yuvec, colortension=colortension, cmap=cmap_tension, plotnodes=plotnodes, plotnodesline=plotnodesline, label=label, alpha=alpha)
            else:
                line.drawLine2d(time, ax, color=color, Xuvec=Xuvec, Yuvec=Yuvec, colortension=colortension, cmap=cmap_tension, plotnodes=plotnodes, plotnodesline=plotnodesline, label=label, alpha=alpha)

            # Add Line labels
            if linelabels == True:
                xloc = np.dot([(line.rA[0]+line.rB[0])/2, (line.rA[1]+line.rB[1])/2, (line.rA[2]+line.rB[2])/2],Xuvec)
                yloc = np.dot([(line.rA[0]+line.rB[0])/2, (line.rA[1]+line.rB[1])/2, (line.rA[2]+line.rB[2])/2],Yuvec)
                ax.text(xloc,yloc,j)
        
        if cbar_tension:
            maxten = max([max(line.getLineTens()) for line in self.lineList])   # find the max tension in the System
            minten = min([min(line.getLineTens()) for line in self.lineList])   # find the min tension in the System
            bounds = range(int(minten),int(maxten), int((maxten-minten)/256)) 
            norm = mpl.colors.BoundaryNorm(bounds, 256)     # set the bounds in a norm object, with 256 being the length of all colorbar strings
            fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap_tension), label='Tension (N)')  # add the colorbar
            fig.tight_layout()
        
        # Add point labels
        i = 0 
        for point in self.pointList:
            i = i + 1
            if pointlabels == True:
                xloc = np.dot([point.r[0], point.r[1], point.r[2]], Xuvec)
                yloc = np.dot([point.r[0], point.r[1], point.r[2]], Yuvec)
                ax.text(xloc, yloc, i, c = 'r')
        
        if isinstance(bathymetry, str):   # or, if it's a string, load in the bathymetry file

            # parse through the MoorDyn bathymetry file
            bathGrid_Xs, bathGrid_Ys, bathGrid = self.readBathymetryFile(bathymetry)
            
            X, Y = np.meshgrid(bathGrid_Xs, bathGrid_Ys)
            Z = -bathGrid
            if rang=='hold':
                rang = (np.min(Z), np.max(Z))
            
            Xind = Xuvec.index(1); Yind = Yuvec.index(1); Zind = int(3-Xind-Yind)
            W = [X,Y,Z]
            
            # plot a contour profile of the bathymetry
            bath = ax.contourf(W[Xind],W[Yind],W[Zind], cmap=cmap_bath, levels=levels, alpha=alpha, vmin=rang[0], vmax=rang[1])
            
            if cbar_bath_aspect!=20 or cbar_bath_ticks!=None:    # make sure the colorbar is turned on just in case it isn't when the other colorbar inputs are used
                cbar_bath=True
            if cbar_bath:
                fig.colorbar(bath, label='depth (m)', aspect=cbar_bath_aspect, ticks=cbar_bath_ticks)
            
        
        ax.axis("equal")
        ax.set_title(title)
        
        return fig, ax  # return the figure and axis object in case it will be used later to update the plot
        
    def updateCoords(self, tStep, colortension, cmap_tension, label, dt): 
        '''Update animation function. This gets called by animateLines every iteration of the animation and 
        redraws the lines and rods in their next positions.'''
        
        for rod in self.rodList:

            if isinstance(rod, mp.Line) and rod.show:  # draw it if MoorPy is representing it as as Rod-Line object, and it's set to be shown
                rod.redrawLine(-tStep)
            
        for line in self.lineList:
            if len(line.Tdata) > 0:
                line.redrawLine(-tStep, colortension=colortension, cmap_tension=cmap_tension)
        
        label.set_text(f'time={np.round(tStep*dt,1)}')
            
        return 

    def animateLines(self, figure = None, axis = None, color=None, interval=200, repeat=True, delay=0, runtime=-1, start_end_ani = False, **kwargs):
        '''
        Parameters
        ----------
        dirname : string
            The name of the directory folder you are in.
        rootname : string
            The name of the front portion of the main file name, like spar_WT1, or DTU_10MW_NAUTILUS_GoM.
        interval : int, optional
            The time between animation frames in milliseconds. The default is 200.
        repeat : bool, optional
            Whether or not to repeat the animation. The default is True.
        delay : int, optional
            The time between consecutive animation runs in milliseconds. The default is 0.
        runtime: int, optional
            The desired time that the animation should run to in seconds. The default is -1, which means to run the full simulation

        Returns
        -------
        line_ani : animation
            an animation of the mooring lines based off of MoorDyn data.
            Needs to be stored, returned, and referenced in a variable
        '''
        
        bathymetry       = kwargs.get('bathymetry'     , False     )     # toggles whether to show the axes or not
        opacity          = kwargs.get('opacity'        , 1.0       )     # the transparency of the bathymetry plot_surface
        hidebox          = kwargs.get('hidebox'        , False     )     # toggles whether to show the axes or not
        rang             = kwargs.get('rang'           , 'hold'    )     # colorbar range: if range not used, set it as a placeholder, it will get adjusted later
        speed            = kwargs.get('speed'          , 10        )     # the resolution of the animation; how fluid/speedy the animation is
        colortension     = kwargs.get("colortension"   , False     )     # toggle to draw the mooring lines in colors based on node tensions
        cmap_tension     = kwargs.get('cmap_tension'   , 'rainbow' )     # the type of color spectrum desired for colortensions
        draw_body        = kwargs.get('draw_body'      , True      )
        # bound kwargs
        xbounds          = kwargs.get('xbounds'        , None      )     # the bounds of the x-axis. The midpoint of these bounds determines the origin point of orientation of the plot
        ybounds          = kwargs.get('ybounds'        , None      )     # the bounds of the y-axis. The midpoint of these bounds determines the origin point of orientation of the plot
        zbounds          = kwargs.get('zbounds'        , None      )     # the bounds of the z-axis. The midpoint of these bounds determines the origin point of orientation of the plot

                # can use any other kwargs that go into self.plot()
        
        if self.qs==1:
            raise ValueError("This System is set to be quasi-static. Import MoorDyn data and make qs=0 to use this method")   

        # create the figure and axes to draw the animation
        fig, ax = self.plot(ax=axis, color=color, draw_body=draw_body, bathymetry=bathymetry, opacity=opacity, hidebox=hidebox, rang=rang, 
                            colortension=colortension, xbounds=xbounds, ybounds=ybounds, zbounds=zbounds)
    
        ax.set_xlabel('x');    ax.set_ylabel('y');      ax.set_zlabel('z');
        label = ax.text(-100, 100, 0, 'time=0', ha='center', va='center', fontsize=10, color="k")
   
        # find idyn, the index of the first line in the lineList that contains a series of Tdata
        idyn = len(self.lineList)-1    # note, the idyn approach is not robust to different Lines having output, or Rods. Should reconsider.
        for line in self.lineList:
            if len(line.Tdata) > 0:
                idyn = line.number-1
                break
        
        if runtime==-1:     # if the full simulation is desired to be animated
            nFrames = len(self.lineList[idyn].Tdata)
        else:               # if only part of the simulation is to be animated, up to a certain runtime in seconds
            itime = int(np.where(self.lineList[idyn].Tdata==runtime)[0])
            nFrames = len(self.lineList[idyn].Tdata[0:itime])
        
        dt = self.lineList[idyn].Tdata[1]-self.lineList[idyn].Tdata[0] # time step length (s) <<< should get this from main MoorDyn output file <<<
        
        # Animation: update the figure with the updated coordinates from update_Coords function
        # NOTE: the animation needs to be stored in a variable, return out of the method, and referenced when calling self.animatelines()
        if start_end_ani:
            line_ani = animation.FuncAnimation(fig, self.updateCoords, np.array([1, nFrames-1]), fargs=(colortension, cmap_tension, label, dt),
                                                interval=1500, repeat=repeat, repeat_delay=delay, blit=False)
                                                # works well when np.arange(...nFrames...) is used. Others iterable ways to do this
        else:
            line_ani = animation.FuncAnimation(fig, self.updateCoords, np.arange(1, nFrames-1, speed), fargs=(colortension, cmap_tension, label, dt),
                                        interval = interval, repeat=repeat, repeat_delay=delay, blit=False)
                                        # works well when np.arange(...nFrames...) is used. Others iterable ways to do this

        return line_ani
 
    def plot_animation(self, fig = None, ax = None, path = "", rootname = "", extension = "", color = "", start_end_ani = False):
        
        if len(path+rootname+extension)==0 or len(rootname)==0:
            raise ValueError("The MoorDyn input file name and the root name of the MoorDyn output files (e.g. the .fst file name without extension) need to be given.")

        for line in self.lineList:
            line.loadData(path, rootname, sep='_')
        
        for rod in self.rodList:
            if isinstance(rod, mp.Line):
                rod.loadData(path, rootname, sep='_')  
            
        animation = self.animateLines(figure=fig, axis=ax, color=color, interval=500, repeat=True, delay=0, runtime=-1, start_end_ani=start_end_ani)
        return animation
    
    def plot_start_end(self, plot2d = False, fig = None, ax= None, path = "", rootname = "", extension = "", md_branch = "", color = None, tMax = 10):
    
        if len(path+rootname+extension)==0 or len(rootname)==0:
            raise ValueError("The MoorDyn input file name and the root name of the MoorDyn output files (e.g. the .fst file name without extension) need to be given.")

        if md_branch != self.md_branch_old or self.first_plot:
            for line in self.lineList:
                line.loadData(path, rootname, sep='_')  
            
            for rod in self.rodList:
                if isinstance(rod, mp.Line):
                    rod.loadData(path, rootname, sep='_')  

            self.first_plot = False
            self.md_branch_old = md_branch
        
        tMax = self.lineList[0].Tdata[-1]
        tMin = self.lineList[0].Tdata[0]

        if type(ax) == type(np.array([])):
            ax1 = ax[0]
            ax2 = ax[1]

            if plot2d:
                self.plot2d(ax=ax1, bounds='default', rbound=0, color=color, time = tMin)
                self.plot2d(ax=ax2, bounds='default', rbound=0, color=color, time = tMax-1)
            else:
                self.plot(ax=ax1, bounds='default', rbound=0, color=color, time = tMin)
                self.plot(ax=ax2, bounds='default', rbound=0, color=color, time = tMax-1)

            ax[0].set_title("tMin = {}".format(tMin))
            ax[1].set_title("tMax-1 = {}".format(tMax-1))

            fig.suptitle("Key: v2old is red, v2new is yellow, dev2 is blue, and v1 is green", fontsize = 12)

        else:
            if plot2d:
                self.plot2d(ax=ax, bounds='default', rbound=0, color = 'r', time = tMin)
                self.plot2d(ax=ax, bounds='default', rbound=0, color= 'b', time = tMax-1)
            else:
                self.plot(ax=ax, bounds='default', rbound=0, color= 'r', time = tMin)
                self.plot(ax=ax, bounds='default', rbound=0, color= 'b', time = tMax-1)
            
            ax.set_title(md_branch)
            fig.suptitle("first and last timesteps. Key: tMin = {} is red, tMax-1 = {} is blue".format(tMin, tMax), fontsize = 12)
    
    def rms_comparison (self, control = "dev2", display = True, save = False, run_v1 = True, run_v2n = True, run_v2o = True, line_rmse = True, ten_rmse = True, plot_ten = True, path = ''):
        
        self.control = self.rootname+control
        self.out_dirname = path

        # TODO: make dirname, rootname, versioname, line, etc all class variables to clear up the warnings
        # TODO: make read_mooring_file outputs all class variables so files are only loaded once for plotting. Best way would probably be to have as seperate function that is called if they havent been loaded (Called from figures or called from rmse_comparison)
        # TODO: add self to both the new functions

        version_list = [self.rootname+'v2new', self.rootname+'v2old', self.rootname+'dev2', 'Line']
        if not run_v1:
            version_list.remove('Line')
        if not run_v2n:
            version_list.remove(self.rootname+'v2new')
        if not run_v2o:
            version_list.remove(self.rootname+'v2old')

        if line_rmse:
            fig1, ax1 = plt.subplots(3,1, sharey=True, sharex= True)
        if ten_rmse:
            fig2, ax2 = plt.subplots(3,1, sharex= True)

        if plot_ten:
            fig3, ax3 = plt.subplots(1,1)

        i = 0
        for version in version_list:
            if control not in version:
                if line_rmse:
                    self.test = version
                    self.compare_all_lines(ax1[i])
                    self.timeseries_not_loaded = True
                    
                if ten_rmse:
                    if version == 'Line':
                        self.test = 'Line_Lines'
                    else:
                        self.test = version
                    self.read_tensions(ax2[i])
                if ten_rmse or line_rmse:
                    i += 1 
            if plot_ten:
                if version == 'Line':
                    version = 'Line_Lines'
                self.plot_ten(version, ax3)
                
        if line_rmse:
            fig1.supxlabel('time (s)')
            fig1.supylabel('rmse (m)')
            fig1.tight_layout()

        if ten_rmse:  
            fig2.supylabel("Tension RMS error (kN)")
            fig2.supxlabel("Time (s)")
            fig2.tight_layout()
        
        if plot_ten:
            fig3.suptitle("Line tensions (kN) for "+self.rootname)
            ax3.legend(['v2n','v2o','dev2','v1'])
            ax3.set_ylabel('Tension (kN)')
            ax3.set_xlabel('Time (s)')
            fig3.tight_layout()

        if display: 
            plt.show()
        if save:
            if line_rmse:
                fig1.savefig(self.rootname+'_line_rmse.png', dpi = 300)
            if ten_rmse:
                fig2.savefig(self.rootname+'_ten_rmse.png', dpi = 300)
            if plot_ten:
                fig3.savefig(self.rootname+'_tens.png', dpi = 300)

    def compare_line (self, line_num):

        self.lineList[line_num].loadData(self.out_dirname, self.test, sep = '_') # remember number starts on 1 rather than 0

        class test_line:
            Tdata = self.lineList[line_num].Tdata
            nNodes = self.lineList[line_num].nNodes
            xp = self.lineList[line_num].xp
            yp = self.lineList[line_num].yp
            zp = self.lineList[line_num].zp
            dummy_list = []

        self.lineList[line_num].loadData(self.out_dirname, self.control, sep = '_') # remember number starts on 1 rather than 0
        
        class control_line:
            Tdata = self.lineList[line_num].Tdata
            nNodes = self.lineList[line_num].nNodes
            xp = self.lineList[line_num].xp
            yp = self.lineList[line_num].yp
            zp = self.lineList[line_num].zp
            dummy_list = []
        
        # Work-around for dev2 dtOut flag error
        if len(control_line.Tdata) > 599:
            print("Faulty time steps in dev2. Removed via workaround:")
            print(control_line.Tdata[1], control_line.Tdata[9], control_line.Tdata[34], control_line.Tdata[259])
            i = 0
            while i < 4:
                if i == 0:
                    control_line.dummy_list = list(control_line.Tdata)
                elif i == 1:
                    control_line.dummy_list = list(control_line.xp)
                elif i == 2:
                    control_line.dummy_list = list(control_line.yp)
                elif i == 3:
                    control_line.dummy_list = list(control_line.zp)
                else:
                    break
                del control_line.dummy_list[1]
                del control_line.dummy_list[8]
                del control_line.dummy_list[32]
                del control_line.dummy_list[256]
                del control_line.dummy_list[-1]
                if i ==0:
                    control_line.Tdata = np.array(control_line.dummy_list)
                elif i == 1:
                    control_line.xp = np.array(control_line.dummy_list)
                elif i == 2:
                    control_line.yp = np.array(control_line.dummy_list)
                elif i == 3:
                    control_line.zp = np.array(control_line.dummy_list)
                else:
                    break
                i+=1
            print('Fixed T, xp, yp, and zp values:')
            # print(len(control_line.Tdata), len(control_line.xp), len(control_line.yp), len(control_line.zp))
            print(control_line.Tdata[0], control_line.Tdata[1], control_line.Tdata[2], control_line.Tdata[7], control_line.Tdata[8], control_line.Tdata[9], control_line.Tdata[31], control_line.Tdata[32], control_line.Tdata[33], control_line.Tdata[255], control_line.Tdata[256], control_line.Tdata[257])
            print('')

        # Work-around for v1 dtOut flag error
        if self.test == 'Line':
            print("Faulty time steps in v1. Removed via workaround:")
            print(test_line.Tdata[1], test_line.Tdata[9], test_line.Tdata[34], test_line.Tdata[259])
            i = 0
            while i < 4:
                if i == 0:
                    test_line.dummy_list = list(test_line.Tdata)
                elif i == 1:
                    test_line.dummy_list = list(test_line.xp)
                elif i == 2:
                    test_line.dummy_list = list(test_line.yp)
                elif i == 3:
                    test_line.dummy_list = list(test_line.zp)
                else:
                    break
                del test_line.dummy_list[1]
                del test_line.dummy_list[8]
                del test_line.dummy_list[32]
                del test_line.dummy_list[256]
                del test_line.dummy_list[-1]
                if i ==0:
                    test_line.Tdata = np.array(test_line.dummy_list)
                elif i == 1:
                    test_line.xp = np.array(test_line.dummy_list)
                elif i == 2:
                    test_line.yp = np.array(test_line.dummy_list)
                elif i == 3:
                    test_line.zp = np.array(test_line.dummy_list)
                else:
                    break
                i+=1
            print('Fixed T, xp, yp, and zp values:')
            # print(len(control_line.Tdata), len(control_line.xp), len(control_line.yp), len(control_line.zp))
            print(test_line.Tdata[0], test_line.Tdata[1], test_line.Tdata[2], test_line.Tdata[7], test_line.Tdata[8], test_line.Tdata[9], test_line.Tdata[31], test_line.Tdata[32], test_line.Tdata[33], test_line.Tdata[255], test_line.Tdata[256], test_line.Tdata[257])
            print('')

        if (len(control_line.Tdata) != len(test_line.Tdata)) or (control_line.nNodes != test_line.nNodes):
            print("Error, make sure both datasets are same length or tMax = 600 becasue of bugs with dtOut. Quitting rms error calcualtion...")
            return

        nNodes = control_line.nNodes
        time = control_line.Tdata
        size = (len(time), int(nNodes-1))
        rmse_array = np.zeros(size)

        for j in range(0, nNodes-1):
            i = 0
            for i in range(0, len(time)):

                position1 = np.array([control_line.xp[i,j], control_line.yp[i,j], control_line.zp[i,j]])
                position2 = np.array([test_line.xp[i,j], test_line.yp[i,j], test_line.zp[i,j]])

                diff = np.subtract(position1, position2)

                rmse = np.sqrt(np.sum(np.square(diff))/len(diff))
                rmse_array[i,j] = rmse

        rmse_totals = np.zeros((len(time),2))
        rmse_totals [:,0] = time
        i = 0
        for i in range(0,len(time)):
            rmse_totals[i,1] = np.sqrt(np.sum(np.square(rmse_array[i,:]))/len(rmse_array[i,:]))

        # rmse_totals = np.flip(rmse_totals)

        return rmse_totals

        # Working! Need to find a good way to smooth the data though, becasue the high time resolution means really messy plotting 

    def compare_all_lines(self, ax):

        for i in range(0,len(self.lineList)):
            rmse_line = self.compare_line(i)
            if i == 0:
                time = rmse_line[:,0]
                rmse = np.zeros((len(time),len(self.lineList)))
            rmse [:,i] = rmse_line[:,1]
        
        rmse_final = np.zeros(len(time))
        for j in range(0,len(time)):
            rmse_final[j] = np.sqrt(np.sum(np.square(rmse[j,:]))/len(rmse[j,:]))
        
        ax.plot(time, rmse_final)
        if self.test == 'Line':
            ax.set_title('RMS position error between {} and {}'.format(self.control, self.rootname+'v1'))
        else:
            ax.set_title('RMS position error between {} and {}'.format(self.control, self.test))
        return

    def read_tensions(self, ax):
        data1, ch1, channels1, units1 = read_mooring_file(self.out_dirname, self.control+'.out')
        data2, ch2, channels2, units2 = read_mooring_file(self.out_dirname, self.test+'.out')
        data3 = np.zeros((599,2))
        data4 = np.zeros((599,2))
        # Work-around for dev2 dtOut flag error
        if len(data1[:,0]) > 599:
            print("Faulty time steps in dev2. Removed via workaround:")
            print(data1[1,0], data1[9,0], data1[34,0], data1[259,0])
            i = 0
            while i < 2:
                if i == 0:
                    dummy_list = list(data1[:,0])
                elif i == 1:
                    dummy_list = list(data1[:,1])
                else:
                    break
                del dummy_list[1]
                del dummy_list[8]
                del dummy_list[32]
                del dummy_list[256]
                del dummy_list[-1]
                if i ==0:
                    data3[:,0] = np.array(dummy_list)
                elif i == 1:
                    data3[:,1] = np.array(dummy_list)
                else:
                    break
                i+=1
            print('Fixed T, xp, yp, and zp values:')
            # print(len(control_line.Tdata), len(control_line.xp), len(control_line.yp), len(control_line.zp))
            print(data3[0,0], data3[1,0], data3[2,0], data3[7,0], data3[8,0], data3[9,0], data3[31,0], data3[32,0], data3[33,0], data3[255,0], data3[256,0], data3[257,0])
            print('')
        else: 
            data3 = data1
        if len(data2[:,0]) > 599:
            print("Faulty time steps in dev2. Removed via workaround:")
            print(data2[1,0], data2[9,0], data2[34,0], data2[259,0])
            i = 0
            while i < 2:
                if i == 0:
                    dummy_list = list(data2[:,0])
                elif i == 1:
                    dummy_list = list(data2[:,1])
                else:
                    break
                del dummy_list[1]
                del dummy_list[8]
                del dummy_list[32]
                del dummy_list[256]
                del dummy_list[-1]
                if i ==0:
                    data4[:,0] = np.array(dummy_list)
                elif i == 1:
                    data4[:,1] = np.array(dummy_list)
                else:
                    break
                i+=1
            print('Fixed T, xp, yp, and zp values:')
            # print(len(control_line.Tdata), len(control_line.xp), len(control_line.yp), len(control_line.zp))
            print(data4[0,0], data4[1,0], data4[2,0], data4[7,0], data4[8,0], data4[9,0], data4[31,0], data4[32,0], data4[33,0], data4[255,0], data4[256,0], data4[257,0])
            print('')
        else: 
            data4 = data2

        rmse = np.zeros(np.size(data4[:,0]))
        rmse = np.abs(np.subtract(data3[:,1],data4[:,1]))/1000
        time = data3[:,0]

        ax.plot(time, rmse) 
        ax.set_title("Absolute difference between {} and {}".format(self.test, self.control))
        
    def plot_ten(self, version, ax):
        data1, ch1, channels1, units1 = read_mooring_file(self.out_dirname,version+'.out')
        data3 = np.zeros((599,len(channels1)))
        if len(data1[:,0]) > 599:
            print("Faulty time steps in dev2. Removed via workaround:")
            print(data1[1,0], data1[9,0], data1[34,0], data1[259,0])
            i = 0
            while i < 2:
                if i == 0:
                    dummy_list = list(data1[:,0])
                elif i == 1:
                    dummy_list = list(data1[:,1])
                else:
                    break
                del dummy_list[1]
                del dummy_list[8]
                del dummy_list[32]
                del dummy_list[256]
                del dummy_list[-1]
                if i ==0:
                    data3[:,0] = np.array(dummy_list)
                elif i == 1:
                    data3[:,1] = np.array(dummy_list)
                else:
                    break
                i+=1
            print('Fixed T, xp, yp, and zp values:')
            # print(len(control_line.Tdata), len(control_line.xp), len(control_line.yp), len(control_line.zp))
            print(data3[0,0], data3[1,0], data3[2,0], data3[7,0], data3[8,0], data3[9,0], data3[31,0], data3[32,0], data3[33,0], data3[255,0], data3[256,0], data3[257,0])
            print('')
        else: 
            data3 = data1
        ax.plot(data3[:,0], data3[:,1]/1000)


class run_infile():

    # TODO: Make file names, fig, ax, etc all class variables 

    def __init__(self):
        # TODO: Move all the boolean flags to here as a couple of dictionaries and reduce passing between functions
        pass
    
    def figures(self, instance, path, rootname, animate_all = False, animate_start_end = False, plot_individual_start_end = False, plot_all_start_end = False , plot2d = False, plot3d = False, from_saved_runs = False, display = True, run_v1 = False, run_v2n = False, run_v2o = False, save = False, show_rmse = False, tMax = None):

        if save:
            if not os.path.isdir("saved_runs/figures/"):
                os.system("mkdir saved_runs/figures/")
            if not os.path.isdir("saved_runs/figures/{}".format(rootname)):
                os.system("mkdir saved_runs/figures/{}/".format(rootname))

        if plot_all_start_end or plot_individual_start_end:
            if not (plot2d or plot3d):
                print("For plotting, please specify 2d or 3d. Ending plot function...")
                return

        if animate_all:
            fig1 = plt.figure() 
            ax1 = plt.axes(projection='3d')
            if not from_saved_runs:
                v2oanmi = instance.plot_animation( fig = fig1, ax = ax1, path = path+"Outputs/", rootname = rootname+"v2old", color = "r", start_end_ani = False)
                if run_v2n:
                    v2nanmi = instance.plot_animation( fig = fig1, ax = ax1, path = path+"Outputs/", rootname = rootname+"v2new", color = "y", start_end_ani = False)
                dev2anmi = instance.plot_animation( fig = fig1, ax = ax1, path = path+"Outputs/", rootname = rootname+"dev2", color = "b", start_end_ani = False)
                if run_v1:
                    v1anmi = instance.plot_animation( fig = fig1, ax = ax1, path = path+"Outputs/", rootname = "Line", color = "g", start_end_ani = False)

                print("v2old is red, v2new is yellow, dev2 is blue, and v1 is green")
            else: 
                v2oanmi = instance.plot_animation(fig = fig1, ax = ax1, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2old", color = "r", start_end_ani = False)
                if run_v2n:
                    v2nanmi = instance.plot_animation( fig = fig1, ax = ax1, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2new", color = "y", start_end_ani = False)
                dev2anmi = instance.plot_animation( fig = fig1, ax = ax1, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"dev2", color = "b", start_end_ani = False)
                if run_v1:
                    v1anmi = instance.plot_animation( fig = fig1, ax = ax1, path = "saved_runs/"+rootname+"_outputs/", rootname = "Line", color = "g", start_end_ani = False)
            fig1.legend(["v2old", "v2new","dev2","v1"])
    
        if animate_start_end:
            fig2 = plt.figure() 
            ax2 = plt.axes(projection='3d')

            if not from_saved_runs: 
                v2oanmi = load_infile.plot_animation(instance, fig = fig2, ax = ax2, path = path+"Outputs/", rootname = rootname+"v2old", color = "r", start_end_ani = True)
                if run_v2n:
                    v2nanmi = load_infile.plot_animation(instance, fig = fig2, ax = ax2, path = path+"Outputs/", rootname = rootname+"v2new", color = "y", start_end_ani = True)
                dev2anmi = load_infile.plot_animation(instance, fig = fig2, ax = ax2, path = path+"Outputs/", rootname = rootname+"dev2", color = "b", start_end_ani = True)
                if run_v1:
                    v1anmi = load_infile.plot_animation(instance, fig = fig2, ax = ax2, path = path+"Outputs/", rootname = "Line", color = "g", start_end_ani = True)
            else: 
                v2oanmi = load_infile.plot_animation(instance, fig = fig2, ax = ax2, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2old", color = "r", start_end_ani = True)
                if run_v2n:
                    v2nanmi = load_infile.plot_animation(instance, fig = fig2, ax = ax2, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2new", color = "y", start_end_ani = True)
                dev2anmi = load_infile.plot_animation(instance, fig = fig2, ax = ax2, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"dev2", color = "b", start_end_ani = True)
                if run_v1:
                    v1anmi = load_infile.plot_animation(instance, fig = fig2, ax = ax2, path = "saved_runs/"+rootname+"_outputs/", rootname = "Line", color = "g", start_end_ani = True)
            fig2.legend(["v2old", "v2new","dev2","v1"])

        if plot_all_start_end and plot_individual_start_end:
            if tMax == None:
                print("Please specify tmax for plotting end timestep. Ending plot function...")
                return
            
            if plot2d and plot3d:
                print("Creating individual 2d and 3d plots...")
                threedim = {'projection':'3d'}
                fig3, ax3 = plt.subplots(2,2, figsize=[6,4])
                fig4, ax4 = plt.subplots(2,2, subplot_kw = threedim)
                print("Creating 3d and 2d plot of all versions...")
                fig5, ax5 = plt.subplots(1,2, figsize = [6,2])
                fig6, ax6 = plt.subplots(1,2, subplot_kw = threedim)
                if not from_saved_runs:
                    v2oplot3 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[0][0], path = path+"Outputs/", rootname = rootname+"v2old", md_branch = "v2old", tMax = tMax)
                    v2oplot4 = load_infile.plot_start_end(instance, plot2d = False, fig = fig4, ax = ax4[0][0], path = path+"Outputs/", rootname = rootname+"v2old", md_branch = "v2old", tMax = tMax)
                    v2oplot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = path+"Outputs/", rootname = rootname+"v2old", md_branch = "v2old", color = "r", tMax = tMax)
                    v2oplot2 = load_infile.plot_start_end(instance, plot2d = False, fig = fig6, ax = ax6, path = path+"Outputs/", rootname = rootname+"v2old", md_branch = "v2old", color = "r", tMax = tMax)
                    if run_v2n:
                        v2nplot3 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[0][1], path = path+"Outputs/", rootname = rootname+"v2new", md_branch = "v2new", tMax = tMax)
                        v2nplot4 = load_infile.plot_start_end(instance, plot2d = False, fig = fig4, ax = ax4[0][1], path = path+"Outputs/", rootname = rootname+"v2new", md_branch = "v2new", tMax = tMax)
                        v2nplot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = path+"Outputs/", rootname = rootname+"v2new", md_branch = "v2new", color = "y", tMax = tMax)
                        v2nplot2 = load_infile.plot_start_end(instance, plot2d = False, fig = fig6, ax = ax6, path = path+"Outputs/", rootname = rootname+"v2new", md_branch = "v2new", color = "y", tMax = tMax)
                    dev2plot3 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[1][0], path = path+"Outputs/", rootname = rootname+"dev2", md_branch = "dev2", tMax = tMax)
                    dev2plot4 = load_infile.plot_start_end(instance, plot2d = False, fig = fig4, ax = ax4[1][0], path = path+"Outputs/", rootname = rootname+"dev2", md_branch = "dev2", tMax = tMax)
                    dev2plot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = path+"Outputs/", rootname = rootname+"dev2", md_branch = "dev2", color = "b", tMax = tMax)
                    dev2plot2 = load_infile.plot_start_end(instance, plot2d = False, fig = fig6, ax = ax6, path = path+"Outputs/", rootname = rootname+"dev2", md_branch = "dev2", color = "b", tMax = tMax)
                    if run_v1:
                        v1plot3 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[1][1], path = path+"Outputs/", rootname = "Line", md_branch = "v1", tMax = tMax)
                        v1plot4 = load_infile.plot_start_end(instance, plot2d = False, fig = fig4, ax = ax4[1][1], path = path+"Outputs/", rootname = "Line", md_branch = "v1", tMax = tMax)
                        v1plot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = path+"Outputs/", rootname = "Line", md_branch = "v1", color = "g", tMax = tMax)
                        v1plot2 = load_infile.plot_start_end(instance, plot2d = False, fig = fig6, ax = ax6, path = path+"Outputs/", rootname = "Line", md_branch = "v1", color = "g", tMax = tMax)

                else:
                    v2oplot3 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[0][0], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2old", md_branch = "v2old", tMax = tMax)
                    v2oplot4 = load_infile.plot_start_end(instance, plot2d = False, fig = fig4, ax = ax4[0][0], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2old", md_branch = "v2old", tMax = tMax)
                    v2oplot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2old", md_branch = "v2old", color = "r", tMax = tMax)
                    v2oplot2 = load_infile.plot_start_end(instance, plot2d = False, fig = fig6, ax = ax6, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2old", md_branch = "v2old", color = "r", tMax = tMax)
                    if run_v2n:
                        v2nplot3 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[0][1], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2new", md_branch = "v2new", tMax = tMax)
                        v2nplot4 = load_infile.plot_start_end(instance, plot2d = False, fig = fig4, ax = ax4[0][1], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2new", md_branch = "v2new", tMax = tMax)
                        v2nplot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2new", md_branch = "v2new", color = "y", tMax = tMax)
                        v2nplot2 = load_infile.plot_start_end(instance, plot2d = False, fig = fig6, ax = ax6, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2new", md_branch = "v2new", color = "y", tMax = tMax)
                    dev2plot3 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[1][0], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"dev2", md_branch = "dev2", tMax = tMax)
                    dev2plot4 = load_infile.plot_start_end(instance, plot2d = False, fig = fig4, ax = ax4[1][0], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"dev2", md_branch = "dev2", tMax = tMax)
                    dev2plot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"dev2", md_branch = "dev2", color = "b", tMax = tMax)
                    dev2plot2 = load_infile.plot_start_end(instance, plot2d = False, fig = fig6, ax = ax6, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"dev2", md_branch = "dev2", color = "b", tMax = tMax)
                    if run_v1:
                        v1plot3 =load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[1][1], path = "saved_runs/"+rootname+"_outputs/", rootname = "Line", md_branch = "v1", tMax = tMax)
                        v1plot4 =load_infile.plot_start_end(instance, plot2d = False, fig = fig4, ax = ax4[1][1], path = "saved_runs/"+rootname+"_outputs/", rootname = "Line", md_branch = "v1", tMax = tMax)
                        v1plot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = "saved_runs/"+rootname+"_outputs/", rootname = "Line", md_branch = "v1", color = "g", tMax = tMax)
                        v1plot2 = load_infile.plot_start_end(instance, plot2d = False, fig = fig6, ax = ax6, path = "saved_runs/"+rootname+"_outputs/", rootname = "Line", md_branch = "v1", color = "g", tMax = tMax)

                fig4.tight_layout()
                fig3.tight_layout()
                fig6.tight_layout()
                fig5.tight_layout()
                if save: 
                    fig4.savefig(rootname+"_individual_3dplot.png", dpi = 300)
                    fig3.savefig(rootname+"_individual_2dplot.png", dpi = 300)
                    fig6.savefig(rootname+"_all_3dplot.png", dpi = 300)  
                    fig5.savefig(rootname+"_all_2dplot.png", dpi = 300)

            elif plot2d:
                print("Creating individual 2d plots...")
                fig3, ax3 = plt.subplots(2,2, figsize=[6,4])
                print("Creating 2d plot of all versions...")
                fig5, ax5 = plt.subplots(1,2, figsize = [6,2])
                if not from_saved_runs:
                    v2oplot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[0][0], path = path+"Outputs/", rootname = rootname+"v2old", md_branch = "v2old", tMax = tMax)
                    v2oplot2 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = path+"Outputs/", rootname = rootname+"v2old", md_branch = "v2old", color = "r", tMax = tMax)
                    if run_v2n:
                        v2nplot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[0][1], path = path+"Outputs/", rootname = rootname+"v2new", md_branch = "v2new", tMax = tMax)
                        v2nplot2 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = path+"Outputs/", rootname = rootname+"v2new", md_branch = "v2new", color = "y", tMax = tMax)
                    dev2plot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[1][0], path = path+"Outputs/", rootname = rootname+"dev2", md_branch = "dev2", tMax = tMax)
                    dev2plot2 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = path+"Outputs/", rootname = rootname+"dev2", md_branch = "dev2", color = "b", tMax = tMax)
                    if run_v1:
                        v1plot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[1][1], path = path+"Outputs/", rootname = "Line", md_branch = "v1", tMax = tMax)
                        v1plot2 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = path+"Outputs/", rootname = "Line", md_branch = "v1", color = "g", tMax = tMax)

                else:
                    v2oplot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[0][0], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2old", md_branch = "v2old", tMax = tMax)
                    v2oplot2 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2old", md_branch = "v2old", color = "r", tMax = tMax)
                    if run_v2n:
                        v2nplot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[0][1], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2new", md_branch = "v2new", tMax = tMax)
                        v2nplot2 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2new", md_branch = "v2new", color = "y", tMax = tMax)
                    dev2plot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[1][0], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"dev2", md_branch = "dev2", tMax = tMax)
                    dev2plot2 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"dev2", md_branch = "dev2", color = "b", tMax = tMax)
                    if run_v1:
                        v1plot1 =load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[1][1], path = "saved_runs/"+rootname+"_outputs/", rootname = "Line", md_branch = "v1", tMax = tMax)
                        v1plot2 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = "saved_runs/"+rootname+"_outputs/", rootname = "Line", md_branch = "v1", color = "g", tMax = tMax)

                fig3.tight_layout()
                fig5.tight_layout()
                if save:
                    fig3.savefig(rootname+"_individual_2dplot.png", dpi = 300)
                    fig5.savefig(rootname+"_all_2dplot.png", dpi = 300)
            
            elif plot3d:
                print("Creating individual 3d plots...")
                threedim = {'projection':'3d'}
                fig4, ax4 = plt.subplots(2,2, subplot_kw = threedim)
                print("Creating 3d plot of all versions...")
                fig6, ax6 = plt.subplots(1,2, subplot_kw = threedim)
                if not from_saved_runs:
                    v2oplot1 = load_infile.plot_start_end(instance, fig = fig4, ax = ax4[0][0], path = path+"Outputs/", rootname = rootname+"v2old", md_branch = "v2old", tMax = tMax)
                    v2oplot2 = load_infile.plot_start_end(instance, fig = fig6, ax = ax6, path = path+"Outputs/", rootname = rootname+"v2old", md_branch = "v2old", color = "r", tMax = tMax)
                    if run_v2n:
                        v2nplot1 = load_infile.plot_start_end(instance, fig = fig4, ax = ax4[0][1], path = path+"Outputs/", rootname = rootname+"v2new", md_branch = "v2new", tMax = tMax)
                        v2nplot2 = load_infile.plot_start_end(instance, fig = fig6, ax = ax6, path = path+"Outputs/", rootname = rootname+"v2new", md_branch = "v2new", color = "y", tMax = tMax)
                    dev2plot1 = load_infile.plot_start_end(instance, fig = fig4, ax = ax4[1][0], path = path+"Outputs/", rootname = rootname+"dev2", md_branch = "dev2", tMax = tMax)
                    dev2plot2 = load_infile.plot_start_end(instance, fig = fig6, ax = ax6, path = path+"Outputs/", rootname = rootname+"dev2", md_branch = "dev2", color = "b", tMax = tMax)
                    if run_v1:
                        v1plot1 = load_infile.plot_start_end(instance, fig = fig4, ax = ax4[1][1], path = path+"Outputs/", rootname = "Line", md_branch = "v1", tMax = tMax)
                        v1plot2 = load_infile.plot_start_end(instance, fig = fig6, ax = ax6, path = path+"Outputs/", rootname = "Line", md_branch = "v1", color = "g", tMax = tMax)

                else:
                    v2oplot1 = load_infile.plot_start_end(instance, fig = fig4, ax = ax4[0][0], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2old", md_branch = "v2old", tMax = tMax)
                    v2oplot2 = load_infile.plot_start_end(instance, fig = fig6, ax = ax6, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2old", md_branch = "v2old", color = "r", tMax = tMax)
                    if run_v2n:
                        v2nplot1 = load_infile.plot_start_end(instance, fig = fig4, ax = ax4[0][1], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2new", md_branch = "v2new", tMax = tMax)
                        v2nplot2 = load_infile.plot_start_end(instance, fig = fig6, ax = ax6, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2new", md_branch = "v2new", color = "y", tMax = tMax)
                    dev2plot1 = load_infile.plot_start_end(instance, fig = fig4, ax = ax4[1][0], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"dev2", md_branch = "dev2", tMax = tMax)
                    dev2plot2 = load_infile.plot_start_end(instance, fig = fig6, ax = ax6, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"dev2", md_branch = "dev2", color = "b", tMax = tMax)
                    if run_v1:
                        v1plot1 =load_infile.plot_start_end(instance, fig = fig4, ax = ax4[1][1], path = "saved_runs/"+rootname+"_outputs/", rootname = "Line", md_branch = "v1", tMax = tMax)
                        v1plot2 = load_infile.plot_start_end(instance, fig = fig6, ax = ax6, path = "saved_runs/"+rootname+"_outputs/", rootname = "Line", md_branch = "v1", color = "g", tMax = tMax)

                fig4.tight_layout()
                fig6.tight_layout()
                if save:
                    fig4.savefig(rootname+"_individual_3dplot.png", dpi = 300)
                    fig6.savefig(rootname+"_all_3dplot.png", dpi = 300)
            
            if save:
                if not os.path.isdir("saved_runs/figures/{}".format(rootname)):
                    os.system("mkdir saved_runs/figures/{}".format(rootname))
                os.system("mv {}*.png saved_runs/figures/{}/".format(rootname,rootname))
          
        elif plot_individual_start_end:
            if tMax == None:
                print("Please specify tmax for plotting end timestep. Ending plot function...")
                return

            if plot2d and plot3d:
                print("Creating individual 2d and 3d plots...")
                fig3, ax3 = plt.subplots(2,2, figsize=[6,4])
                threedim = {'projection':'3d'}
                fig4, ax4 = plt.subplots(2,2, subplot_kw = threedim)
                if not from_saved_runs:
                    v2oplot3 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[0][0], path = path+"Outputs/", rootname = rootname+"v2old", md_branch = "v2old", tMax = tMax)
                    v2oplot4 = load_infile.plot_start_end(instance, plot2d = False, fig = fig4, ax = ax4[0][0], path = path+"Outputs/", rootname = rootname+"v2old", md_branch = "v2old", tMax = tMax)
                    if run_v2n:
                        v2nplot3 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[0][1], path = path+"Outputs/", rootname = rootname+"v2new", md_branch = "v2new", tMax = tMax)
                        v2nplot4 = load_infile.plot_start_end(instance, plot2d = False, fig = fig4, ax = ax4[0][1], path = path+"Outputs/", rootname = rootname+"v2new", md_branch = "v2new", tMax = tMax)
                    dev2plot3 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[1][0], path = path+"Outputs/", rootname = rootname+"dev2", md_branch = "dev2", tMax = tMax)
                    dev2plot4 = load_infile.plot_start_end(instance, plot2d = False, fig = fig4, ax = ax4[1][0], path = path+"Outputs/", rootname = rootname+"dev2", md_branch = "dev2", tMax = tMax)
                    if run_v1:
                        v1plot3 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[1][1], path = path+"Outputs/", rootname = "Line", md_branch = "v1", tMax = tMax)
                        v1plot4 = load_infile.plot_start_end(instance, plot2d = False, fig = fig4, ax = ax4[1][1], path = path+"Outputs/", rootname = "Line", md_branch = "v1", tMax = tMax)

                else:
                    v2oplot3 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[0][0], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2old", md_branch = "v2old", tMax = tMax)
                    v2oplot4 = load_infile.plot_start_end(instance, plot2d = False, fig = fig4, ax = ax4[0][0], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2old", md_branch = "v2old", tMax = tMax)
                    if run_v2n:
                        v2nplot3 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[0][1], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2new", md_branch = "v2new", tMax = tMax)
                        v2nplot4 = load_infile.plot_start_end(instance, plot2d = False, fig = fig4, ax = ax4[0][1], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2new", md_branch = "v2new", tMax = tMax)
                    dev2plot3 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[1][0], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"dev2", md_branch = "dev2", tMax = tMax)
                    dev2plot4 = load_infile.plot_start_end(instance, plot2d = False, fig = fig4, ax = ax4[1][0], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"dev2", md_branch = "dev2", tMax = tMax)
                    if run_v1:
                        v1plot3 =load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[1][1], path = "saved_runs/"+rootname+"_outputs/", rootname = "Line", md_branch = "v1", tMax = tMax)
                        v1plot4 =load_infile.plot_start_end(instance, plot2d = False, fig = fig4, ax = ax4[1][1], path = "saved_runs/"+rootname+"_outputs/", rootname = "Line", md_branch = "v1", tMax = tMax)

                fig4.tight_layout()
                fig3.tight_layout()
                if save:
                    fig4.savefig(rootname+"_individual_3dplot.png", dpi = 300)
                    fig3.savefig(rootname+"_individual_2dplot.png", dpi = 300)

            elif plot2d:
                print("Creating individual 2d plots...")
                fig3, ax3 = plt.subplots(2,2, figsize=[6,4])
                if not from_saved_runs:
                    v2oplot = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[0][0], path = path+"Outputs/", rootname = rootname+"v2old", md_branch = "v2old", tMax = tMax)
                    if run_v2n:
                        v2nplot = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[0][1], path = path+"Outputs/", rootname = rootname+"v2new", md_branch = "v2new", tMax = tMax)
                    dev2plot = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[1][0], path = path+"Outputs/", rootname = rootname+"dev2", md_branch = "dev2", tMax = tMax)
                    if run_v1:
                        v1plot = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[1][1], path = path+"Outputs/", rootname = "Line", md_branch = "v1", tMax = tMax)

                else:
                    v2oplot = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[0][0], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2old", md_branch = "v2old", tMax = tMax)
                    if run_v2n:
                        v2nplot = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[0][1], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2new", md_branch = "v2new", tMax = tMax)
                    dev2plot = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[1][0], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"dev2", md_branch = "dev2", tMax = tMax)
                    if run_v1:
                        v1plot =load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig3, ax = ax3[1][1], path = "saved_runs/"+rootname+"_outputs/", rootname = "Line", md_branch = "v1", tMax = tMax)
                
                fig3.tight_layout()
                if save:
                    fig3.savefig(rootname+"_individual_2dplot.png", dpi = 300)
            
            elif plot3d:
                print("Creating individual 3d plots...")
                threedim = {'projection':'3d'}
                fig4, ax4 = plt.subplots(2,2, subplot_kw = threedim)
                if not from_saved_runs:
                    v2oplot = load_infile.plot_start_end(instance, fig = fig4, ax = ax4[0][0], path = path+"Outputs/", rootname = rootname+"v2old", md_branch = "v2old", tMax = tMax)
                    if run_v2n:
                        v2nplot = load_infile.plot_start_end(instance, fig = fig4, ax = ax4[0][1], path = path+"Outputs/", rootname = rootname+"v2new", md_branch = "v2new", tMax = tMax)
                    dev2plot = load_infile.plot_start_end(instance, fig = fig4, ax = ax4[1][0], path = path+"Outputs/", rootname = rootname+"dev2", md_branch = "dev2", tMax = tMax)
                    if run_v1:
                        v1plot = load_infile.plot_start_end(instance, fig = fig4, ax = ax4[1][1], path = path+"Outputs/", rootname = "Line", md_branch = "v1", tMax = tMax)

                else:
                    v2oplot = load_infile.plot_start_end(instance, fig = fig4, ax = ax4[0][0], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2old", md_branch = "v2old", tMax = tMax)
                    if run_v2n:
                        v2nplot = load_infile.plot_start_end(instance, fig = fig4, ax = ax4[0][1], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2new", md_branch = "v2new", tMax = tMax)
                    dev2plot = load_infile.plot_start_end(instance, fig = fig4, ax = ax4[1][0], path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"dev2", md_branch = "dev2", tMax = tMax)
                    if run_v1:
                        v1plot =load_infile.plot_start_end(instance, fig = fig4, ax = ax4[1][1], path = "saved_runs/"+rootname+"_outputs/", rootname = "Line", md_branch = "v1", tMax = tMax)
                fig4.tight_layout()
                if save:
                    fig4.savefig(rootname+"_individual_3dplot.png", dpi = 300)
            
        elif plot_all_start_end: 
            if tMax == None:
                print("Please specify tmax for plotting end timestep. Ending plot function...")
                return

            if plot2d and plot3d:
                print("Creating 3d and 2d plot of all versions...")
                threedim = {'projection':'3d'}
                fig5, ax5 = plt.subplots(1,2, figsize = [6,2])
                fig6, ax6 = plt.subplots(1,2, subplot_kw = threedim)
                if not from_saved_runs:
                    v2oplot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = path+"Outputs/", rootname = rootname+"v2old", md_branch = "v2old", color = "r", tMax = tMax)
                    v2oplot2 = load_infile.plot_start_end(instance, plot2d = False, fig = fig6, ax = ax6, path = path+"Outputs/", rootname = rootname+"v2old", md_branch = "v2old", color = "r", tMax = tMax)
                    if run_v2n:
                        v2nplot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = path+"Outputs/", rootname = rootname+"v2new", md_branch = "v2new", color = "y", tMax = tMax)
                        v2nplot2 = load_infile.plot_start_end(instance, plot2d = False, fig = fig6, ax = ax6, path = path+"Outputs/", rootname = rootname+"v2new", md_branch = "v2new", color = "y", tMax = tMax)
                    dev2plot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = path+"Outputs/", rootname = rootname+"dev2", md_branch = "dev2", color = "b", tMax = tMax)
                    dev2plot2 = load_infile.plot_start_end(instance, plot2d = False, fig = fig6, ax = ax6, path = path+"Outputs/", rootname = rootname+"dev2", md_branch = "dev2", color = "b", tMax = tMax)
                    if run_v1:
                        v1plot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = path+"Outputs/", rootname = "Line", md_branch = "v1", color = "g", tMax = tMax)
                        v1plot2 = load_infile.plot_start_end(instance, plot2d = False, fig = fig6, ax = ax6, path = path+"Outputs/", rootname = "Line", md_branch = "v1", color = "g", tMax = tMax)

                else:
                    v2oplot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2old", md_branch = "v2old", color = "r", tMax = tMax)
                    v2oplot2 = load_infile.plot_start_end(instance, plot2d = False, fig = fig6, ax = ax6, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2old", md_branch = "v2old", color = "r", tMax = tMax)
                    if run_v2n:
                        v2nplot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2new", md_branch = "v2new", color = "y", tMax = tMax)
                        v2nplot2 = load_infile.plot_start_end(instance, plot2d = False, fig = fig6, ax = ax6, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2new", md_branch = "v2new", color = "y", tMax = tMax)
                    dev2plot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"dev2", md_branch = "dev2", color = "b", tMax = tMax)
                    dev2plot2 = load_infile.plot_start_end(instance, plot2d = False, fig = fig6, ax = ax6, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"dev2", md_branch = "dev2", color = "b", tMax = tMax)
                    if run_v1:
                        v1plot1 = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = "saved_runs/"+rootname+"_outputs/", rootname = "Line", md_branch = "v1", color = "g", tMax = tMax)
                        v1plot2 = load_infile.plot_start_end(instance, plot2d = False, fig = fig6, ax = ax6, path = "saved_runs/"+rootname+"_outputs/", rootname = "Line", md_branch = "v1", color = "g", tMax = tMax)

                fig6.tight_layout()
                fig5.tight_layout()
                if save:
                    fig6.savefig(rootname+"_all_3dplot.png", dpi = 300)  
                    fig5.savefig(rootname+"_all_2dplot.png", dpi = 300)
            
            elif plot2d: 
                print("Creating 2d plot of all versions...")
                fig5, ax5 = plt.subplots(1,2, figsize = [6,2])
                if not from_saved_runs:
                    v2oplot = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = path+"Outputs/", rootname = rootname+"v2old", md_branch = "v2old", color = "r", tMax = tMax)
                    if run_v2n:
                        v2nplot = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = path+"Outputs/", rootname = rootname+"v2new", md_branch = "v2new", color = "y", tMax = tMax)
                    dev2plot = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = path+"Outputs/", rootname = rootname+"dev2", md_branch = "dev2", color = "b", tMax = tMax)
                    if run_v1:
                        v1plot = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = path+"Outputs/", rootname = "Line", md_branch = "v1", color = "g", tMax = tMax)

                else:
                    v2oplot = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2old", md_branch = "v2old", color = "r", tMax = tMax)
                    if run_v2n:
                        v2nplot = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2new", md_branch = "v2new", color = "y", tMax = tMax)
                    dev2plot = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"dev2", md_branch = "dev2", color = "b", tMax = tMax)
                    if run_v1:
                        v1plot = load_infile.plot_start_end(instance, plot2d = plot2d, fig = fig5, ax = ax5, path = "saved_runs/"+rootname+"_outputs/", rootname = "Line", md_branch = "v1", color = "g", tMax = tMax)

                fig5.tight_layout()
                if save:
                    fig5.savefig(rootname+"_all_2dplot.png", dpi = 300)

            elif plot3d:
                print("Creating 3d plot of all versions...")
                threedim = {'projection':'3d'}
                fig6, ax6 = plt.subplots(1,2, subplot_kw = threedim)
                if not from_saved_runs:
                    v2oplot = load_infile.plot_start_end(instance, fig = fig6, ax = ax6, path = path+"Outputs/", rootname = rootname+"v2old", md_branch = "v2old", color = "r", tMax = tMax)
                    if run_v2n:
                        v2nplot = load_infile.plot_start_end(instance, fig = fig6, ax = ax6, path = path+"Outputs/", rootname = rootname+"v2new", md_branch = "v2new", color = "y", tMax = tMax)
                    dev2plot = load_infile.plot_start_end(instance, fig = fig6, ax = ax6, path = path+"Outputs/", rootname = rootname+"dev2", md_branch = "dev2", color = "b", tMax = tMax)
                    if run_v1:
                        v1plot = load_infile.plot_start_end(instance, fig = fig6, ax = ax6, path = path+"Outputs/", rootname = "Line", md_branch = "v1", color = "g", tMax = tMax)

                else:
                    v2oplot = load_infile.plot_start_end(instance, fig = fig6, ax = ax6, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2old", md_branch = "v2old", color = "r", tMax = tMax)
                    if run_v2n:
                        v2nplot = load_infile.plot_start_end(instance, fig = fig6, ax = ax6, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"v2new", md_branch = "v2new", color = "y", tMax = tMax)
                    dev2plot = load_infile.plot_start_end(instance, fig = fig6, ax = ax6, path = "saved_runs/"+rootname+"_outputs/", rootname = rootname+"dev2", md_branch = "dev2", color = "b", tMax = tMax)
                    if run_v1:
                        v1plot = load_infile.plot_start_end(instance, fig = fig6, ax = ax6, path = "saved_runs/"+rootname+"_outputs/", rootname = "Line", md_branch = "v1", color = "g", tMax = tMax)
                fig6.tight_layout()
                if save:
                    fig6.savefig(rootname+"_all_3dplot.png", dpi = 300)  
        
        if show_rmse:
            if run_v2o or run_v1 or run_v2n:
                if from_saved_runs:
                    location = "saved_runs/"+rootname+"_outputs/"
                    load_infile.rms_comparison(instance, display = display, save = save, run_v1 = run_v1, run_v2n= run_v2n, run_v2o = run_v2o, path = location)
                else:
                    location = path+"Outputs/"
                    load_infile.rms_comparison(instance, display = display, save = save, run_v1 = run_v1, run_v2n= run_v2n, run_v2o = run_v2o, path = location)
            else:
                print("RMS error calcaulation requires at least two datasets. Ending plot function...")
                return
            
        if not (show_rmse or plot_all_start_end or plot_individual_start_end or animate_all or animate_start_end):
            print("Spececify an option to make an annimation or plot. Ending plot function...")
            return

        if display:
            print("Displaying figure...")
            plt.show()
       
        if save:
            os.system("mv {}*.png saved_runs/figures/{}/".format(rootname,rootname))
            print("{}_individual_*.png saved to saved_runs/figures/".format(rootname))

        return

    def save_outputs(self, remove_version_inputs = True, del_logs = True):
        os.system("mkdir saved_runs/")
        os.system("mkdir saved_runs/{}_outputs/".format(self.rootname)) # TODO: Add a 'if is directory, then skip'
        os.system("mv MooringTest/*.out saved_runs/{}_outputs/".format(self.rootname))
        os.system("mv MooringTest/*.log saved_runs/{}_outputs/".format(self.rootname))
        os.system("mv MooringTest/Outputs/* saved_runs/{}_outputs/".format(self.rootname))
        os.system("./clean_outputs")
        if del_logs:
            os.system("rm saved_runs/{}_outputs/*.log".format(self.rootname))

        print("-------------------------")
        print("Output files from versions saved to saved_runs/{}_outputs/".format(self.rootname))
        print("Files of the form Line_Line*.out are v1 output files")

        if remove_version_inputs:
            os.system("rm {}".format(self.path+self.rootname+"v2old"+self.extension))
            os.system("rm {}".format(self.path+self.rootname+"v2new"+self.extension))
            os.system("rm {}".format(self.path+self.rootname+"dev2"+self.extension))

    def get_positions(self):
        scaler = [1., 1., 1., np.pi/180., np.pi/180., np.pi/180.]  # for scaling platform position inputs
        outFileName = "PtfmMotions.dat"
        i=0  # file line number
        t_in = []
        Xp_in = []
        myfile2 = open(outFileName, 'r');     # open an input stream to the line data input file
        if myfile2:
            print(outFileName, " opened.") 
            for line2 in myfile2:
                # skip data in first two lines (headers)
                if (i < 2):
                    i+=1
                    continue
                
                #  split line by tabs
                datarow = list(line2.split())
                
                if (len(datarow) < 7): 
                    print("Seems like we've hit a bad line or end of file. ")
                    break;                  # break if we're on a last empty line
                
                t_in.append(float(datarow[0]))
                scaled_data = [0,0,0,0,0,0]			
                for j in range(5):
                    scaled_data[j] = float(datarow[j+1])*scaler[j] # add platform positions
                Xp_in.append(scaled_data)
                i += 1
        else:
            print("ERROR: could not open ", outFileName)
            return

        myfile2.close()
        print("Done reading PtfmMotions.dat. Last line read: ", i)

        Xp_in = np.array(Xp_in)

        # interpolator for platform positions: t_in is vector of time steps from position input file. xp_in is dof
        ts = 0
        self.xp = np.zeros((len(self.time),len(Xp_in[0])))
        for its in range(0, len(self.time)):

            t = its*self.dtC
            
            # interpolate platform positions from .out file data, and approximate velocities
            while ts < (len(t_in)-1):  # search through platform motion data time steps (from .out file)	
                if (t_in[ts+1] > t):				
                    frac = ( t - t_in[ts] )/( t_in[ts+1] - t_in[ts] )		# get interpolation fraction
                    for j in range(0, len(Xp_in[0])):
                        self.xp[its][j] = Xp_in[ts][j] + frac*( Xp_in[ts+1][j] - Xp_in[ts][j] ) # interpolate for each platform DOF
                    break
                ts += 1

        self.xdp = np.zeros((len(self.time),len(Xp_in[0])))
        xold = np.zeros(len(Xp_in[0]))
        # calculate velocities using finite difference
        for i in range(len(self.time)):
            self.xdp [i] = (self.xp[i] - xold)/self.dtC
            xold =  self.xp[i]
        
        return
        
    def sin (self, axis = 0):
        # axis 0 -> x, 1 -> y, 3 -> z
        self.xp = np.zeros((len(self.time),6))
        
        # Harmonic properties
        omega = (2*np.pi)/(len(self.time)/2)
        A = 50

        for i in range(len(self.time)):
            self.xp[i,axis] = A * np.sin(i*omega)

        self.xdp = np.zeros((len(self.time),6))
        xold = np.zeros(6)
        # calculate velocities using finite difference
        for i in range(len(self.time)):
            self.xdp [i] = (self.xp[i] - xold)/self.dtC
            xold =  self.xp[i]

    def load_dynamics(self, static, from_file = False, y_sin = True):

       # initializing
        self.time = np.arange(0, self.tMax, self.dtC)
        size = (len(self.time), self.vector_size)
        self.x = np.zeros(size)
        self.xd = np.zeros(size)

        if static:
            self.xdp = np.zeros(len(self.time),6)
            self.xp = np.zeros(len(self.time),6)
            for i in range(len(self.time)):
                self.x[i,:] = self.xi

        elif from_file:
            self.get_positions()
            for i in range(len(self.time)):
                if i == 0:
                    self.x[i,:] = self.xi
                else:
                    j = 0
                    while j < self.vector_size:
                        self.x[i,j:j+self.dof] = self.x[i-1,j:j+self.dof] + self.xdp[i, 0:self.dof] * self.dtC
                        self.xd[i,j:j+self.dof] = self.xdp[i, 0:self.dof]
                        j += self.dof
        elif y_sin:
            self.sin(axis = 1) 
            for i in range(len(self.time)):
                if i == 0:
                    self.x[i,:] = self.xi
                else:
                    j = 0
                    while j < self.vector_size:
                        self.x[i,j:j+self.dof] = self.x[i-1,j:j+self.dof] + self.xdp[i, 0:self.dof] * self.dtC
                        self.xd[i,j:j+self.dof] = self.xdp[i, 0:self.dof]
                        j += self.dof

        return

    def run_v1 (self, dylib = None): 
        #Double vector pointer data type
        double_p = ctypes.POINTER(ctypes.c_double)
        # -------------------- load the MoorDyn DLL ---------------------

        # Make MoorDyn function prototypes and parameter lists (remember, first entry is return type, rest are args)
        MDInitProto = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(ctypes.c_double*6), ctypes.POINTER(ctypes.c_double*6)) #need to add filename option here, maybe this c_char works? #need to determine char size 
        MDStepProto = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(ctypes.c_double*6), ctypes.POINTER(ctypes.c_double*6), ctypes.POINTER(ctypes.c_double*6), double_p, double_p)
        MDClosProto = ctypes.CFUNCTYPE(ctypes.c_int)

        MDInitParams = (1, "x"), (1, "xd")
        MDStepParams = (1, "x"), (1, "xd"), (2, "f"), (1, "t"), (1, "dtC") 

        if dylib == None:
            dylib = 'compileDYLIB/MoorDyn.dylib'
        MDdylib = ctypes.CDLL(dylib) #load moordyn dylib

        MDInit = MDInitProto(("LinesInit", MDdylib), MDInitParams)
        MDStep = MDStepProto(("LinesCalc", MDdylib), MDStepParams)
        MDClose= MDClosProto(("LinesClose", MDdylib))
        # ------------------------ run MoorDyn ---------------------------
        print("==================================================")    
        # initialize some arrays for communicating with MoorDyn
        t  = double_p()    # pointer to t
        # Converting to ctypes
        dtC = ctypes.pointer(ctypes.c_double(self.dtC))

        # initialize MoorDyn at origin
        MDInit((self.xp[0,:]).ctypes.data_as(ctypes.POINTER(ctypes.c_double*6)),(self.xdp[0,:]).ctypes.data_as(ctypes.POINTER(ctypes.c_double*6)))

        # loop through coupling time steps
        print("MoorDyn initialized - now performing calls to MoorDynStep...")
        for i in range(len(self.time)):
            t = ctypes.pointer(ctypes.c_double(self.time[i]))
            MDStep((self.xp[i,:]).ctypes.data_as(ctypes.POINTER(ctypes.c_double*6)), (self.xdp[i,:]).ctypes.data_as(ctypes.POINTER(ctypes.c_double*6)), t, dtC)    

        print("Succesffuly simulated for {} seconds - now closing MoorDyn...".format(self.tMax))  
            
        # close MoorDyn simulation (clean up the internal memory, hopefully) when finished
        MDClose()   
        print("v1 script executed successfully")
        print("++++++++++++++++++++++++++++++++++++++++++++++++++")
        del MDdylib

    def run_old_API (self, version = None, dylib = None): 
        print("==================================================")
        # -------------------- load the MoorDyn DLL ---------------------

        #Double vector pointer data type
        double_p = ctypes.POINTER(ctypes.c_double)

        # Make MoorDyn function prototypes and parameter lists (remember, first entry is return type, rest are args)
        MDInitProto = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(ctypes.c_double*self.vector_size), ctypes.POINTER(ctypes.c_double*self.vector_size), ctypes.c_char_p) #need to add filename option here, maybe this c_char works? #need to determine char size 
        MDStepProto = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.POINTER(ctypes.c_double*self.vector_size), ctypes.POINTER(ctypes.c_double*self.vector_size), ctypes.POINTER(ctypes.c_double*self.vector_size), double_p, double_p)
        MDClosProto = ctypes.CFUNCTYPE(ctypes.c_int)

        MDInitParams = (1, "x"), (1, "xd"), (1, "infilename") #guessing the 1 flag here means input?
        MDStepParams = (1, "x"), (1, "xd"), (2, "f"), (1, "t"), (1, "dtC") 

        if dylib == None:
            if version == "dev2":
                dylib_path = "dev2_DYLIB/MoorDyn2.dylib"
            elif version == "v2old":
                dylib_path = "v2_DYLIB/libmoordyn2.dylib"
            else:
                print("Please specify dev2 or v2 as versions for old_API run")
                return
        else:
            dylib_path = dylib

        print("dylib path is ", dylib_path)
        MDdylib = ctypes.CDLL(dylib_path) #load moordyn dylib

        filename = self.path+self.rootname+version+self.extension

        MDInit = MDInitProto(("MoorDynInit", MDdylib), MDInitParams)
        MDStep = MDStepProto(("MoorDynStep", MDdylib), MDStepParams)
        MDClose= MDClosProto(("MoorDynClose", MDdylib))  
    # ------------------------ run MoorDyn ---------------------------
        # initialize some arrays for communicating with MoorDyn
        t  = double_p()    # pointer to t

        # parameters
        dtC = ctypes.pointer(ctypes.c_double(self.dtC))

        infile = ctypes.c_char_p(bytes(filename, encoding='utf8'))

        # initialize MoorDyn at origin
        MDInit((self.x[0,:]).ctypes.data_as(ctypes.POINTER(ctypes.c_double*self.vector_size)),(self.xd[0,:]).ctypes.data_as(ctypes.POINTER(ctypes.c_double*self.vector_size)),infile)
        print("MoorDyn initialized - now performing calls to MoorDynStep...")

        # loop through coupling time steps
        for i in range(len(self.time)):
            t = ctypes.pointer(ctypes.c_double(self.time[i]))
            MDStep((self.x[i,:]).ctypes.data_as(ctypes.POINTER(ctypes.c_double*self.vector_size)), (self.xd[i,:]).ctypes.data_as(ctypes.POINTER(ctypes.c_double*self.vector_size)), t, dtC)    
        print("Succesffuly simulated for {} seconds - now closing MoorDyn...".format(self.tMax))  

        # close MoorDyn simulation (clean up the internal memory, hopefully) when finished
        MDClose()   
        print("Old API {} script executed successfully".format(version))
        print("++++++++++++++++++++++++++++++++++++++++++++++++++")
        del MDdylib

    def run_v2new (self): 

        print("==================================================")
        print("This runs the python wrapper of MoorDynV2, it does not reference the local MoorDyn copy that is being edited in ../MoorDyn")
        system = moordyn.Create(self.path+self.rootname+"v2new"+self.extension)
        moordyn.Init(system, self.x[0,:], self.xd[0,:])
        # loop through coupling time steps
        print("MoorDyn initialized - now performing calls to MoorDynStep...")
        for i in range(len(self.time)):
            # call the MoorDyn step function
            moordyn.Step(system, self.x[i,:], self.xd[i,:], self.time[i], self.dtC)    #force value returned here in array

        print("Successfuly simulated for {} seconds - now closing MoorDyn...".format(self.tMax))  

        # close MoorDyn simulation (clean up the internal memory, hopefully) when finished
        moordyn.Close(system)   

        print("New API v2 script executed successfully")
        print("++++++++++++++++++++++++++++++++++++++++++++++++++")
        del system

    def simulate_all (self, run_v1, run_v2n, run_v2o):
        print("System initialized, now simulating versions of MoorDyn")
        if os.path.isfile("MooringTest/*.out") or os.path.isfile("Mooring/*.out"):
            os.system('rm Mooring*/*.out')
            os.system('echo "Removed files with .out extension from Mooring*"')

        if os.path.isfile("MooringTest/*.log") or os.path.isfile("Mooring/*.log"):
            os.system('rm Mooring*/*.log')
            os.system('echo "Removed files with .log extension from Mooring*"')

        if len(os.listdir("MooringTest/Outputs")) == 0 :
            os.system('rm MooringTest/Outputs/*')
            os.system('echo "Removed files with .out and .log extension from MooringTest/Outputs/"')
        os.system("cp {} {}".format(self.path+self.rootname+self.extension, self.path+self.rootname+"v2old"+self.extension))
        os.system("cp {} {}".format(self.path+self.rootname+self.extension, self.path+self.rootname+"dev2"+self.extension))

        self.run_old_API(version = "dev2")
       
        if run_v2o:
            self.run_old_API(version = "v2old")

        if run_v2n:
            os.system("cp {} {}".format(self.path+self.rootname+self.extension, self.path+self.rootname+"v2new"+self.extension))
            self.run_v2new()

        if run_v1: 
            self.run_v1()
            os.system("./namechange_v1")

        os.system("mv Mooring/Line_Line* MooringTest/Outputs/")
        os.system("mv MooringTest/*.out MooringTest/Outputs/")
        os.system("mv MooringTest/*.log MooringTest/Outputs/")

    def run(self, rootname, extension, path, tMax, dof, debug = False, run_v1 = False, run_v2n = False, run_v2o = False, simulate = False, plot = False, static = False, fig_options = None):
        #------------------- Set up Mooring line conditions -----------------------------

        print("Loading input data...")
        self.path = path
        self.rootname = rootname
        self.extension = extension
        self.dof = int(dof)
        from_saved_runs = True

        if from_saved_runs:
            out_dirname = "saved_runs/"+rootname+"_outputs/"
        else:
            out_dirname = path+"Outputs/"

        inputs = load_infile(in_dirname = path, out_dirname = out_dirname, rootname = rootname, extension = extension)

        self.vector_size = int(inputs.numfair*self.dof)

        # parameters
        self.dtC = float(inputs.MDoptions["dtM"])
        self.tMax = tMax

        # Inital fairlead locations
        i=0
        self.xi = np.zeros(self.vector_size)
        for point in inputs.pointList:
            if point.type == -1:  
                self.xi[i]=(point.r[0])
                self.xi[i+1]=(point.r[1])
                self.xi[i+2]=(point.r[2])
                i += dof

        self.load_dynamics(static)
        
        print("------------------------------------------------------")
        
        if not debug: 
            if simulate: 
                if run_v1:
                    print("Converting v2 input into v1...")
                    inputs.v1_build(outfile = "Mooring/lines.txt")

                self.simulate_all (run_v1, run_v2n, run_v2o)
                
                if plot:
                    if fig_options['save_fig']:
                        if not os.path.isdir("saved_runs/"):
                            os.system("mkdir saved_runs/")
                    self.figures(inputs, path, rootname, animate_all = fig_options['animate_all'], animate_start_end =  fig_options['animate_start_end'], plot_individual_start_end = fig_options['plot_individual_start_end'], plot_all_start_end = fig_options['plot_all_start_end'], plot2d = fig_options['plot2d'], plot3d = fig_options['plot3d'], save = fig_options['save_fig'], from_saved_runs = False, display = fig_options['display'], run_v1 = run_v1, run_v2n = run_v2n, run_v2o = run_v2o, show_rmse = fig_options['show_rmse'], tMax = self.tMax-1)

                self.save_outputs (remove_version_inputs = True, del_logs = True)

            elif plot:
                self.figures(inputs, path, rootname, animate_all = fig_options['animate_all'], animate_start_end =  fig_options['animate_start_end'], plot_individual_start_end = fig_options['plot_individual_start_end'], plot_all_start_end = fig_options['plot_all_start_end'], plot2d = fig_options['plot2d'], plot3d = fig_options['plot3d'], save = fig_options['save_fig'], from_saved_runs = True, display = fig_options['display'], run_v1 = run_v1, run_v2n = run_v2n, run_v2o = run_v2o, show_rmse = fig_options['show_rmse'], tMax = self.tMax)

            else:
                print("Please specify appropriate boolean flags. Quitting...")
                exit()
        
        else: # debugging space
            inputs.rms_comparison (control = "dev2", display = True, save = False, run_v1 = True, run_v2n = True, run_v2o = True, line_rmse = False, ten_rmse = True, path = "MooringTest/Outputs/")

        plt.close('all')

        print("------------end--------------")

if __name__ == "__main__":
    #------------------- Run All Scripts -----------------------------

    # Flags for running
   
    run_v1 = True
    run_v2n = True
    run_v2o = True
    static = False

    debug = False
    simulate = True
    plot = True

    fig_options = {}
    if plot:
    # Note: This only apply if plot boolean is true
    # Flags for plotting: can only plot indidvual or all, not both. If both are true, only individual runs
        fig_options = {'animate_all' : True, 'animate_start_end' : False, 'plot2d' : False, 'plot3d' : False, 'plot_individual_start_end' : False, 'plot_all_start_end' : False, 'display' : True, 'save_fig' : False, 'show_rmse' : True}
    
    #TODO: figure out why plot_animation and plot_start_end are giving different pictures for the final timestep
        # plot_start_end is the correct output, verified by MoorPy
        # plot_animation is wierd when plotting multiple files. Only the last one is animated, all others are stuck at time step 0

    rootname = "cable" 
    extension = ".dat"
    path = "MooringTest/"

    tMax = 600.0  # simulation duration (s)
    
    dof = 3 # Size of X and XD vectors: 3 DOF for lines, points, connections, 6 DOF for bodies and rods. Ex for three points, size should be 9. 

    instance = run_infile()
    run_infile.run(instance, rootname = rootname, extension = extension, path = path, tMax = tMax, dof = dof, debug = debug, run_v1 = run_v1, run_v2n = run_v2n, run_v2o = run_v2o, simulate = simulate, plot = plot, fig_options = fig_options, static = static)