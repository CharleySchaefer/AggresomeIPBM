import sys
import numpy as np
import math
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
import os
from itertools import product, combinations
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from mpl_toolkits import mplot3d
from matplotlib.widgets import Slider, Button, RadioButtons
import sys
import pandas as pd
import time


#-------------------------------------------------------------------------------------------------------------
# GIVES EQUAL RATIO FOR AXES IN 3D PLOT
# https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
#-------------------------------------------------------------------------------------------------------------
def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
      
    Source: From "Karlo"
    https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
    
    
    
    
    
#-------------------------------------------------------------------------------------------------------------    
#DRAWS MY CELL BOX with microtubuli if wanted
#-------------------------------------------------------------------------------------------------------------
def draw_box(ax, s_p, microtub, fixed_aggresomes, X0, NX, Y0, NY, Z0, NZ, xnuL, xnuU, ynuL, ynuU, znuL, znuU):
    """Draws the cubic cell box with desired structures using the current matplotlib backend.
    
    The cubic cell can be drawn with or without structures and with or without the microtubule.
    The convention of the arguments X0,NX etc. is inconsistend compared to the one used for
    xnuL, xnU.
    X0 is part of the cell, NX is not. 
    xnuL is part of the cell, xnuU is part of the cell.
    
    Args:
        ax (Axes): Instance of Axes class from matplotlib that box is drawn on.
        s_p (bool): Single phase parameter. If True, a cell with only cytoplasm is 
            drawn.
        fixed_aggresomes (bool): Determines if fixed cubic aggresomes are drawn.
        microtub (bool): Microtubule parameter. If True, microtubule are drawn.
        X0 (int): x-value of the lower boundary of the cell. (X0,y,z) must be part of cell.
        NX (int): x-value of the upper boundary of the cell. (NX,y,z) must NOT be part of cell.
        ...
        xnuL (int): x-value of the lower boundary of the nucleus. (xnuL,ynuL,znuL) is part of nucleus.
        xnuU (int): x-value of the upper boundary of the nucleus. (xnuU,ynuU,znuU) is part of nucleus too.
        ...
        
    """
    
    
    if(s_p == False):
        points_nu = np.array(list(product([xnuL-0.5,xnuU+0.5],[ynuL-0.5,ynuU+0.5],[znuL-0.5,znuU+0.5])))
        
        for s,e in combinations(points_nu, 2):
            if np.sum(np.abs(s-e)) in {xnuU-xnuL+1,ynuU-ynuL+1,znuU-znuL+1}:
                ax.plot3D(*zip(s, e), color="black", alpha=0.25)
            
        sides_nu = [[[xnuL-0.5,ynuL-0.5,znuL-0.5],[xnuU+0.5,ynuL-0.5,znuL-0.5],[xnuU+0.5,ynuU+0.5,znuL-0.5],[xnuL-0.5,ynuU+0.5,znuL-0.5]],
                [[xnuL-0.5,ynuL-0.5,znuL-0.5],[xnuU+0.5,ynuL-0.5,znuL-0.5],[xnuU+0.5,ynuL-0.5,znuU+0.5],[xnuL-0.5,ynuL-0.5,znuU+0.5]],
                 [[xnuU+0.5,ynuL-0.5,znuL-0.5],[xnuU+0.5,ynuU+0.5,znuL-0.5],[xnuU+0.5,ynuU+0.5,znuU+0.5],[xnuU+0.5,ynuL-0.5,znuU+0.5]],
                 [[xnuU+0.5,ynuU+0.5,znuL-0.5],[xnuL-0.5,ynuU+0.5,znuL-0.5],[xnuL-0.5,ynuU+0.5,znuU+0.5],[xnuU+0.5,ynuU+0.5,znuU+0.5]],
                 [[xnuL-0.5,ynuL-0.5,znuL-0.5],[xnuL-0.5,ynuU+0.5,znuL-0.5],[xnuL-0.5,ynuU+0.5,znuU+0.5],[xnuL-0.5,ynuL-0.5,znuU+0.5]],
                 [[xnuL-0.5,ynuL-0.5,znuU+0.5],[xnuU+0.5,ynuL-0.5,znuU+0.5],[xnuU+0.5,ynuU+0.5,znuU+0.5],[xnuL-0.5,ynuU+0.5,znuU+0.5]]
                ]
        
        collection_nu = Poly3DCollection(sides_nu, alpha=.1)
        
        
        collection_nu.set_facecolor("green")
        
        ax.add_collection3d(collection_nu)
    
    points = np.array(list(product([X0-0.5,NX-0.5],[Y0-0.5,NY-0.5],[Z0-0.5,NZ-0.5])))
    

    for s,e in combinations(points, 2):
        if np.sum(np.abs(s-e)) in {NX-X0, NY-Y0, NZ-Z0}:
            ax.plot3D(*zip(s, e), color="black", alpha=0.25)

            

    sides = [[[X0-0.5,Y0-0.5,Z0-0.5],[NX-0.5,Y0-0.5,Z0-0.5],[NX-0.5,NY-0.5,Z0-0.5],[X0-0.5,NY-0.5,Z0-0.5]],
            [[X0-0.5,Y0-0.5,Z0-0.5],[NX-0.5,Y0-0.5,Z0-0.5],[NX-0.5,Y0-0.5,NZ-0.5],[X0-0.5,Y0-0.5,NZ-0.5]],
             [[NX-0.5,Y0-0.5,Z0-0.5],[NX-0.5,NY-0.5,Z0-0.5],[NX-0.5,NY-0.5,NZ-0.5],[NX-0.5,Y0-0.5,NZ-0.5]],
             [[NX-0.5,NY-0.5,Z0-0.5],[X0-0.5,NY-0.5,Z0-0.5],[X0-0.5,NY-0.5,NZ-0.5],[NX-0.5,NY-0.5,NZ-0.5]],
             [[X0-0.5,Y0-0.5,Z0-0.5],[X0-0.5,NY-0.5,Z0-0.5],[X0-0.5,NY-0.5,NZ-0.5],[X0-0.5,Y0-0.5,NZ-0.5]],
             [[X0-0.5,Y0-0.5,NZ-0.5],[NX-0.5,Y0-0.5,NZ-0.5],[NX-0.5,NY-0.5,NZ-0.5],[X0-0.5,NY-0.5,NZ-0.5]]
            ]


    collection = Poly3DCollection(sides, alpha=.05)
    collection.set_facecolor("blue")
    ax.add_collection3d(collection)
    
    
    #draw microtubuli
    if(microtub):
        var_tmp = 20
        
        ax.plot3D((46-var_tmp,46+var_tmp),(49,49),(49,49), color="orange", alpha=1)
        ax.plot3D((46,46),(49-var_tmp,49+var_tmp),(49,49), color="orange", alpha=1)
        ax.plot3D((46,46),(49,49),(49-var_tmp,49+var_tmp), color="orange", alpha=1)
        
        ax.plot3D((251-var_tmp,251+var_tmp),(49,49),(49,49), color="orange", alpha=1)
        ax.plot3D((251,251),(49-var_tmp,49+var_tmp),(49,49), color="orange", alpha=1)
        ax.plot3D((251,251),(49,49),(49-var_tmp,49+var_tmp), color="orange", alpha=1)
        
        
        
        var_tmp = 12
        
        ax.plot3D((46-var_tmp,46+var_tmp),(49-var_tmp,49+var_tmp),(49-var_tmp,49+var_tmp), color="orange", alpha=1)
        ax.plot3D((46-var_tmp,46+var_tmp),(49-var_tmp,49+var_tmp),(49+var_tmp,49-var_tmp), color="orange", alpha=1)
        ax.plot3D((46-var_tmp,46+var_tmp),(49+var_tmp,49-var_tmp),(49-var_tmp,49+var_tmp), color="orange", alpha=1)
        ax.plot3D((46-var_tmp,46+var_tmp),(49+var_tmp,49-var_tmp),(49+var_tmp,49-var_tmp), color="orange", alpha=1)
       
        ax.plot3D((251-var_tmp,251+var_tmp),(49-var_tmp,49+var_tmp),(49-var_tmp,49+var_tmp), color="orange", alpha=1)
        ax.plot3D((251-var_tmp,251+var_tmp),(49-var_tmp,49+var_tmp),(49+var_tmp,49-var_tmp), color="orange", alpha=1)
        ax.plot3D((251-var_tmp,251+var_tmp),(49+var_tmp,49-var_tmp),(49-var_tmp,49+var_tmp), color="orange", alpha=1)
        ax.plot3D((251-var_tmp,251+var_tmp),(49+var_tmp,49-var_tmp),(49+var_tmp,49-var_tmp), color="orange", alpha=1)
        
        
        
        var_tmp = 7
        
        ax.plot3D((251-var_tmp,251+var_tmp),(49-var_tmp,49+var_tmp),(49,49), color="orange", alpha=1)
        ax.plot3D((251-var_tmp,251+var_tmp),(49+var_tmp,49-var_tmp),(49,49), color="orange", alpha=1)
        ax.plot3D((251-var_tmp,251+var_tmp),(49,49),(49+var_tmp,49-var_tmp), color="orange", alpha=1)
        ax.plot3D((251-var_tmp,251+var_tmp),(49,49),(49-var_tmp,49+var_tmp), color="orange", alpha=1)
        ax.plot3D((251,251),(49-var_tmp,49+var_tmp),(49-var_tmp,49+var_tmp), color="orange", alpha=1)
        ax.plot3D((251,251),(49-var_tmp,49+var_tmp),(49+var_tmp,49-var_tmp), color="orange", alpha=1)
        
        ax.plot3D((46-var_tmp,46+var_tmp),(49-var_tmp,49+var_tmp),(49,49), color="orange", alpha=1)
        ax.plot3D((46-var_tmp,46+var_tmp),(49+var_tmp,49-var_tmp),(49,49), color="orange", alpha=1)
        ax.plot3D((46-var_tmp,46+var_tmp),(49,49),(49+var_tmp,49-var_tmp), color="orange", alpha=1)
        ax.plot3D((46-var_tmp,46+var_tmp),(49,49),(49-var_tmp,49+var_tmp), color="orange", alpha=1)
        ax.plot3D((46,46),(49-var_tmp,49+var_tmp),(49-var_tmp,49+var_tmp), color="orange", alpha=1)
        ax.plot3D((46,46),(49-var_tmp,49+var_tmp),(49+var_tmp,49-var_tmp), color="orange", alpha=1)


        #add rest
        collection = Poly3DCollection(sides, alpha=.05)
        collection.set_facecolor("blue")
        ax.add_collection3d(collection)
        
        
        
    
    if(fixed_aggresomes):
        
        #Aggresome coordinates
        aggresome_size=30
        xggL1, yggL1, zggL1 = NX/4-aggresome_size/2 , NY/4-aggresome_size/2 , NZ/4-aggresome_size/2
        xggU1, yggU1, zggU1 = NX/4+aggresome_size/2-1 , NY/4+aggresome_size/2-1 , NZ/4+aggresome_size/2-1

        xggL2, yggL2, zggL2 = 3*NX/4-aggresome_size/2 , NY/4-aggresome_size/2 , NZ/4-aggresome_size/2
        xggU2, yggU2, zggU2 = 3*NX/4+aggresome_size/2-1 , NY/4+aggresome_size/2-1 , NZ/4+aggresome_size/2-1
        
        
        points_agg1 = np.array(list(product([xggL1-0.5,xggU1+0.5],[yggL1-0.5,yggU1+0.5],[zggL1-0.5,zggU1+0.5])))
        points_agg2 = np.array(list(product([xggL2-0.5,xggU2+0.5],[yggL2-0.5,yggU2+0.5],[zggL2-0.5,zggU2+0.5])))
            
        for s,e in combinations(points_agg1, 2):
            if np.sum(np.abs(s-e)) in {xggU1-xggL1+1,yggU1-yggL1+1,zggU1-zggL1+1}:
                ax.plot3D(*zip(s, e), color="black", alpha=0.25)
            
        for s,e in combinations(points_agg2, 2):
            if np.sum(np.abs(s-e)) in {xggU2-xggL2+1,yggU2-yggL2+1,zggU2-zggL2+1}:
                ax.plot3D(*zip(s, e), color="black", alpha=0.25)
        
        
        sides_agg1 = [[[xggL1-0.5,yggL1-0.5,zggL1-0.5],[xggU1+0.5,yggL1-0.5,zggL1-0.5],[xggU1+0.5,yggU1+0.5,zggL1-0.5],[xggL1-0.5,yggU1+0.5,zggL1-0.5]],
                [[xggL1-0.5,yggL1-0.5,zggL1-0.5],[xggU1+0.5,yggL1-0.5,zggL1-0.5],[xggU1+0.5,yggL1-0.5,zggU1+0.5],[xggL1-0.5,yggL1-0.5,zggU1+0.5]],
                 [[xggU1+0.5,yggL1-0.5,zggL1-0.5],[xggU1+0.5,yggU1+0.5,zggL1-0.5],[xggU1+0.5,yggU1+0.5,zggU1+0.5],[xggU1+0.5,yggL1-0.5,zggU1+0.5]],
                 [[xggU1+0.5,yggU1+0.5,zggL1-0.5],[xggL1-0.5,yggU1+0.5,zggL1-0.5],[xggL1-0.5,yggU1+0.5,zggU1+0.5],[xggU1+0.5,yggU1+0.5,zggU1+0.5]],
                 [[xggL1-0.5,yggL1-0.5,zggL1-0.5],[xggL1-0.5,yggU1+0.5,zggL1-0.5],[xggL1-0.5,yggU1+0.5,zggU1+0.5],[xggL1-0.5,yggL1-0.5,zggU1+0.5]],
                 [[xggL1-0.5,yggL1-0.5,zggU1+0.5],[xggU1+0.5,yggL1-0.5,zggU1+0.5],[xggU1+0.5,yggU1+0.5,zggU1+0.5],[xggL1-0.5,yggU1+0.5,zggU1+0.5]]
                ]

        sides_agg2 = [[[xggL2-0.5,yggL2-0.5,zggL2-0.5],[xggU2+0.5,yggL2-0.5,zggL2-0.5],[xggU2+0.5,yggU2+0.5,zggL2-0.5],[xggL2-0.5,yggU2+0.5,zggL2-0.5]],
                [[xggL2-0.5,yggL2-0.5,zggL2-0.5],[xggU2+0.5,yggL2-0.5,zggL2-0.5],[xggU2+0.5,yggL2-0.5,zggU2+0.5],[xggL2-0.5,yggL2-0.5,zggU2+0.5]],
                 [[xggU2+0.5,yggL2-0.5,zggL2-0.5],[xggU2+0.5,yggU2+0.5,zggL2-0.5],[xggU2+0.5,yggU2+0.5,zggU2+0.5],[xggU2+0.5,yggL2-0.5,zggU2+0.5]],
                 [[xggU2+0.5,yggU2+0.5,zggL2-0.5],[xggL2-0.5,yggU2+0.5,zggL2-0.5],[xggL2-0.5,yggU2+0.5,zggU2+0.5],[xggU2+0.5,yggU2+0.5,zggU2+0.5]],
                 [[xggL2-0.5,yggL2-0.5,zggL2-0.5],[xggL2-0.5,yggU2+0.5,zggL2-0.5],[xggL2-0.5,yggU2+0.5,zggU2+0.5],[xggL2-0.5,yggL2-0.5,zggU2+0.5]],
                 [[xggL2-0.5,yggL2-0.5,zggU2+0.5],[xggU2+0.5,yggL2-0.5,zggU2+0.5],[xggU2+0.5,yggU2+0.5,zggU2+0.5],[xggL2-0.5,yggU2+0.5,zggU2+0.5]]
                ]
        
        collection_agg1 = Poly3DCollection(sides_agg1, alpha=.05)
        collection_agg2 = Poly3DCollection(sides_agg2, alpha=.05)
        
        
        collection_agg1.set_facecolor("yellow")
        collection_agg2.set_facecolor("yellow")
        
        
        ax.add_collection3d(collection_agg1)
        ax.add_collection3d(collection_agg2)
        
        

        
        
#-------------------------------------------------------------------------------------------------------------        
# FUNCTION FOR UPDATING MARKERS while using the slider
# https://stackoverflow.com/questions/48474699/marker-size-alpha-scaling-with-window-size-zoom-in-plot-scatter
# Thomas Kühn
#-------------------------------------------------------------------------------------------------------------

class MarkerUpdater:
    def __init__(self):
        ##for storing information about Figures and Axes
        self.figs = {}

        ##for storing timers
        self.timer_dict = {}

    def add_ax(self, ax, features=[]):
        update_now = True
        ax_dict = self.figs.setdefault(ax.figure,dict())
        ax_dict[ax] = {
            'xlim' : ax.get_xlim(),
            'ylim' : ax.get_ylim(),
            'zlim' : ax.get_zlim(),
            'figw' : ax.figure.get_figwidth(),
            'figh' : ax.figure.get_figheight(),
            'scale_s' : 1.0,
            'scale_a' : 1.0,
            'features' : [features] if isinstance(features,str) else features,
        }
        ax.figure.canvas.mpl_connect('draw_event', self.update_axes)

    def update_axes(self, event):

        for fig,axes in self.figs.items():
            if fig is event.canvas.figure:

                for ax, args in axes.items():
                    ##make sure the figure is re-drawn
                    #update = True

                    fw = fig.get_figwidth()
                    fh = fig.get_figheight()
                    fac1 = min(fw/args['figw'], fh/args['figh'])


                    xl = ax.get_xlim()
                    yl = ax.get_ylim()
                    zl = ax.get_zlim()
                    fac2 = min(
                        abs(args['xlim'][1]-args['xlim'][0])/abs(xl[1]-xl[0]),
                        abs(args['ylim'][1]-args['ylim'][0])/abs(yl[1]-yl[0]),
                        abs(args['zlim'][1]-args['zlim'][0])/abs(zl[1]-zl[0]),
                    )

                    ##factor for marker size
                    facS = (fac1*fac2)/args['scale_s']

                    ##factor for alpha -- limited to values smaller 1.0
                    facA = min(1.0,fac1*fac2)/args['scale_a']

                    ##updating the artists
                    if facS != 1.0:
                        for line in ax.lines:
                            if 'size' in args['features']:
                                line.set_markersize(line.get_markersize()*facS)

                            if 'alpha' in args['features']:
                                alpha = line.get_alpha()
                                if alpha is not None:
                                    line.set_alpha(alpha*facA)


                        for path in ax.collections:
                            if 'size' in args['features']:
                                path.set_sizes([s*facS**2 for s in path.get_sizes()])

                            if 'alpha' in args['features']:
                                alpha = path.get_alpha()
                                if alpha is not None:
                                    path.set_alpha(alpha*facA)

                        args['scale_s'] *= facS
                        args['scale_a'] *= facA

                self._redraw_later(fig)



    def _redraw_later(self, fig):
        if not plt.fignum_exists(fig):
            return 0
        timer = fig.canvas.new_timer(interval=10)
        timer.single_shot = True
        timer.add_callback(lambda : fig.canvas.draw_idle())
        timer.start()

        ##stopping previous timer
        if fig in self.timer_dict:
            self.timer_dict[fig].stop()

        ##storing a reference to prevent garbage collection
        self.timer_dict[fig] = timer

        
        
        
        
#-------------------------------------------------------------------------------------------------------------   
# FINDS CLUSTERS IN CELL & CALCULATES SIZES
# ADAPTED FROM: Charley Schaefer, University of York, 2020
#-------------------------------------------------------------------------------------------------------------
def get_cluster_list(matrix):
    M=len(matrix)
    N=len(matrix[0])
    G=len(matrix[0][0])
    analysed=np.zeros([M,N,G], int) # keep track of what matrix entries are analysed
    clusters=np.empty((0,4), int)
    Nclusters=0
    for i in range(M):
        for j in range(N):
            for p in range(G):
                if ((analysed[i,j,p]==0) and (matrix[i,j,p]>0)):  # New cluster
                    analysed[i,j,p]=1
                    Nclusters+=1
                    clustersize=1;
                    nn_list=np.empty((0,3), int);
                    Nnn=0 # length of nn_list

                    # get neighbour sites
                    if (p>0):
                        k=i; l=j; b=p-1;
                        if ( analysed[k,l,b]==0):
                            nn_list=np.append(nn_list, [[k,l,b]], axis=0); Nnn=Nnn+1;
                    if ( p<G-1 ):
                        k=i; l=j; b=p+1;
                        if ( analysed[k,l,b]==0):
                            nn_list=np.append(nn_list, [[k,l,b]], axis=0); Nnn=Nnn+1;
                            
                    
                    if (j>0):
                        k=i; l=j-1; b=p;
                        if ( analysed[k,l,b]==0):
                            nn_list=np.append(nn_list, [[k,l,b]], axis=0); Nnn=Nnn+1;
                    if ( j<N-1 ):
                        k=i; l=j+1; b=p;
                        if ( analysed[k,l,b]==0):
                            nn_list=np.append(nn_list, [[k,l,b]], axis=0); Nnn=Nnn+1;
                            
                            
                    if ( i>0 ):
                        k=i-1; l=j; b=p
                        if ( analysed[k,l,b]==0):
                            nn_list=np.append(nn_list, [[k,l,b]], axis=0); Nnn=Nnn+1;
                    if ( i<M-1 ):
                        k=i+1; l=j; b=p;
                        if (analysed[k,l,b]==0):
                            nn_list=np.append(nn_list, [[k,l,b]], axis=0); Nnn=Nnn+1;
                            

                            
                            
                    # analyse neighbours
                    while (Nnn>0): # analyse neighbours
                        m=nn_list[Nnn-1,0]; # analyse neighbour
                        n=nn_list[Nnn-1,1]; # analyse neighbour
                        g=nn_list[Nnn-1,2]; # analyse neighbour
                        
                        nn_list=np.delete( nn_list, Nnn-1, axis=0) ; Nnn=Nnn-1;
                        if ( analysed[m,n,g]!=1):
                            analysed[m,n,g]=1;
                            if ( matrix[m,n,g]>0): # new element in cluster
                                clustersize=clustersize+1;
                                
                                # Get new neighbours
                                # get neighbour sites
                                if ( g>0 ):
                                    k=m; l=n; b=g-1;
                                    if ( analysed[k,l,b]==0):
                                        nn_list=np.append(nn_list, [[k,l,b]], axis=0); Nnn=Nnn+1;
                                        
                                if ( g<G-1 ):
                                    k=m; l=n; b=g+1;
                                    if (analysed[k,l,b]==0):
                                        nn_list=np.append(nn_list, [[k,l,b]], axis=0); Nnn=Nnn+1;
                                        
                                        
                                if ( n>0 ):
                                    k=m; l=n-1; b=g;
                                    if ( analysed[k,l,b]==0):
                                        nn_list=np.append(nn_list, [[k,l,b]], axis=0); Nnn=Nnn+1;
                                        
                                if ( n<N-1 ):
                                    k=m; l=n+1; b=g;
                                    if ( analysed[k,l,b]==0):
                                        nn_list=np.append(nn_list, [[k,l,b]], axis=0); Nnn=Nnn+1;
                                        
                                        
                                if ( m>0 ):
                                    k=m-1; l=n; b=g;
                                    if ( analysed[k,l,b]==0):
                                        nn_list=np.append(nn_list, [[k,l,b]], axis=0); Nnn=Nnn+1;
                                        
                                if ( m<M-1 ):
                                    k=m+1; l=n; b=g;
                                    if (analysed[k,l,b]==0):
                                        nn_list=np.append(nn_list, [[k,l,b]], axis=0); Nnn=Nnn+1;
                                        
                                        
                                        
                                    
                    # end while loop: cluster completed
                    clusters=np.append(clusters, [[i,j,p,clustersize]], axis=0)
    return clusters