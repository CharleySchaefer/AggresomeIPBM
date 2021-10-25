# Preparations
#-------------------------------------------------------------------------

import numpy as np
import matplotlib as mpl
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
import pickle
import tkinter as tk
import time
import tkinter.filedialog as dialog
import tkinter.messagebox as msg_box

import mylib

#-------------------------------------------------------------------------

def main():
    
    # Set matplotlib backend
    mpl.use("Qt5Agg")
    
    # Default Cell coordinates & Sizes
    NX, NY, NZ = 300, 100, 100
    NYZ = NY*NZ
    X0 ,Y0, Z0 = 0,0,0
    
    # Default Nucleus coordinates
    xnuL, ynuL, znuL = 93,11,11
    xnuU, ynuU, znuU = 205,87,87
    
    
    
    
    
    # GUI for loading the simulation data file and for starting plot
    class gui_sim:
        def __init__(self):
            
            self.data = 'nan'
            self.filename=""
            
            
            window = tk.Tk()
            window.title("Visualisation Manager")
            
            
            # Checkbox preparations for single phase, microtubule & cubic aggresomes
            self.single_phase = tk.IntVar()
            self.microtubule = tk.IntVar()
            self.fixed_aggresomes = tk.IntVar()
            
            
            load_label = tk.Label(window, text="Load Simulation Data File")
            load_button = tk.Button(window, text="Load...", command=self.load_data, height=2, width=6)
            self.load_label_2 = tk.Label(window, text="No file chosen.")
            
            vis_label = tk.Label(window, text="Open Visualisation Window")
            vis_button = tk.Button(window, text="Vis", command=self.visualize, height=2, width=6)
            
            
            load_label.grid(row=0, column=0)
            load_button.grid(row=0, column=1)
            self.load_label_2.grid(row=1,column=0)
            vis_label.grid(row=2, column=0)
            vis_button.grid(row=2, column=1)
            
            single_phase_checkbox = tk.Checkbutton(window, text="Single Phase", variable=self.single_phase).grid(row=3, column=0)
            microtubule_checkbox = tk.Checkbutton(window, text="Microtubule", variable=self.microtubule).grid(row=3, column=1)
            fixed_aggresomes_checkbox = tk.Checkbutton(window, text="Fixed Aggresomes", variable=self.fixed_aggresomes).grid(row=4, column=0)
            
            
            window.mainloop()
    
    
        def load_data(self):
            """Outpacks simulation data

            Opens a dialog window in that you can select the simulation output file.
            It then returns a pandas data file with the simulation data.
            """

            filename_tmp=dialog.askopenfilename()
            if(filename_tmp!=""):
                self.filename = filename_tmp
                self.load_label_2['text'] = "Wait, loading!"
                self.data = pd.read_csv(self.filename, sep="\t", header=4, dtype={0: np.float64, 1: np.int32, 2:np.int32,3:np.int32,4:np.int32})
                self.load_label_2['text'] = "Finished loading."
                
                
                
        def visualize(self):
            if(self.filename==""):
                msg_box.showerror("Error", "Please load a file first.")
                return
            
                
            plt.close('all')
            particle_num = int((len(self.data.iloc[0])-1)/3)

            #Initialising Plotting
            self.fig = plt.figure("VIS2_normal")
            self.ax = self.fig.add_subplot(111,projection='3d')
            self.ax.set_aspect('equal')

            last_ts = 0
            cur_time = 0

            self.sct, = self.ax.plot([], [], [], "o", markersize=1, color="red", alpha=0.3)

            #Draw the Box
            mylib.draw_box(self.ax, self.single_phase.get(), self.microtubule.get(), self.fixed_aggresomes.get(), X0, NX, Y0, NY, Z0, NZ, xnuL, xnuU, ynuL, ynuU, znuL, znuU)

            len_timesteps = len(self.data["Time"])

            #Slider
            self.ax_time = plt.axes([0.25, 0.05, 0.65, 0.02], facecolor="#ECECEC")
            self.slide_time = Slider(self.ax_time, 'Timestep', 0, len_timesteps-1, valinit=0,valfmt="%d")
            plt.xlabel("Time: 0.0")

            cur_ts = int(self.slide_time.val)
            cur_time = self.data["Time"][cur_ts]
            self.sct.set_data(np.array(self.data.iloc[0][1:particle_num+1]), np.array(self.data.iloc[0][particle_num+1:2*particle_num+1]))
            self.sct.set_3d_properties(np.array(self.data.iloc[0][2*particle_num+1:3*particle_num+1]))


            def update(val):
                cur_ts = int(self.slide_time.val)
                cur_time = self.data["Time"][cur_ts]

                self.sct.set_data(np.array(self.data.iloc[cur_ts][1:particle_num+1]), np.array(self.data.iloc[cur_ts][particle_num+1:2*particle_num+1]))
                self.sct.set_3d_properties(np.array(self.data.iloc[cur_ts][2*particle_num+1:3*particle_num+1]))
                plt.xlabel("Time: "+"{0:.2f}".format(cur_time))

            my_updater = mylib.MarkerUpdater()

            #setting up the updater
            my_updater.add_ax(self.ax, ['size', 'alpha']) ##scatter plot, marker size and alpha    

            #updating if slider value changes
            self.slide_time.on_changed(update)

            #Using equal axis function from above
            mylib.set_axes_equal(self.ax)
            
                        
            
            
    
    window_instance = gui_sim()
    
    
    
if __name__ == "__main__":
    main()
