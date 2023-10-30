import tkinter as tk
from matplotlib.backend_bases import key_press_handler
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
import glob
import os


from nano_load_days import Impact
from nano_load_days import Day
from nano_load_days import load_all_days
from nano_load_days import load_list

impacts_location = os.path.join("998_generated","impacts","")



class MyGUI:

    def __init__(self):
        self.myfont = ("Arial",12)
        self.root = tk.Tk()
        self.root.geometry("1000x700")
        self.root.title = "Waveform viewer"
        self.myframe = tk.Frame(self.root)
        self.myframe.columnconfigure(0, weight = 3)
        self.myframe.columnconfigure(1, weight = 3)
        self.myframe.columnconfigure(2, weight = 1)

        #daybox
        self.daybox = tk.Text(self.myframe,
                              font = self.myfont,
                              height=35)
        self.daybox.grid(row=0,column=0,rowspan=2,sticky="we",padx=20,pady=20)
        self.daybox.bind("<Return>",self.daybox_click)

        #eventbox
        self.eventbox = tk.Text(self.myframe,
                                font = self.myfont,
                                height=35)
        self.eventbox.grid(row=0,column=1,rowspan=2,sticky="we",padx=20,pady=20)
        self.eventbox.bind("<Return>",self.eventbox_click)

        #figures
        self.fig, ax = construct_subplots(dpi=100)
        self.canvas1 = FigureCanvasTkAgg(self.fig, master=self.myframe)
        self.canvas1.draw()
        self.canvas1.get_tk_widget().grid(row=0,column=2,padx=20,pady=20)

        #infobox
        self.infobox = tk.Text(self.myframe,
                                font = self.myfont,
                                height=9)
        self.infobox.grid(row=1,column=2,sticky="we",padx=20,pady=20)

        self.myframe.pack(fill="x")

        #list days
        days = load_all_days()
        YYYYMMDDs = [day.YYYYMMDD for day in days]
        for YYYYMMDD in YYYYMMDDs:
            self.daybox.insert(tk.END,YYYYMMDD+"\n")
        self.daybox.config(state=tk.DISABLED)



        self.root.mainloop()


    def list_events(self,YYYYMMDD):
        self.eventbox.config(state=tk.NORMAL)
        self.eventbox.delete('1.0', tk.END)
        files = glob.glob(impacts_location+"*.pkl")
        file = [f for f in files if YYYYMMDD in f][0]
        self.impacts = load_list(file,"")
        impact_indices = [impact.index for impact in self.impacts]
        impact_indices.sort()
        self.eventbox.insert(tk.END,YYYYMMDD+"\n")
        self.eventbox.insert(tk.END,"---------------"+"\n")
        for i in impact_indices:
            self.eventbox.insert(tk.END,str(i)+"\n")
        self.eventbox.config(state=tk.DISABLED)


    def info_impact(self,impact):
        self.infobox.config(state=tk.NORMAL)
        self.infobox.delete('1.0', tk.END)
        self.infobox.insert(tk.END,str(impact.datetime)+"\n")
        self.infobox.insert(tk.END,str(impact.sampling_rate)+"\n")
        self.infobox.insert(tk.END,str(impact.amplitude)+"\n")
        self.infobox.insert(tk.END,str(impact.symmetry)+"\n")
        self.infobox.insert(tk.END,str(impact.polarity)+"\n")
        self.infobox.config(state=tk.DISABLED)



    def daybox_click(self, event):
        YYYYMMDD = self.daybox.get(tk.SEL_FIRST,tk.SEL_LAST)
        print(YYYYMMDD)
        self.list_events(YYYYMMDD)


    def eventbox_click(self, event):
        index = self.eventbox.get(tk.SEL_FIRST,tk.SEL_LAST)
        self.impact = [impact for impact in self.impacts if impact.index == int(index)][0]
        self.info_impact(self.impact)





def construct_subplots(dpi=100):
    fig = plt.figure(figsize=(7,4),dpi=dpi)
    gs = fig.add_gridspec(3,2,width_ratios=[1, 2],hspace=0.1,wspace=0.1)
    ax = gs.subplots()
    ax[0,0].xaxis.set_ticklabels([])
    ax[1,0].xaxis.set_ticklabels([])
    ax[0,1].xaxis.set_ticklabels([])
    ax[1,1].xaxis.set_ticklabels([])
    ax[0,1].yaxis.set_ticklabels([])
    ax[1,1].yaxis.set_ticklabels([])
    ax[2,1].yaxis.set_ticklabels([])
    ax[0,0].set_ylabel("V1 [V]")
    ax[1,0].set_ylabel("V2 [V]")
    ax[2,0].set_ylabel("V3 [V]")
    ax[2,0].set_xlabel("time [$\mu s$]")
    ax[2,1].set_xlabel("time [$ms$]")
    return fig,ax




MyGUI()















