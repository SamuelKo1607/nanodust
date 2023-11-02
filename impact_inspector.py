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


class MyGUI:

    def __init__(self):
        self.myfont = ("Arial",12)
        self.root = tk.Tk()
        self.root.geometry("1200x750")
        self.root.title = "Waveform viewer"
        self.myframe = tk.Frame(self.root)
        self.myframe.columnconfigure(0, weight = 2)
        self.myframe.columnconfigure(1, weight = 2)
        self.myframe.columnconfigure(2, weight = 1)

        #daybox
        self.daybox = tk.Text(self.myframe,
                              font = self.myfont,
                              height=40)
        self.daybox.grid(row=0,column=0,rowspan=2,sticky="we",padx=10,pady=10)
        self.daybox.bind("<Return>",self.daybox_click)

        #eventbox
        self.eventbox = tk.Text(self.myframe,
                                font = self.myfont,
                                height=40)
        self.eventbox.grid(row=0,column=1,rowspan=2,sticky="we",padx=10,pady=10)
        self.eventbox.bind("<Return>",self.eventbox_click)

        #figures
        self.fig, self.ax = construct_subplots(dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.myframe)
        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0,column=2,padx=10,pady=10)

        #infobox
        self.infobox = tk.Text(self.myframe,
                                font = self.myfont,
                                height=10)
        self.infobox.grid(row=1,column=2,sticky="we",padx=10,pady=10)

        self.myframe.pack(fill="x")

        #list days
        days = load_all_days()
        YYYYMMDDs = [day.YYYYMMDD for day in days]
        for YYYYMMDD in YYYYMMDDs:
            self.daybox.insert(tk.END,add_dash(YYYYMMDD)+"\n")
        self.daybox.config(state=tk.DISABLED)



        self.root.mainloop()


    def list_events(self,
                    YYYYMMDD,
                    impacts_location=os.path.join("998_generated","impacts","")):
        """
        A procedure to list all the events on the selected day 
        into the self.eventbox.

        Parameters
        ----------
        YYYYMMDD : str
            The date of interest.

        impacts_location : str
            The location of the impact pickles.

        Returns
        -------
        None.

        """
        self.eventbox.config(state=tk.NORMAL)
        self.eventbox.delete('1.0', tk.END)
        files = glob.glob(impacts_location+"*.pkl")
        file = [f for f in files if YYYYMMDD in f][0]
        self.impacts = load_list(file,"")
        impact_indices = [impact.index for impact in self.impacts]
        impact_indices.sort()
        self.eventbox.insert(tk.END,YYYYMMDD+"\n")
        self.eventbox.insert(tk.END,"---------------"+"\n")
        for i,index in enumerate(impact_indices):
            polarity = [impact.polarity for impact in self.impacts if impact.index == index][0]
            if np.sign(polarity)>0:
                sign = "+"
            elif np.sign(polarity)<0:
                sign = "-"
            else:
                sign = "0"
            antenna = [impact.antenna_hit for impact in self.impacts if impact.index == index][0]
            ant = antenna*"*"
            self.eventbox.insert(tk.END,f"{index} : {sign} {ant}\n")
        self.eventbox.config(state=tk.DISABLED)


    def info_impact(self,impact):
        """
        A procedure to print impact info in the self.infobox.

        Parameters
        ----------
        impact : Impact
            The impact of interest.

        Returns
        -------
        None.

        """
        self.infobox.config(state=tk.NORMAL)
        self.infobox.delete('1.0', tk.END)
        self.infobox.insert(tk.END,f"Impact time: {impact.datetime}\n")
        self.infobox.insert(tk.END,f"Sampling rate: {impact.sampling_rate}\n")
        self.infobox.insert(tk.END,f"Amplitude: {impact.amplitude:.4f} V\n")
        self.infobox.insert(tk.END,f"Symmetry: {impact.symmetry:.3f}\n")
        self.infobox.insert(tk.END,f"Polarity: {int(impact.polarity)}\n")
        self.infobox.insert(tk.END,f"Antenna hit: {int(impact.antenna_hit)}\n")
        self.infobox.config(state=tk.DISABLED)


    def plot_impact(self,
                    impact,
                    xrange=(-900,1600),
                    colors=["red","blue"],
                    labels=["filtered","orig"]):
        """
        Plot the impact in self.canvas. 

        Parameters
        ----------
        impact : Impact
            An Impact object containing all the relevant data.
        xrange : tuple, optional
            The definition of the zoom in the LHS panels, 
            microseconds since impact. The default is [-900,1600].
        colors : list of str, optional
            The colors for the waveforms. The default is ["blue","red"].
        labels : list of str, optional
            The list of labels for the waveforms. 
            The default is ["filtered","orig"].

        Returns
        -------
        None.

        """
        self.fig, self.ax = construct_subplots(dpi=100)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.myframe)

        ylim = 1.1*max([max(np.abs(impact.wf_filtered[0])),
                        max(np.abs(impact.wf_filtered[1])),
                        max(np.abs(impact.wf_filtered[2]))])
        xmin = xrange[0]
        xmax = xrange[1]

        #time windows is centered about the most important peak
        epoch_center = impact.epoch[impact.extreme_index]
        epoch = impact.epoch

        for i,wf in enumerate([impact.wf_orig,impact.wf_filtered]):

            self.ax[0,0].plot((epoch-epoch_center)/1000,wf[0],color=colors[i],lw=0.5)
            self.ax[1,0].plot((epoch-epoch_center)/1000,wf[1],color=colors[i],lw=0.5)
            self.ax[2,0].plot((epoch-epoch_center)/1000,wf[2],color=colors[i],lw=0.5)
        
            self.ax[0,1].plot((epoch-epoch_center)/1000000,wf[0],color=colors[i],lw=0.5,
                         label = labels[i])
            self.ax[1,1].plot((epoch-epoch_center)/1000000,wf[1],color=colors[i],lw=0.5)
            self.ax[2,1].plot((epoch-epoch_center)/1000000,wf[2],color=colors[i],lw=0.5)

        self.ax[0,1].legend(fontsize="small")
        for axis in self.ax[:,0]:
            axis.vlines(0,-ylim,ylim,color="grey",zorder=0,alpha=0.2)
            axis.hlines(0,xmin,xmax,color="grey",zorder=0,alpha=0.2)
            axis.set_ylim(-ylim,ylim)
            axis.set_xlim(xmin,xmax)
        for axis in self.ax[:,1]:
            axis.vlines(0,-ylim,ylim,color="grey",zorder=0,alpha=0.2)
            axis.hlines(0,xmin/1000,xmax/1000,color="grey",zorder=0,alpha=0.2)
            axis.set_ylim(-ylim,ylim)
            axis.axvspan(xmin/1000,xmax/1000,alpha=0.3,color='lightgrey')

        self.canvas.draw()
        self.canvas.get_tk_widget().grid(row=0,column=2,padx=10,pady=10)


    def daybox_click(self, event):
        """
        The action called on <Return> while a date is selected in daybox.

        Parameters
        ----------
        event : tkinter.event
            The event that pointed to the action.

        Returns
        -------
        None.

        """
        YYYYMMDD_dashed = self.daybox.get(tk.SEL_FIRST,tk.SEL_LAST)
        self.list_events(remove_dash(YYYYMMDD_dashed))


    def eventbox_click(self, event):
        """
        The action called on <Return> while an index is selected in eventbox.

        Parameters
        ----------
        event : tkinter.event
            The event that pointed to the action.

        Returns
        -------
        None.

        """
        index = self.eventbox.get(tk.SEL_FIRST,tk.SEL_LAST)
        self.impact = [impact for impact in self.impacts if impact.index == int(index)][0]
        self.info_impact(self.impact)
        #self.fig.close()
        self.plot_impact(self.impact)



def add_dash(YYYYMMDD):
    """
    To make the date a bit more readable.

    Parameters
    ----------
    YYYYMMDD : str
        Date.

    Returns
    -------
    YYYYMMDD_dashed
        YYYY-MM-DD.
    """
    YYYYMMDD_dashed = YYYYMMDD[:4]+"-"+YYYYMMDD[4:6]+"-"+YYYYMMDD[6:]
    return YYYYMMDD_dashed


def remove_dash(YYYYMMDD_dashed):
    """
    To remove the dashes from a string.

    Parameters
    ----------
    YYYYMMDD_dashed : str
        YYYY-MM-DD.

    Returns
    -------
    YYYYMMDD : str
        Date.
    """
    YYYYMMDD = YYYYMMDD_dashed.replace("-","")
    return YYYYMMDD


def construct_subplots(dpi=100):
    """
    To prerare empy axes for plotting waveforms.

    Parameters
    ----------
    dpi : float, optional
        The resolution. The default is 100.

    Returns
    -------
    fig : matplotblib.figure
        An empty figure.
    ax : matplotblib.Figure.Axes
        Empty axes.
    """
    fig = plt.figure(figsize=(12,5),dpi=dpi)
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



if __name__ == "__main__":
    MyGUI()















