import numpy as np
import matplotlib.pyplot as plt

x = np.arange(0,20,0.01)
t = 12*x/(2*np.pi)

fig = plt.figure(figsize=(3,3))
gs = fig.add_gridspec(3, hspace=.05)
axs = gs.subplots()
axs[0].text(2.25,22,"a)",backgroundcolor="white",zorder=10)
axs[1].text(2.25,22,"b)",backgroundcolor="white",zorder=10)
axs[2].text(2.25,22,"c)",backgroundcolor="white",zorder=10)
for ax in axs:
    ax.fill_between(t,16*np.cos(x)**2+7,color=u"#FFFD82",ec=u"#2D3047",lw=0.5,label="Observed rate")
    ax.set_ylim(0,28)
    ax.set_xlim(0,37)
    ax.xaxis.set_ticks(np.arange(0, max(t), 12))
    ax.xaxis.set_label_coords(0.5, -0.15)
    ax.set_aspect(37/(28*2*1.618))

axs[1].set_ylabel("Rate [/day]")
axs[2].set_xlabel("Month")
axs[0].tick_params(axis='x',labelbottom=False)
axs[1].tick_params(axis='x',labelbottom=False)

axs[0].fill_between(t,4,color=u"#1B998B",ec="black",lw=0.3,label="Non-hyperb.")
axs[1].fill_between(t,4*np.cos(x)**2+2,color=u"#1B998B",ec=u"#2D3047",lw=0.5)
axs[2].fill_between(t,4*np.sin(x)**2+2,color=u"#1B998B",ec=u"#2D3047",lw=0.5)

axs[0].plot(t,16*np.cos(x)**2+7-4,                 lw=1.5,ls="dashed",color=u"#F46036",label="Hyperbolic")
axs[1].plot(t,16*np.cos(x)**2+7-(4*np.cos(x)**2+2),lw=1.5,ls="dashed",color=u"#F46036")
axs[2].plot(t,16*np.cos(x)**2+7-(4*np.sin(x)**2+2),lw=1.5,ls="dashed",color=u"#F46036")

#get handles and labels
handles, labels = axs[0].get_legend_handles_labels()
#specify order of items in legend
order = [1,2,0]
#add legend to plot
axs[0].legend([handles[idx] for idx in order],[labels[idx] for idx in order],edgecolor="none",framealpha=0.9,fontsize="x-small",loc=1)

fig.savefig("C:\\Users\\skoci\\Disk Google\\000 škola\\UIT\\presentations\\beta velocity\\background_profiles.pdf", format='pdf', dpi=600, bbox_inches="tight")
