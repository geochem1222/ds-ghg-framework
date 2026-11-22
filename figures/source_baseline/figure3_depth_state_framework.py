"""Editable conceptual figure; all trajectories are prescribed hypotheses."""
from pathlib import Path
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'figure3_depth_state_framework'
WIDTH, HEIGHT = 183, 154
INK, BLUE, GREEN, GREY = '#263238', '#257C9E', '#3B7E68', '#777F84'
LIGHT, HAIR = '#F0F4F5', '#CCD3D6'
mpl.rcParams.update({'font.family':'sans-serif', 'font.sans-serif':['Arial','DejaVu Sans'],
    'font.size':7, 'svg.fonttype':'none', 'pdf.fonttype':42, 'text.color':INK})
fig = plt.figure(figsize=(WIDTH/25.4,HEIGHT/25.4), facecolor='white')
ax = fig.add_axes([0,0,1,1])
ax.set(xlim=(0,WIDTH),ylim=(0,HEIGHT),aspect='equal')
ax.axis('off')
labels=[]
def label(x,y,s,size=7,color=INK,ha='left',bold=False,**kw):
    obj=ax.text(x,y,s,fontsize=size,color=color,ha=ha,va='center',
        weight='bold' if bold else 'normal',linespacing=1.2,**kw)
    labels.append(obj)
    return obj
def arrow(a,b,color=GREY,lw=.85,style='-|>'):
    ax.add_patch(FancyArrowPatch(a,b,arrowstyle=style,mutation_scale=7,
        lw=lw,color=color,shrinkA=1,shrinkB=1))

label(3,149,'a',9,bold=True)
label(9,149,'Depth history changes the conditions of the response',8.4,bold=True)
label(5,139,'External forcing $U(t)$',7.5,BLUE)
label(5,133.5,'Climate, light and resource inputs',6.8)
label(5,123,r'Depth history $d(t)$, $\dot d(t)$',7.5,BLUE)
label(5,116,'Atmospheric pressure',6.8,GREY)
arrow((48,135),(60,135),BLUE)
arrow((48,124),(60,124),BLUE)
ax.add_patch(Rectangle((61,115),66,26,facecolor=LIGHT,edgecolor='none'))
label(94,136.8,'Temperature, light, oxidants, resources',7.1,ha='center')
label(94,131.4,'Gas stores, mixing and legacy',7.1,ha='center')
label(72,121.5,'Physics\n'+r'$P_{\mathrm{phys}}$',6.8,BLUE,ha='center')
label(94,121.5,'Reactions\n'+r'$R_{\mathrm{rxn}}$',6.8,GREEN,ha='center')
label(116,121.5,'Transport\n'+r'$J_{\mathrm{trans}}$',6.8,GREY,ha='center')
arrow((128,133),(139,133),BLUE)
label(159,136,'Local flux $f_g$ (3)',7.1,ha='center')
arrow((158,132),(158,126),BLUE,lw=.7)
label(159,122.5,'Area-integrated $F_g$ (1)',7.1,ha='center')
label(159,116.5,'CO$_2$  CH$_4$  N$_2$O',7.5,ha='center',bold=True)
label(94,111.8,'State budget (2)',6.7,GREY,ha='center')
ax.plot([48,54,54,135,135],[123,123,107.2,107.2,130],color=GREY,lw=.65)
arrow((48,116),(54,116),GREY,lw=.65)
arrow((135,130),(140,130),GREY,lw=.65)
label(91,107.2,'Pressure-sensitive release',6.5,GREY,ha='center',
      bbox=dict(facecolor='white',edgecolor='none',pad=1))
ax.plot([4,179],[103,103],color=HAIR,lw=.6)

ax.plot([6,13],[98.5,98.5],color=GREY,lw=1.2,ls=(0,(3,2)))
label(15,98.5,'Shallower / stable depth',6.5,GREY)
ax.plot([101,108],[98.5,98.5],color=BLUE,lw=1.2)
label(110,98.5,'Deeper / depth cycle',6.5,BLUE)
xs=[8,68,128]
titles=['Thermal buffering','Oxidant limitation','Retention and release']
subtitles=['Oxidation fraction maintained','Oxidation loss offsets buffering','Excess net CH$_4$ supply matched']
for i,(x,title,sub) in enumerate(zip(xs,titles,subtitles)):
    label(x-4,91.5,'bcd'[i],9,bold=True)
    label(x+2,91.5,title,7.3,bold=True)
    label(x,85.7,sub,6.2,GREY)
for x in (62,122):
    ax.plot([x,x],[8,93],color=HAIR,lw=.5)

# Smooth shapes are specified solely to draw conceptual hypotheses. They do not
# estimate environmental parameters, and are not produced by the DS-GHG model.
t=np.linspace(0,1,2001)
dt=t[1]-t[0]
def pulse(start,end,height=1):
    v=np.zeros_like(t)
    m=(t>=start)&(t<=end)
    v[m]=height*np.sin(np.pi*(t[m]-start)/(end-start))**2
    return v
def integrate(v):
    return np.r_[0,np.cumsum((v[1:]+v[:-1])*.5*dt)]
heat=pulse(.13,.32)
temp_shallow=pulse(.15,.43)
temp_deep=pulse(.20,.62,.50)
flux_ref=pulse(.17,.48)
flux_thermal=pulse(.22,.64,.48)
flux_redox=pulse(.22,.83,1.22)
flux_delayed=pulse(.50,.92)
flux_delayed*=np.trapezoid(flux_ref,t)/np.trapezoid(flux_delayed,t)
depth=np.interp(t,[0,.10,.28,.47,.67,1],[0,0,1,1,0,0])
inventory=integrate(flux_ref-flux_delayed)
assert inventory.min()>-1e-10
assert abs(inventory[-1])<1e-10
assert np.trapezoid(flux_thermal,t)<np.trapezoid(flux_ref,t)<np.trapezoid(flux_redox,t)
ref_area=integrate(flux_ref)
cum_scale=integrate(flux_redox)[-1]

def plot_row(x,y,h,v1,v2,scale,shade=False):
    w=48
    ax.plot([x,x+w],[y,y],color=HAIR,lw=.6)
    if shade:
        ax.add_patch(Rectangle((x+w*.13,y),w*.19,h,facecolor='#E7EAEB',
            edgecolor='none',alpha=.65,zorder=0))
    ax.plot(x+w*t,y+h*v1/scale,color=GREY,lw=1.2,ls=(0,(3,2)))
    ax.plot(x+w*t,y+h*v2/scale,color=BLUE,lw=1.35)

for j,x in enumerate(xs):
    label(x,79.8,'Bed-temperature anomaly' if j<2 else 'Water-depth change',6.7)
    if j<2:
        plot_row(x,66,10,temp_shallow,temp_deep,1,shade=True)
    else:
        plot_row(x,66,10,np.zeros_like(t),depth,1,shade=True)
    label(x,59.3,'CH$_4$ emission anomaly',6.7)
    flux=[flux_thermal,flux_redox,flux_delayed][j]
    plot_row(x,37,17,flux_ref,flux,1.25,shade=True)
    label(x,30.5,'Cumulative excess CH$_4$',6.7)
    plot_row(x,10,16,ref_area,integrate(flux),cum_scale)
    arrow((x,7),(x+48,7),GREY,lw=.6)
    label(x+24,3,'Time',6.4,GREY,ha='center')

fig.canvas.draw()
renderer=fig.canvas.get_renderer()
outside=[]
collisions=[]
for i,obj in enumerate(labels):
    bounds=obj.get_window_extent(renderer)
    if bounds.x0<0 or bounds.y0<0 or bounds.x1>fig.bbox.width or bounds.y1>fig.bbox.height:
        outside.append(obj.get_text())
    for other in labels[i+1:]:
        if bounds.overlaps(other.get_window_extent(renderer)):
            collisions.append([obj.get_text(),other.get_text()])
assert not outside,outside
assert not collisions,collisions
for ext in ('svg','pdf'):
    fig.savefig(str(OUT)+'.'+ext)
fig.savefig(str(OUT)+'.png',dpi=600)
fig.savefig(str(OUT)+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
fig.savefig(str(OUT)+'_preview.png',dpi=220)
(ROOT/'figure3_layout_checks.json').write_text(json.dumps(dict(
    width_mm=WIDTH,height_mm=HEIGHT,off_canvas_labels=outside,text_collisions=collisions,
    empirical_data=False,supplementary_model_output=False,curves='Prescribed conceptual hypotheses',
    thermal_emission_peak_smaller=bool(flux_thermal.max()<flux_ref.max()),
    redox_case_emission_larger=bool(flux_redox.max()>flux_ref.max()),
    storage_excess_inventory_min=float(inventory.min()),
    storage_final_excess_inventory=float(inventory[-1]),
    storage_cumulative_release_equal=bool(abs(inventory[-1])<1e-10),
    comparison='Perturbed minus unperturbed under each corresponding depth history',
    arbitrary_time_and_magnitude=True),indent=2),encoding='utf-8')
plt.close(fig)
