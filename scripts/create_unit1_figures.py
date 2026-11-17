# Optional source for the three original diagrams in Unit 1.
# Requires matplotlib and numpy; NOT required to render the book.
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
out=Path(__file__).resolve().parents[1]/'fig';out.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'savefig.facecolor':'white','figure.facecolor':'white'})
# Figure 1: data structures
fig, axs = plt.subplots(2,2,figsize=(10,5.7));fig.subplots_adjust(hspace=.5,wspace=.26,left=.06,right=.97,top=.86,bottom=.09)
fig.suptitle('Four ways to organize data in R',fontsize=17,fontweight='bold')
def cell(ax,x,y,w,h,text,color='#DCE8EC',fontsize=11):
 p=FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.025,rounding_size=0.06',facecolor=color,edgecolor='#52636D',linewidth=1.15)
 ax.add_patch(p);ax.text(x+w/2,y+h/2,text,va='center',ha='center',fontsize=fontsize)
for ax in axs.ravel():ax.set_xlim(0,5);ax.set_ylim(0,3);ax.axis('off')
ax=axs[0,0];ax.set_title('Vector',loc='left',fontweight='bold');
for i,s in enumerate(['2','4','6','8']):cell(ax,.28+i*1.14,1.1,1.02,.72,s,color='#CFE5DB')
ax.text(.28,.60,'One dimension; one atomic type',fontsize=9.4,color='#475569')
ax=axs[0,1];ax.set_title('Matrix',loc='left',fontweight='bold');
for row in range(2):
 for col in range(3):cell(ax,.65+col*1.16,1.87-row*.84,1.02,.72,str(row*3+col+1),color='#DCE4F7')
ax.text(.65,.42,'Rows × columns; one atomic type',fontsize=9.4,color='#475569')
ax=axs[1,0];ax.set_title('Data frame',loc='left',fontweight='bold')
for j,(hdr,items,colr) in enumerate([('species',['Bombus','Apis'],'#CFE5DB'),('count',['12','8'],'#DCE4F7'),('native',['TRUE','FALSE'],'#F3DFC6')]):
 x=.20+j*1.56
 cell(ax,x,2.15,1.44,.55,hdr,color=colr,fontsize=9.6)
 for i,item in enumerate(items):cell(ax,x,1.46-i*.69,1.44,.57,item,color='#F7F9FB',fontsize=9.4)
ax.text(.20,.10,'Columns can have different types; equal lengths',fontsize=9.0,color='#475569')
ax=axs[1,1];ax.set_title('List',loc='left',fontweight='bold')
cell(ax,.30,1.86,1.23,.69,'species',color='#CFE5DB')
cell(ax,1.64,1.86,2.85,.69,'c("Bombus", "Apis")',color='#F7F9FB',fontsize=10)
cell(ax,.30,1.01,1.23,.69,'counts',color='#DCE4F7')
cell(ax,1.64,1.01,2.85,.69,'matrix(1:6, 2, 3)',color='#F7F9FB',fontsize=10)
ax.text(.3,.33,'Components may differ in type and size',fontsize=9.4,color='#475569')
fig.savefig(out/'unit1-data-structures.png',dpi=180,bbox_inches='tight');plt.close(fig)
# Figure 2: logistic trajectory conceptual
fig,ax=plt.subplots(figsize=(9.5,4.4));K=100;n0=10;ts=np.arange(30)
for r,color,lab in [(0.6,'#206D79','r = 0.6: gradual approach'),(1.9,'#C77831','r = 1.9: overshoots, then oscillates')]:
 N=[float(n0)]
 for _ in range(29):N.append(N[-1]+r*N[-1]*(1-N[-1]/K))
 ax.plot(ts,N,label=lab,color=color,linewidth=2.3)
ax.axhline(K,color='#555C66',ls='--',linewidth=1.5,label='Carrying capacity, K = 100')
ax.set(xlabel='Time step',ylabel='Population size',xlim=(0,29),ylim=(0,140));ax.set_title('Density dependence in a discrete-time logistic model',loc='left',fontweight='bold',fontsize=13)
ax.legend(frameon=False,loc='lower right',fontsize=9.2);ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.14)
fig.savefig(out/'unit1-logistic-trajectories.png',dpi=180,bbox_inches='tight');plt.close(fig)
# Figure 3: conceptual R0 vs r difference
fig,axs=plt.subplots(1,2,figsize=(9.8,4.1),gridspec_kw={'width_ratios':[.94,1.06]});fig.subplots_adjust(wspace=.32,left=.05,right=.97,top=.81,bottom=.17)
fig.suptitle('Two different quantities in early epidemic growth',fontweight='bold',fontsize=15)
ax=axs[0];ax.set_xlim(0,5);ax.set_ylim(0,4);ax.axis('off');ax.set_title(r'$R_0$: average secondary infections',fontsize=11.2)
cell(ax,.23,1.72,1.18,.73,'1 case',color='#DCE4F7')
for y in [2.6,1.73,.85]:
 ax.add_patch(FancyArrowPatch((1.52,2.08),(2.65,y+.31),arrowstyle='->',mutation_scale=14,color='#657581',linewidth=1.4))
 cell(ax,2.75,y,1.52,.65,'new case',color='#F3DFC6',fontsize=10)
ax.text(.22,.2,'Illustrative: R₀ = 3',fontsize=9.6,color='#475569')
ax=axs[1];t=np.linspace(0,30,300);rr=np.log(3)/12;N=np.exp(rr*t)
ax.plot(t,N,color='#206D79',lw=2.5);ax.scatter([0],[1],color='#206D79',s=30,zorder=5)
ax.set_title(r'$r$: exponential growth rate per day',fontsize=11.2);ax.set_xlabel('Days');ax.set_ylabel('Early-stage cases (illustrative)');ax.grid(alpha=.15);ax.spines[['top','right']].set_visible(False)
fig.text(.50,.015,r'Using a fixed generation interval of 12 days: $r \approx \ln(R_0)/12$. This approximation does not model epidemic saturation.',ha='center',fontsize=9.5,color='#475569')
fig.savefig(out/'unit1-reproduction-and-growth.png',dpi=180,bbox_inches='tight');plt.close(fig)
print('Created figures',[(p.name,p.stat().st_size) for p in out.iterdir()])
