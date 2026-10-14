import sys, pickle; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from util import *; from curves import *
from scipy import ndimage as ndi
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
img, finaldir = sys.argv[1:3]
a = load(img); R,G,B = a[...,0],a[...,1],a[...,2]
red=(R>150)&(G<130)&(B<130)
blk=(a.mean(2)<100)&~red&(R<120)
white=a.mean(2)>160
XT=np.array([85,204,322.5,441,560,678.]); CT=np.array([0,20,40,60,80,100.]); pc=np.polyfit(XT,CT,1)
YT=np.array([10,105,200,295,389.5,484.]); VT=np.array([3,2,1,0,-1,-2.]); pm=np.polyfit(YT,VT,1)
cx=lambda x: np.polyval(pc,x); my=lambda y: np.polyval(pm,y); xc=lambda c:(np.asarray(c)-pc[1])/pc[0]; ym=lambda m:(np.asarray(m)-pm[1])/pm[0]
print('c fit',pc,CT-cx(XT)); print('m fit',pm,VT-my(YT))
# legend box: x 251-459, y 122-281 -> exclude
blk_l=blk.copy(); blk_l[118:286,248:463]=False; red_l=red.copy(); red_l[118:286,248:463]=False
blk_l[:, :88]=False; blk_l[:, 676:]=False; blk_l[:14,:]=False
# ticks inward on left/right frame (to ~x 102 / 662) -> remove short horizontal tick runs: exclude tick rows near frames
for yt in list(YT)+[57,152,247,342,437]:
    blk_l[int(yt)-2:int(yt)+3, 86:103]=False; blk_l[int(yt)-2:int(yt)+3, 662:677]=False
for xt in list(XT)+[145,263,381,500,619]:
    blk_l[12:26, int(xt)-2:int(xt)+3]=False; blk_l[470:483, int(xt)-2:int(xt)+3]=False
xs=list(range(88,676))
fe = track(blk_l, xs, 14, 200, y_start=ym(2.25), win=6, maxrun=8, maxgap=80)
cr = track(blk_l, xs, 280, 430, y_start=ym(-1.32), win=6, maxrun=14, maxgap=80)
# DLM symbols: white cores enclosed by red rings
lab,n=ndi.label(white)
sym=[]
for i,sl in enumerate(ndi.find_objects(lab)):
    ys,xs_=sl; h=ys.stop-ys.start; w=xs_.stop-xs_.start; area=(lab[sl]==i+1).sum()
    if not (2<=h<=14 and 2<=w<=14 and area>=15): continue
    cy=(ys.start+ys.stop-1)/2; cxp=(xs_.start+xs_.stop-1)/2
    if 248<cxp<463 and 118<cy<286: continue
    patch=red[max(0,int(cy)-11):int(cy)+12, max(0,int(cxp)-11):int(cxp)+12]
    if patch.mean()<0.15: continue
    note='white core'
    if cxp<92: cxp=85.0; note='core cut by left frame; x set to frame (c=0)'
    if cxp>670: continue
    kind='M_Fe DLM (open red circle)' if cy<250 else 'M_Cr DLM (open red square)'
    sym.append((kind,cxp,cy,w,h,area,note))
sym.append(('M_Cr DLM (open red square)',678.0,295.0,0,0,0,'core split by right frame (pieces at x 674.5 and 680.5); x set to frame (c=100)'))
sym.sort(key=lambda t:(t[0],t[1]))
for s in sym: print(s[0], round(float(cx(s[1])),2), round(float(my(s[2])),3), s[3:])
pickle.dump(dict(fe=fe,cr=cr,sym=sym,pc=pc,pm=pm), open(sys.argv[3],'wb'))
print('fe',len(fe),'cr',len(cr))
