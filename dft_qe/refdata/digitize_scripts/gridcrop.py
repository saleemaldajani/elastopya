import sys; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from util import *
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
img,out = sys.argv[1],sys.argv[2]; x0,x1,y0,y1 = map(int, sys.argv[3:7]); st=int(sys.argv[7]) if len(sys.argv)>7 else 5
a = load(img).astype(np.uint8)
fig,ax=plt.subplots(figsize=(12,12*(y1-y0)/(x1-x0)+0.5),dpi=100)
ax.imshow(a[y0:y1,x0:x1], extent=(x0-0.5,x1-0.5,y1-0.5,y0-0.5), interpolation='nearest')
ax.set_xticks(range(x0 - x0%st + st, x1, st)); ax.set_yticks(range(y0 - y0%st + st, y1, st))
ax.tick_params(labelsize=6); ax.grid(color='c', lw=0.3, alpha=0.6)
plt.xticks(rotation=90); fig.tight_layout(); fig.savefig(out)
