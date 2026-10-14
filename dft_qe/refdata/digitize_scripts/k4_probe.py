import sys; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from holes import *
a=np.asarray(Image.open(sys.argv[1]).convert('L')); ink=a>128
hs=find_holes(ink,(100,796,14,547))
for h in sorted(hs,key=lambda h:(h['cx'],h['cy'])):
    print(f"{h['cx']:6.1f} {h['cy']:6.1f} w{h['w']:2d} h{h['h']:2d} a{h['area']:3d} fill{h['fill']:.2f} dot{int(h['dot'])}")
