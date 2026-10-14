import sys; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from util import *
a = load(sys.argv[1]); x0,x1,y0,y1 = map(int, sys.argv[2:6])
R,G,B = a[...,0],a[...,1],a[...,2]
print('     '+''.join(str((x//10)%10) if x%10==0 else ' ' for x in range(x0,x1)))
print('     '+''.join(str(x%10) for x in range(x0,x1)))
for y in range(y0,y1):
    s=''
    for x in range(x0,x1):
        r,g,b = R[y,x],G[y,x],B[y,x]
        if r>150 and g<110 and b<110: c='r'
        elif g>150 and r<120 and b<120: c='g'
        elif b>150 and r<120 and g<150: c='b'
        elif (r+g+b)/3<90: c='#'
        elif (r+g+b)/3<170: c='+'
        else: c='.'
        s+=c
    print(f'{y:4d} '+s)
