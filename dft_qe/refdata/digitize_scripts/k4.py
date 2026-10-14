import sys, pickle; sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from tmatch import *
from PIL import Image
a=np.asarray(Image.open(sys.argv[1]).convert('L')); ink=a>128
# templates (centre pixel chosen from isolated markers; circle centre 331.5,102.5 -> use 331/102 with even/odd handling)
Tc = make_template(ink, 331, 102, 10, 10)   # Fe circle (x 322-341, y 94-111)
Ts = make_template(ink, 535, 257, 10, 10)   # average square (x 528-543, y 249-265)
Td = make_template(ink, 447, 281, 10, 10)   # Cr diamond (x 439-456, y 273-290)
res={}
for name,T in [('circle',Tc),('square',Ts),('diamond',Td)]:
    sc=match(ink,T); res[name]=sc
    print(name, 'self', sc[102,331] if name=='circle' else (sc[257,535] if name=='square' else sc[281,447]))
pickle.dump(res, open(sys.argv[2],'wb'))
