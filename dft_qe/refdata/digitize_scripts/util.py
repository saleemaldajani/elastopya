import numpy as np
from PIL import Image
def load(fn):
    return np.asarray(Image.open(fn).convert('RGB')).astype(float)
def crop(fn, box, out, scale=3):
    im = Image.open(fn).convert('RGB').crop(box)
    im = im.resize((im.width*scale, im.height*scale), Image.NEAREST)
    im.save(out)
