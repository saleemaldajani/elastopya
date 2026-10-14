import sys, fitz, os
pdf, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
d = fitz.open(pdf)
for pno, page in enumerate(d):
    for img in page.get_images(full=True):
        xref = img[0]
        info = d.extract_image(xref)
        rects = page.get_image_rects(xref)
        fn = f"{out}/p{pno+1}_x{xref}.{info['ext']}"
        open(fn,'wb').write(info['image'])
        print(pno+1, xref, info['width'], info['height'], info['ext'], img[2:5], [tuple(round(v) for v in r) for r in rects])
