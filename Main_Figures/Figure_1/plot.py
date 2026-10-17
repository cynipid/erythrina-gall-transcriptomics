"""Rebuild the preferred Figure 1 layout and its matching 600 dpi PNG.

Panel A retains the exact supplied photograph and crop. No analyses are run.
Dependencies: reportlab, Pillow, PyMuPDF. Fonts: standard PDF Helvetica.
Run from any working directory; the adjacent photo is the only input asset.
"""
from pathlib import Path
from PIL import Image
import fitz
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics

import argparse
ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--outdir", type=Path, default=ROOT)
args=parser.parse_args()
args.outdir.mkdir(parents=True, exist_ok=True)
SOURCE=ROOT/'Figure_1_photo.png'
W,H=2046,1980
PDF_W=180/25.4*72
SCALE=PDF_W/W
PDF=args.outdir/'Figure_1.pdf'
PNG=args.outdir/'Figure_1.png'
FONTS={'regular':'Helvetica','bold':'Helvetica-Bold','italic':'Helvetica-Oblique'}

c=canvas.Canvas(str(PDF),pagesize=(PDF_W,H*SCALE),pageCompression=1)
c.setTitle('Figure 1 Study system and experimental design')
c.setSubject('Preferred original layout; larger labels; one-tree and petiole replication explicit; attribution retained in caption')
c.scale(SCALE,SCALE)
INK='#27313e';GRAY='#697585';LINE='#5f6a78'
def fill(color):c.setFillColor(HexColor(color))
def stroke(color,width=3):c.setStrokeColor(HexColor(color));c.setLineWidth(width)
def rr(x,y,w,h,r,face,edge,width=3):
    fill(face);stroke(edge,width);c.roundRect(x,H-y-h,w,h,r,stroke=1,fill=1)
def line(points,color=LINE,width=4):
    stroke(color,width);p=c.beginPath();p.moveTo(points[0][0],H-points[0][1])
    for x,y in points[1:]:p.lineTo(x,H-y)
    c.drawPath(p,stroke=1)
def label(text,x,y,size=28,font='regular',color=INK,maxwidth=None):
    font=FONTS.get(font,font)
    if maxwidth:
        size=min(size,size*maxwidth/pdfmetrics.stringWidth(text,font,size))
    fill(color);c.setFont(font,size);c.drawCentredString(x,H-y,text)
def ellipse(cx,cy,rx,ry,face,edge,width=3):
    fill(face);stroke(edge,width);c.ellipse(cx-rx,H-cy-ry,cx+rx,H-cy+ry,stroke=1,fill=1)
def path(points,edge,face=None,width=3,close=False):
    stroke(edge,width)
    if face:fill(face)
    p=c.beginPath()
    for op,pts in points:
        xy=[]
        for i in range(0,len(pts),2):xy.extend([pts[i],H-pts[i+1]])
        getattr(p,op)(*xy)
    if close:p.close()
    c.drawPath(p,stroke=1,fill=bool(face))

# Preserve the exact photographic panel from the supplied composite.
with Image.open(SOURCE) as source:
    photo=source.convert('RGB')
    c.drawImage(ImageReader(photo),0,H-820,width=W,height=820,mask='auto')
# Panel marker only; the underlying photographic crop is unchanged.
rr(30,30,84,84,11,'#ffffff',INK,3)
label('a',72,94,58,'bold')
line([(0,838),(W,838)],'#d6dde5',2)
rr(30,884,84,84,11,'#ffffff',INK,3)
label('b',72,946,58,'bold')
rr(466,895,1114,93,18,'#f3f5f7','#d0d7df',2.6)
label('Single host tree',1023,955,44,'bold',maxwidth=1040)
line([(1023,988),(1023,1025)])
line([(511,1025),(1534,1025)])
for x in [511,1534]:
    path([('moveTo',[x-26,1020]),('lineTo',[x+26,1020]),('lineTo',[x,1067])],LINE,LINE,2,True)

for x in [62,1085]:rr(x,1071,899,786,25,'#fafbfc','#c4cbd3',3.3)
label('Insect developmental contrast',511,1133,38,'bold')
label('Quadrastichus erythrinae',511,1176,33,'italic',GRAY)
label('Host-petiole contrast',1534,1133,38,'bold')
label('Erythrina variegata',1534,1176,33,'italic',GRAY)

cards=[(108,'#e8f6f1','#009e73','Larva',True),(534,'#faedf5','#cc79a7','Pupa',True),(1130,'#e8f2f8','#0072b2','Galled petiole',False),(1556,'#fceee6','#d55e00','Ungalled petiole',False)]
for x,bg,edge,name,insect in cards:
    rr(x,1213,381,398,23,bg,edge,4.2)
    cx=x+190.5
    label(name,cx,1438,36,'bold')
    label('3 pooled samples' if insect else '3 petiole replicates',cx,1493,33,'bold')
    label('10 individuals' if insect else '1 whole petiole',cx,1541,33)
    label('per pool' if insect else 'per replicate',cx,1584,33,color=GRAY)

# Larval schematic retained as overlapping curved segments.
for x,y,rx,ry in [(247,1351,57,49),(304,1345,54,47),(357,1350,51,45),(403,1360,44,40)]:
    ellipse(x,y,rx,ry,'#ffffff','#008669',4.1)

# Pupal schematic is slightly reduced and raised: its bottom is 1378,
# leaving about 30 design units before the top of the aligned Pupa label.
cx,cy,k=724.5,1308,.85
def local(coords):
    return [v*k+(cx if i%2==0 else cy) for i,v in enumerate(coords)]
path([('moveTo',local([0,-71])),('curveTo',local([-54,-73,-72,-27,-59,12])),('curveTo',local([-48,53,-22,74,0,82])),('curveTo',local([22,74,48,53,59,12])),('curveTo',local([72,-27,54,-73,0,-71]))],'#b34c91','#ffffff',3.8,True)
ellipse(cx,cy-44*k,24*k,18*k,'#faedf5','#b34c91',3.4)
for y,half in [(-22,38),(7,45),(36,39)]:
    path([('moveTo',local([-half,y])),('curveTo',local([-half/3,y+13,half/3,y+13,half,y]))],'#b34c91',None,3.3)

def petiole(cx,galled):
    line([(cx-84,1397),(cx+79,1281)],'#2b8252',10)
    path([('moveTo',[cx,1337]),('curveTo',[cx-32,1292,cx-62,1271,cx-91,1275]),('curveTo',[cx-88,1318,cx-51,1351,cx,1337])],'#42864e','#90c979',3.3,True)
    path([('moveTo',[cx+20,1320]),('curveTo',[cx+49,1275,cx+82,1259,cx+126,1267]),('curveTo',[cx+111,1308,cx+74,1332,cx+20,1320])],'#42864e','#90c979',3.3,True)
    if galled:
        for dx,dy,r in [(-44,1373,16),(-23,1353,24),(0,1337,24),(23,1328,18)]:ellipse(cx+dx,dy,r,r,'#b6d77e','#5a843b',3.2)
petiole(1320.5,True);petiole(1746.5,False)

for x,center in [(108,511),(534,511),(1130,1534),(1556,1534)]:
    cx=x+190.5
    line([(cx,1611),(cx,1662),(center+(-82 if cx<center else 82),1688)],LINE,4)
for x,center,top,sub in [(193,511,'6 insect RNA-seq libraries','Genome-guided analysis'),(1216,1534,'6 plant RNA-seq libraries','De novo transcriptome')]:
    rr(x,1688,636,123,20,'#ffffff','#7c8794',3.4)
    label(top,center,1740,35,'bold')
    label(sub,center,1783,33,color=GRAY,maxwidth=582)
line([(511,1811),(511,1879),(750,1879)],LINE,4)
line([(1534,1811),(1534,1879),(1296,1879)],LINE,4)
rr(750,1845,546,88,20,INK,INK,1)
label('12 RNA-seq libraries',1023,1902,36,'bold','#ffffff')
c.showPage();c.save()

doc=fitz.open(PDF);page=doc[0]
assert len(page.get_drawings())>=35
assert len(page.search_for('Pupa'))==1
assert '12 RNA-seq libraries' in page.get_text()
pix=page.get_pixmap(dpi=600,alpha=False)
PNG.write_bytes(pix.tobytes('png'))
assert PNG.stat().st_size>1000000
print('PDF:',PDF)
print('PNG:',PNG,pix.width,pix.height,'600 dpi')
print('Vector drawing objects:',len(page.get_drawings()),'embedded photo images:',len(page.get_images()))
