"""Rebuild the preview with clipped artwork below its heading/caption layer."""
import io
import json
from pathlib import Path

import fitz
import matplotlib
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen import canvas
from reportlab.platypus import Paragraph

ROOT = Path(__file__).resolve().parent
DATA = json.loads((ROOT / 'captions.json').read_text())
FONTROOT = Path(matplotlib.get_data_path()) / 'fonts/ttf'
for name, filename in [('Serif', 'DejaVuSerif.ttf'), ('Serif-Bold', 'DejaVuSerif-Bold.ttf'),
                       ('Serif-Italic', 'DejaVuSerif-Italic.ttf'), ('Serif-BoldItalic', 'DejaVuSerif-BoldItalic.ttf')]:
    pdfmetrics.registerFont(TTFont(name, str(FONTROOT / filename)))
pdfmetrics.registerFontFamily('Serif', normal='Serif', bold='Serif-Bold', italic='Serif-Italic', boldItalic='Serif-BoldItalic')
W, H = 595.276, 841.89
MARGIN, WIDTH = 48, W-96
STYLE = ParagraphStyle('caption', fontName='Serif', fontSize=10.5, leading=14.1, textColor='black')
TITLE = ParagraphStyle('title', parent=STYLE, fontSize=14, leading=20, alignment=TA_CENTER)
CENTER = ParagraphStyle('center', parent=STYLE, fontSize=12, leading=17, alignment=TA_CENTER)


def paragraph(c, text, style, top, width=WIDTH, x=MARGIN):
    p = Paragraph(text, style)
    _, h = p.wrap(width, H)
    p.drawOn(c, x, top-h)
    return h


def footer(c, n):
    c.setFont('Serif', 9)
    c.drawCentredString(W/2, 30, str(n))


buffer = io.BytesIO()
c = canvas.Canvas(buffer, pagesize=(W, H))
c.setTitle(DATA['title'])
c.setAuthor(DATA['author'])
paragraph(c, '<b>'+DATA['title']+'</b>', TITLE, H-100)
paragraph(c, DATA['paper_title_html'], TITLE, H-155)
paragraph(c, DATA['author'], CENTER, H-250)
paragraph(c, DATA['range_label'], CENTER, H-292)
top = H-360
for figure in DATA['figures']:
    top -= paragraph(c, figure['index_label'], STYLE, top, WIDTH-84, MARGIN+42) + 8
footer(c, 1)
c.showPage()

placements = []
for i, figure in enumerate(DATA['figures'], 1):
    asset = ROOT / figure['folder'] / (figure['stem']+'.pdf')
    with fitz.open(asset) as source:
        iw, ih = source[0].rect.width, source[0].rect.height
    cap = Paragraph(figure['caption_html'], STYLE)
    _, ch = cap.wrap(WIDTH, H)
    maxheight = min(530, H-95-50-24-ch)
    scale = min(WIDTH/iw, maxheight/ih)
    fw, fh = iw*scale, ih*scale
    left, top = (W-fw)/2, H-75
    bottom = top-fh
    c.setFont('Serif-Bold', 12)
    c.drawString(MARGIN, H-48, figure['label'])
    cap.drawOn(c, MARGIN, bottom-20-ch)
    assert bottom-20-ch > 48
    placements.append((i, asset, [left, H-top, left+fw, H-bottom]))
    footer(c, i+1)
    c.showPage()
c.save()
pdf = fitz.open(stream=buffer.getvalue(), filetype='pdf')
for page, asset, rect in placements:
    with fitz.open(asset) as source:
        target = fitz.Rect(rect)
        assert target.y0 >= 75, 'Artwork must remain below the heading area'
        pdf[page].show_pdf_page(
            target, source, 0, overlay=False, clip=source[0].rect
        )
output = ROOT / DATA['output_name']
pdf.save(output, garbage=4, deflate=True)
print(output.name, len(pdf), 'pages')
