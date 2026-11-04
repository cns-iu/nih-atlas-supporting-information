from pathlib import Path
from xml.sax.saxutils import escape

from PIL import Image, ImageDraw, ImageFont

ROOT = Path('/Users/yongxinkong/dev/nih-map-supporting-information')
OUT = ROOT / 'outputs' / 'tables'
SVG_PATH = OUT / 'table_4_conditional_bars.svg'
PNG_PATH = OUT / 'table_4_conditional_bars.png'

headers = [
    'Atlas Program',
    'Internal Collaboration Activity Level',
    'Percent of PIs',
    'Percent of Researchers',
]
rows = [
    ('HTAN', 'Highly Active (10+)', 3.17, 0.06),
    ('HTAN', 'Very Active (5-9)', 0.00, 0.32),
    ('HTAN', 'Moderately Active (2-4)', 22.22, 1.56),
    ('HTAN', 'Slightly Active (1)', 22.22, 7.96),
    ('HTAN', 'Inactive (0)', 52.38, 90.10),
    ('HuBMAP', 'Highly Active (10+)', 5.71, 0.06),
    ('HuBMAP', 'Very Active (5-9)', 18.57, 0.90),
    ('HuBMAP', 'Moderately Active (2-4)', 30.00, 6.14),
    ('HuBMAP', 'Slightly Active (1)', 24.29, 16.55),
    ('HuBMAP', 'Inactive (0)', 21.43, 76.35),
    ('SenNet', 'Highly Active (10+)', 18.84, 0.22),
    ('SenNet', 'Very Active (5-9)', 17.39, 0.97),
    ('SenNet', 'Moderately Active (2-4)', 33.33, 10.99),
    ('SenNet', 'Slightly Active (1)', 23.19, 31.00),
    ('SenNet', 'Inactive (0)', 7.25, 56.82),
]

col_widths = [128, 278, 225, 225]
header_h = 27
row_h = 27
width = sum(col_widths)
height = header_h + row_h * len(rows)
max_pis = max(r[2] for r in rows)
max_researchers = max(r[3] for r in rows)

def esc(value):
    return escape(str(value))

def svg_text(x, y, text, anchor='start', bold=False, size=13):
    weight = '700' if bold else '400'
    return f'<text x="{x}" y="{y}" dominant-baseline="middle" text-anchor="{anchor}" font-family="Calibri, Arial, sans-serif" font-size="{size}px" font-weight="{weight}" fill="#000000">{esc(text)}</text>'

svg = [
    f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
    '<defs>',
    '<linearGradient id="piGradient" x1="0" y1="0" x2="1" y2="0">',
    '<stop offset="0%" stop-color="#2F75B5"/>',
    '<stop offset="100%" stop-color="#FFFFFF"/>',
    '</linearGradient>',
    '<linearGradient id="researcherGradient" x1="0" y1="0" x2="1" y2="0">',
    '<stop offset="0%" stop-color="#2E8B57"/>',
    '<stop offset="100%" stop-color="#FFFFFF"/>',
    '</linearGradient>',
    '</defs>',
    '<rect width="100%" height="100%" fill="#FFFFFF"/>',
]

x = 0
for w in col_widths:
    svg.append(f'<rect x="{x}" y="0" width="{w}" height="{header_h}" fill="#FFFFFF" stroke="#D9E1DC"/>')
    x += w
svg.append(svg_text(col_widths[0] / 2, header_h / 2, headers[0], anchor='middle', bold=True))
svg.append(svg_text(col_widths[0] + col_widths[1] / 2, header_h / 2, 'Internal Collaboration Activity Level', anchor='middle', bold=True))
svg.append(svg_text(sum(col_widths[:2]) + col_widths[2] / 2, header_h / 2, headers[2], anchor='middle', bold=True))
svg.append(svg_text(sum(col_widths[:3]) + col_widths[3] / 2, header_h / 2, headers[3], anchor='middle', bold=True))

for i, row in enumerate(rows):
    y = header_h + i * row_h
    x = 0
    for w in col_widths:
        svg.append(f'<rect x="{x}" y="{y}" width="{w}" height="{row_h}" fill="#FFFFFF" stroke="#D9E1DC"/>')
        x += w
    svg.append(svg_text(10, y + row_h / 2, row[0]))
    svg.append(svg_text(col_widths[0] + 10, y + row_h / 2, row[1]))
    for j, (value, max_value) in enumerate(((row[2], max_pis), (row[3], max_researchers))):
        cell_x = sum(col_widths[:2 + j])
        cell_w = col_widths[2 + j]
        bar_w = max(0, (cell_w - 8) * value / max_value)
        if bar_w > 0:
            gradient_id = 'piGradient' if j == 0 else 'researcherGradient'
            svg.append(f'<rect x="{cell_x + 4}" y="{y + 4}" width="{bar_w:.2f}" height="{row_h - 8}" fill="url(#{gradient_id})"/>')
        svg.append(svg_text(cell_x + cell_w - 10, y + row_h / 2, f'{value:.2f}', anchor='end'))

svg.append('</svg>')
SVG_PATH.write_text(''.join(svg), encoding='utf-8')

# Render a PNG using the same geometry. Calibri is specified in the SVG;
# this runtime does not ship Calibri, so Arial is used as the closest local fallback.
scale = 2
img = Image.new('RGB', (width * scale, height * scale), 'white')
draw = ImageDraw.Draw(img)
font_path = '/System/Library/Fonts/Supplemental/Arial.ttf'
bold_path = '/System/Library/Fonts/Supplemental/Arial Bold.ttf'
font = ImageFont.truetype(font_path, 13 * scale)
bold_font = ImageFont.truetype(bold_path, 13 * scale)

def rect(box, fill='white', outline='#D9E1DC', width_px=1):
    draw.rectangle(tuple(int(v * scale) for v in box), fill=fill, outline=outline, width=max(1, int(width_px * scale)))

def text(x, y, value, font_obj, anchor='la'):
    draw.text((int(x * scale), int(y * scale)), value, fill='black', font=font_obj, anchor=anchor)

for i, w in enumerate(col_widths):
    x0 = sum(col_widths[:i])
    rect((x0, 0, x0 + w, header_h))
text(col_widths[0] / 2, header_h / 2, headers[0], bold_font, 'mm')
text(col_widths[0] + col_widths[1] / 2, header_h / 2, 'Internal Collaboration Activity Level', bold_font, 'mm')
text(sum(col_widths[:2]) + col_widths[2] / 2, header_h / 2, headers[2], bold_font, 'mm')
text(sum(col_widths[:3]) + col_widths[3] / 2, header_h / 2, headers[3], bold_font, 'mm')

def blend(c0, c1, t):
    return tuple(round(a + (b - a) * t) for a, b in zip(c0, c1))

pi_start = (47, 117, 181)
researcher_start = (46, 139, 87)
white = (255, 255, 255)
for i, row in enumerate(rows):
    y = header_h + i * row_h
    x = 0
    for w in col_widths:
        rect((x, y, x + w, y + row_h))
        x += w
    text(10, y + row_h / 2, row[0], font, 'lm')
    text(col_widths[0] + 10, y + row_h / 2, row[1], font, 'lm')
    for j, (value, max_value) in enumerate(((row[2], max_pis), (row[3], max_researchers))):
        cell_x = sum(col_widths[:2 + j])
        cell_w = col_widths[2 + j]
        bar_w = max(0, (cell_w - 8) * value / max_value)
        for px in range(int(bar_w * scale)):
            t = px / max(1, int(bar_w * scale) - 1)
            start = pi_start if j == 0 else researcher_start
            color = blend(start, white, t)
            draw.line((int((cell_x + 4) * scale + px), int((y + 4) * scale), int((cell_x + 4) * scale + px), int((y + row_h - 4) * scale)), fill=color, width=1)
        text(cell_x + cell_w - 10, y + row_h / 2, f'{value:.2f}', font, 'rm')

img.save(PNG_PATH)
print(SVG_PATH)
print(PNG_PATH)
