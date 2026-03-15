"""Generate analysis flow diagram for the hitch load calculator."""

import math
from PIL import Image, ImageDraw, ImageFont

# --- Config ---
SCALE = 2
TARGET_W, TARGET_H = 900, 520

def s(val):
    return int(val * SCALE)

# Colors - muted professional palette on dark background
BG = (18, 22, 28)
CARD_BG = (28, 34, 42)
CARD_ACCENT = (36, 42, 52)
BORDER = (55, 62, 72)
TEXT = (200, 208, 216)
TEXT_DIM = (130, 138, 148)
TITLE_COLOR = (235, 240, 245)
BLUE = (80, 150, 230)
TEAL = (70, 180, 170)
AMBER = (220, 160, 50)
CORAL = (210, 95, 80)
BLUE_DIM = (50, 90, 140)
TEAL_DIM = (40, 110, 105)
AMBER_DIM = (135, 100, 30)
CORAL_DIM = (130, 60, 50)

# Fonts
FONT_BOLD = "C:/Windows/Fonts/segoeuib.ttf"
FONT_REG = "C:/Windows/Fonts/segoeui.ttf"

title_font = ImageFont.truetype(FONT_BOLD, s(22))
heading_font = ImageFont.truetype(FONT_BOLD, s(16))
body_font = ImageFont.truetype(FONT_REG, s(13))
small_font = ImageFont.truetype(FONT_REG, s(11))
label_font = ImageFont.truetype(FONT_BOLD, s(11))

# Canvas
img = Image.new("RGB", (s(TARGET_W), s(TARGET_H)), BG)
draw = ImageDraw.Draw(img)


def make_node(x, y, w, h):
    return {
        "rect": (s(x), s(y), s(x + w), s(y + h)),
        "top": (s(x + w / 2), s(y)),
        "bottom": (s(x + w / 2), s(y + h)),
        "left": (s(x), s(y + h / 2)),
        "right": (s(x + w), s(y + h / 2)),
        "center": (s(x + w / 2), s(y + h / 2)),
    }


def draw_arrow(start, end, color, width=2, head_size=9):
    w = s(width)
    hs = s(head_size)
    draw.line([start, end], fill=color, width=w)
    angle = math.atan2(end[1] - start[1], end[0] - start[0])
    x, y = end
    for da in [2.6, -2.6]:
        ax = x - hs * math.cos(angle + da)
        ay = y - hs * math.sin(angle + da)
        draw.line([(x, y), (int(ax), int(ay))], fill=color, width=w)


def draw_card(node, fill, border_color, radius=10):
    draw.rounded_rectangle(node["rect"], radius=s(radius), fill=fill, outline=border_color, width=s(1))


def draw_text_centered(x, y, text, font, color):
    draw.text((s(x), s(y)), text, font=font, fill=color, anchor="mm")


def draw_color_bar(node, color, side="top"):
    x1, y1, x2, y2 = node["rect"]
    if side == "top":
        draw.rounded_rectangle((x1, y1, x2, y1 + s(4)), radius=s(2), fill=color)
    elif side == "left":
        draw.rounded_rectangle((x1, y1 + s(5), x1 + s(4), y2 - s(5)), radius=s(2), fill=color)


# === Layout ===

# Row 1: Input Parameters (centered)
input_node = make_node(310, 30, 280, 65)

# Row 2: Static Analysis (two cards side by side)
weld_node = make_node(80, 155, 230, 95)
torsion_node = make_node(590, 155, 230, 95)

# Row 3: Deflection (wide card, center)
defl_node = make_node(175, 310, 550, 95)

# Row 4: Dynamic Analysis (centered)
dyn_node = make_node(310, 455, 280, 50)

# --- Draw connecting lines first ---

# Input -> Weld
draw_arrow(input_node["bottom"], weld_node["top"], BLUE_DIM)
# Input -> Torsion
draw_arrow(input_node["bottom"], torsion_node["top"], BLUE_DIM)

# Weld -> Deflection
draw_arrow(weld_node["bottom"], (s(195), defl_node["top"][1]), TEAL_DIM)
# Torsion -> Deflection
draw_arrow(torsion_node["bottom"], (s(705), defl_node["top"][1]), TEAL_DIM)

# Deflection -> Dynamic
draw_arrow(defl_node["bottom"], dyn_node["top"], AMBER_DIM)


# --- Draw cards ---

# Input Parameters
draw_card(input_node, CARD_BG, BORDER)
draw_color_bar(input_node, BLUE)
draw_text_centered(450, 52, "Input Parameters", heading_font, TITLE_COLOR)
draw_text_centered(450, 73, "Load  |  Material  |  Geometry", small_font, TEXT_DIM)

# Static Analysis label
draw_text_centered(450, 130, "STATIC ANALYSIS", label_font, TEXT_DIM)

# Weld Shear
draw_card(weld_node, CARD_BG, BORDER)
draw_color_bar(weld_node, TEAL)
draw_text_centered(195, 178, "Weld Shear", heading_font, TITLE_COLOR)
draw_text_centered(195, 200, "Weld-as-a-line method", small_font, TEXT_DIM)
draw_text_centered(195, 220, "E70 electrode, 3/16\" fillet", small_font, TEXT_DIM)
draw_text_centered(195, 238, "FOS check vs. allowable", small_font, TEAL)

# Crossbar Torsion
draw_card(torsion_node, CARD_BG, BORDER)
draw_color_bar(torsion_node, TEAL)
draw_text_centered(705, 178, "Crossbar Torsion", heading_font, TITLE_COLOR)
draw_text_centered(705, 200, "Bredt-Batho thin-wall theory", small_font, TEXT_DIM)
draw_text_centered(705, 220, "2x2\" square tube, 0.25\" wall", small_font, TEXT_DIM)
draw_text_centered(705, 238, "Shear stress vs. Von Mises", small_font, TEAL)

# Deflection
draw_card(defl_node, CARD_BG, BORDER)
draw_color_bar(defl_node, AMBER)
draw_text_centered(450, 333, "Deflection Analysis (Superposition)", heading_font, TITLE_COLOR)
draw_text_centered(450, 358, "Twist contribution (torsion angle)  +  Bending contribution (cantilever)", body_font, TEXT_DIM)
draw_text_centered(450, 380, "Combined tip angle and vertical drop at bike rack COG", small_font, AMBER)

# Dynamic Analysis
draw_card(dyn_node, CARD_ACCENT, BORDER)
draw_color_bar(dyn_node, CORAL)
draw_text_centered(450, 475, "Dynamic Analysis  (4G Pothole)", heading_font, TITLE_COLOR)
draw_text_centered(450, 495, "All checks re-evaluated at impact load", small_font, CORAL)

# --- Subtle title ---
draw_text_centered(450, 510, "Hitch Carrier Load & Deflection Calculator", small_font, (70, 75, 82))

# --- Save ---
final = img.resize((TARGET_W, TARGET_H), Image.Resampling.LANCZOS)
final.save("C:/Users/aaron/Documents/c-projects/_temp-hitch-calc/docs/analysis-diagram.png", quality=95)
print("Saved analysis-diagram.png")
