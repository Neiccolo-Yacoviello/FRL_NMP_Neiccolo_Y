"""Draws rocket.png (the sprite used by live_grapher.py). Run once, or swap in your own PNG
(nose pointing up, transparent background)."""
from PIL import Image, ImageDraw

S = 4                      # supersample factor for smooth edges
W, H = 64, 240             # final size in pixels


def make_rocket(path="rocket.png"):
    img = Image.new("RGBA", (W * S, H * S), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    p = lambda pts: [(x * S, y * S) for x, y in pts]
    outline = (30, 30, 40, 255)
    lw = 2 * S

    # fins (behind body)
    d.polygon(p([(18, 165), (2, 205), (2, 228), (18, 212)]), fill=(200, 40, 40, 255), outline=outline, width=lw)
    d.polygon(p([(46, 165), (62, 205), (62, 228), (46, 212)]), fill=(200, 40, 40, 255), outline=outline, width=lw)
    # body tube
    d.rectangle(p([(18, 58), (46, 214)]), fill=(235, 235, 240, 255), outline=outline, width=lw)
    # body stripe + window
    d.rectangle(p([(18, 120), (46, 128)]), fill=(40, 70, 160, 255))
    d.ellipse(p([(26, 80), (38, 92)]), fill=(90, 170, 230, 255), outline=outline, width=lw)
    # nose cone (ogive-ish)
    nose =[(18, 58), (20, 38), (24, 22), (28, 10), (32, 4), (36, 10), (40, 22), (44, 38), (46, 58)]
    d.polygon(p(nose), fill=(200, 40, 40, 255), outline=outline, width=lw)
    # nozzle
    d.polygon(p([(26, 214), (38, 214), (40, 222), (24, 222)]), fill=(70, 70, 80, 255), outline=outline, width=lw)

    img.resize((W, H), Image.LANCZOS).save(path)
    return path


if __name__ == "__main__":
    print("wrote", make_rocket())
