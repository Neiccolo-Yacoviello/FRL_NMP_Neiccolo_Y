"""
Live grapher for the 3DOF rocket sim.

The sim steps forward in (scaled) real time while the plots fill in, and a rocket
sprite on the left shows what the rocket is doing: position, pitch, motor flame,
air brake deployment and parachutes.

Usage:
    python live_grapher.py               # uses RUNTIME / PLAYBACK_SPEED below
    python live_grapher.py 30            # graph 30 s of flight
    python live_grapher.py 30 2          # graph 30 s of flight at 2x speed
    python live_grapher.py full 5        # whole flight (until touchdown) at 5x

Needs rocket.png in the same folder (make_rocket_png.py draws one; if it's
missing it gets generated automatically). Put this file next to
model_functions.py, model_constants.py, controller.py and the CSVs.
"""
import sys
import time as wall
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation 
from matplotlib.offsetbox import OffsetImage, AnnotationBbox 
from PIL import Image, ImageDraw 

# model.py runs the sim when imported, so pull from model_functions instead
from model_functions import *

# ---------------- SETTINGS ----------------
RUNTIME = 190.0           # seconds of flight to simulate & graph (None = until touchdown)
PLAYBACK_SPEED = 3    # 1.0 = real time, 2.0 = twice as fast, 0.5 = slow-mo
INITIAL_STEP = 0.01      # starting RK45 step size (same as model.py)
MAX_STEP = 0.1           # cap on step size so descent animates smoothly (ascent is already capped at controller_step)
FPS = 120                 # animation frame rate
ROCKET_PNG = Path(__file__).parent / "rocket.png"
SPRITE_ZOOM = 0.38       # on-screen size of the rocket sprite
BRAKE_VISUAL_SCALE = 2   # air brakes max out at 20 deg; drawn 3x so they're visible
# ------------------------------------------

LOG_KEYS = ("t", "y", "z", "v_y", "v_z", "pitch", "v_pitch", "m",
            "speed", "control_value", "predicted_apogee", "thrust")
MAX_THRUST = float(thrust_data[:, 1].max())


class LiveSim:
    """Same integration loop as model() in model.py, but advanced a bit at a time."""

    def __init__(self, max_time):
        reset_controller()
        self.max_time = max_time
        self.state = np.array([0, 0, 0, 0, 0, 0, m0 + MCV["mass_offset"]], dtype=float)
        self.prev = self.state
        self.step_size = INITIAL_STEP
        self.t = 0.0
        self.done = False
        self.end_reason = ""
        self.log = {k: [] for k in LOG_KEYS}
        self._record()

    def _record(self):
        y, z, v_y, v_z, pitch, v_pitch, m = self.state
        pa = controller_state["predicted_apogee"]
        for k, v in zip(LOG_KEYS, (self.t, y, z, v_y, v_z, np.degrees(pitch), np.degrees(v_pitch), m,
                                    np.hypot(v_y, v_z), controller_state["control_value"],
                                    np.nan if pa is None else pa, thrust_LUT(self.t))):
            self.log[k].append(v)

    def advance_to(self, target_t):
        """Integrate until sim time reaches target_t (or the flight ends)."""
        while not self.done and self.t < target_t:
            if self.max_time is not None:
                self.step_size = min(self.step_size, self.max_time - self.t)

            b4, b5 = rk_45_step(self.t, self.step_size, self.state)
            accepted, err = check_step(b4, b5, self.prev)
            if accepted:
                self.prev, self.state = self.state, b5
                self.t += self.step_size
                update_controller(self.t, self.state)
                self._record()

            self.step_size = min(adjust_step(self.step_size, err), MAX_STEP)
            if self.state[3] > 0:
                self.step_size = min(self.step_size, controller_step)

            if self.state[1] < 0:
                self.done, self.end_reason = True, "touchdown"
            elif abs(self.state[4]) > 1 and self.t < burnout_time:
                self.done, self.end_reason = True, "over-rotation"
            elif self.max_time is not None and self.t >= self.max_time - 1e-9:
                self.done, self.end_reason = True, "runtime reached"


class RocketSprite:
    """Builds a frame of the rocket from rocket.png + flame / air brakes / parachute overlays."""

    CW, CH = 240, 500            # canvas size (px)
    RX, RY = 88, 170             # where the rocket png's top-left goes on the canvas

    def __init__(self, png_path):
        if not Path(png_path).exists():
            from make_rocket_png import make_rocket
            make_rocket(str(png_path))
        self.rocket = Image.open(png_path).convert("RGBA")
        rw, rh = self.rocket.size
        self.cx = self.RX + rw / 2                 # rocket centerline
        self.nose_y = self.RY + 2
        self.body_half = rw * 14 / 64              # body tube half-width (matches make_rocket_png)
        self.tail_y = self.RY + rh * 222 / 240     # nozzle exit
        self.brake_y = self.RY + rh * 150 / 240    # air brake hinge line
        self.rng = np.random.default_rng(0)

    def frame(self, thrust, control_value, pitch_deg, descending, main_out):
        img = Image.new("RGBA", (self.CW, self.CH), (0, 0, 0, 0))
        d = ImageDraw.Draw(img)
        cx = self.cx

        # parachute(s), drawn behind the rocket
        if descending:
            r = 95 if main_out else 38
            top = self.nose_y - (150 if main_out else 110)
            box = [cx - r, top, cx + r, top + r * 1.1]
            color = (255, 140, 30, 255) if main_out else (240, 200, 40, 255)
            d.pieslice(box, 180, 360, fill=color, outline=(60, 40, 10, 255), width=2)
            for k in np.linspace(-r, r, 5):
                d.line([(cx + k, top + r * 0.55), (cx, self.nose_y)], fill=(80, 80, 80, 255), width=1)

        # motor flame
        if thrust > 1:
            frac = thrust / MAX_THRUST
            L = 25 + 110 * frac * (0.85 + 0.3 * self.rng.random())
            w = 10 + 6 * frac
            ty = self.tail_y
            d.polygon([(cx - w, ty), (cx + w, ty), (cx, ty + L)], fill=(255, 120, 0, 230))
            d.polygon([(cx - w * 0.5, ty), (cx + w * 0.5, ty), (cx, ty + L * 0.6)], fill=(255, 235, 120, 255))

        img.alpha_composite(self.rocket, (self.RX, self.RY))

        # air brake panels (hinged at the body, swing outward)
        if control_value > 0.001:
            ang = np.radians(20 * control_value * BRAKE_VISUAL_SCALE)
            plen = 32
            for side in (-1, 1):
                hx = cx + side * self.body_half
                hy = self.brake_y
                ex, ey = hx + side * plen * np.sin(ang), hy + plen * np.cos(ang)
                d.line([(hx, hy), (ex, ey)], fill=(40, 40, 50, 255), width=6)
                d.line([(hx, hy), (ex, ey)], fill=(160, 170, 190, 255), width=3)

        # PIL rotates counter-clockwise; positive pitch tilts toward +y (right), i.e. clockwise
        return np.asarray(img.rotate(-pitch_deg, resample=Image.BICUBIC, center=(cx, self.tail_y)))

    @property
    def anchor(self):
        """Fraction of the canvas (x from left, y from bottom) that sits on the rocket's (y, z) point: the nozzle."""
        return (self.cx / self.CW, 1 - self.tail_y / self.CH)


def build_figure(max_time):
    fig = plt.figure(figsize=(15, 9))
    gs = fig.add_gridspec(4, 3, width_ratios=[1.1, 1, 1], hspace=0.35, wspace=0.3)

    view = fig.add_subplot(gs[:, 0])
    view.set_title("Flight view")
    view.set_xlabel("Downrange y (m)")
    view.set_ylabel("Altitude z (m)")
    view.set_ylim(-150, target_apogee * 1.35)   # headroom so the sprite + chute aren't clipped at apogee
    view.set_xlim(-400, 400)
    view.axhspan(-150, 0, color="#6b8e4e", alpha=0.5)
    view.axhline(target_apogee, color="red", ls=":", lw=0.8)
    view.text(0.98, target_apogee + 30, f"target {target_apogee} m", color="red", fontsize=8,
              ha="right", transform=view.get_yaxis_transform())
    view.grid(alpha=0.25)

    panels = [
        ("z", "Altitude (m)"), ("v_z", "v_z (m/s)"),
        ("speed", "Speed |v| (m/s)"), ("m", "Mass (kg)"),
        ("control_value", "Air brake (0–1)"), ("predicted_apogee", "Pred. apogee (m)"),
        ("pitch", "Pitch (deg)"), ("y", "Downrange y (m)"),
    ]
    axes, lines = {}, {}
    for i, (key, label) in enumerate(panels):
        ax = fig.add_subplot(gs[i // 2, 1 + i % 2])
        ax.set_title(label, fontsize=9)
        ax.grid(alpha=0.3)
        ax.set_xlim(0, max_time if max_time else 30)
        if i // 2 == 3:
            ax.set_xlabel("Time (s)")
        lines[key], = ax.plot([], [], lw=1.4)
        axes[key] = ax
    axes["z"].axhline(target_apogee, color="red", ls=":", lw=0.8)
    axes["predicted_apogee"].axhline(target_apogee, color="red", ls=":", lw=0.8)
    axes["control_value"].set_ylim(-0.05, 1.05)
    return fig, view, axes, lines


def run(max_time=RUNTIME, speed=PLAYBACK_SPEED, clock=None, save=None):
    sim = LiveSim(max_time)
    sprite = RocketSprite(ROCKET_PNG)
    fig, view, axes, lines = build_figure(max_time)

    trail, = view.plot([], [], color="tab:blue", lw=1, alpha=0.6)
    img_box = OffsetImage(sprite.frame(0, 0, 0, False, False), zoom=SPRITE_ZOOM)
    ab = AnnotationBbox(img_box, (0, 0), frameon=False, zorder=5, box_alignment=sprite.anchor)
    view.add_artist(ab)
    hud = view.text(0.03, 0.97, "", transform=view.transAxes, va="top", family="monospace",
                    fontsize=9, bbox=dict(boxstyle="round", fc="white", alpha=0.8))

    start = {"t0": None}
    if clock is None:
        clock = wall.perf_counter

    def update(_):
        now = clock()
        if start["t0"] is None:
            start["t0"] = now
        sim.advance_to((now - start["t0"]) * speed)

        L = {k: np.asarray(v) for k, v in sim.log.items()}
        t = L["t"]
        for key, line in lines.items():
            line.set_data(t, L[key])
            ax = axes[key]
            if max_time is None and t[-1] > ax.get_xlim()[1] * 0.95:
                ax.set_xlim(0, t[-1] * 1.5)
            if key != "control_value":
                ax.relim()
                ax.autoscale_view(scalex=False)

        y, z = L["y"][-1], L["z"][-1]
        v_z, thrust, cv = L["v_z"][-1], L["thrust"][-1], L["control_value"][-1]
        descending = sim.t > burnout_time and v_z < 0
        main_out = descending and z < 300
        img_box.set_data(sprite.frame(thrust, cv, L["pitch"][-1], descending, main_out))
        ab.xybox = ab.xy = (y, z)
        trail.set_data(L["y"], L["z"])
        half = max(400, 1.3 * np.abs(L["y"]).max())
        view.set_xlim(-half, half)

        if sim.t < burnout_time:
            phase = "BOOST"
        elif not descending:
            phase = "COAST" + (" + AIR BRAKES" if cv > 0.001 else "")
        else:
            phase = "MAIN CHUTE" if main_out else "DROGUE"
        if sim.done:
            phase += f"  [{sim.end_reason}]"
        hud.set_text(f"t      {sim.t:7.2f} s\nalt    {z:7.1f} m\nspeed  {L['speed'][-1]:7.1f} m/s\n"
                     f"max    {L['z'].max():7.1f} m\n{phase}")

        if sim.done and save is None:
            anim.event_source.stop()
            print(f"Done: t = {sim.t:.3f} s ({sim.end_reason}), max altitude {L['z'].max():.3f} m "
                  f"at {t[np.argmax(L['z'])]:.3f} s")
        return []

    n_frames = None
    if save is not None:  # offline render (fake clock): enough frames to cover the run
        total = max_time if max_time else 200
        n_frames = int(total / speed * FPS) + FPS
    anim = FuncAnimation(fig, update, frames=n_frames, interval=1000 / FPS,
                         blit=False, cache_frame_data=False, repeat=False)
    if save is not None:
        anim.save(save, fps=FPS, dpi=80)
    else:
        plt.show()
    return sim


if __name__ == "__main__":
    max_time, speed = RUNTIME, PLAYBACK_SPEED
    if len(sys.argv) > 1:
        max_time = None if sys.argv[1].lower() in ("full", "none", "all") else float(sys.argv[1])
    if len(sys.argv) > 2:
        speed = float(sys.argv[2])
    run(max_time, speed)
