"""Dear ImGui parameter editor panel."""

import math
from imgui_bundle import imgui

# Vector diagram constants
_VD_SIZE = 250       # diagram height/width in pixels
_VD_RADIUS = 80      # max arrow length (pixels)
_VD_MAX_SPEED = 100  # speed (deg/s) that maps to full radius
_VD_COLORS = {
    "S1":  (66, 165, 245, 255),
    "M1":  (30, 100, 200, 255),
    "S2":  (255, 112, 67, 255),
    "M2":  (255, 167, 130, 255),
    "VA1": (76, 175, 80, 255),
    "VA2": (139, 195, 74, 255),
    "IOC1": (255, 235, 59, 255),
    "IOC2": (255, 213, 79, 255),
    "CL":  (255, 255, 255, 50),
}

_GROUP_PARAMS = [
    # (suffix, label, drag_speed, format)
    ("GratingContrast",          "Contrast",           0.01,  "%.3f"),
    ("GratingDirectionDeg",      "Direction (deg)",     1.0,   "%.1f"),
    ("GratingSpatialFreqCPD",    "Spatial Freq (CPD)",  0.001, "%.4f"),
    ("GratingTemporalFreqCPS",   "Temporal Freq (CPS)", 0.1,   "%.2f"),
    ("GratingPhaseDeg",          "Phase (deg)",         1.0,   "%.1f"),
    ("GratingDiameterDeg",       "Diameter (deg)",      1.0,   "%.1f"),
    ("GratingAzimuthDeg",        "Azimuth (deg)",       1.0,   "%.1f"),
    ("GratingElevationDeg",      "Elevation (deg)",     1.0,   "%.1f"),
]
_GROUP_PREFIXES = ["stimOne", "maskOne", "stimTwo", "maskTwo"]


class ImguiPanel:
    """Renders the parameter editor sidebar."""

    def __init__(self, model):
        self.model = model
        self.search_text = ""
        self.show_modified_only = False
        self.save_flash_timer = 0.0  # countdown for "Saved!" flash
        self.panel_width = 420       # default sidebar width
        self.show_stim_two = False
        self.show_mask_two = False
        self.overlay_enabled = False
        self.overlay_opacity = 0.5

        # Group Controls state
        self._group_linked = {"stimOne": True, "maskOne": True, "stimTwo": False, "maskTwo": False}
        self._group_deltas = {suffix: 0.0 for suffix, *_ in _GROUP_PARAMS}
        self._group_param_enabled = {suffix: True for suffix, *_ in _GROUP_PARAMS}

    def render(self, dt: float):
        """Render the full imgui panel. Called each frame."""
        self.save_flash_timer = max(0.0, self.save_flash_timer - dt)

        io = imgui.get_io()
        display_w = io.display_size.x
        display_h = io.display_size.y

        # Position panel on the right side
        imgui.set_next_window_pos(imgui.ImVec2(display_w - self.panel_width, 0))
        imgui.set_next_window_size(imgui.ImVec2(self.panel_width, display_h))

        flags = (imgui.WindowFlags_.no_move | imgui.WindowFlags_.no_resize |
                 imgui.WindowFlags_.no_collapse)

        expanded, _ = imgui.begin("Parameters", True, flags)
        if expanded:
            imgui.begin_child("##scroll", imgui.ImVec2(0, 0))
            self._render_top_bar()
            imgui.separator()
            self._render_grating_params()
            imgui.separator()
            self._render_group_controls()
            imgui.separator()
            self._render_vector_diagram()
            imgui.separator()
            self._render_all_variables()
            imgui.end_child()
        imgui.end()

        # Keyboard shortcut: Ctrl+S to save
        if io.key_ctrl and imgui.is_key_pressed(imgui.Key.s):
            self._do_save()

    def _do_save(self):
        self.model.save()
        self.save_flash_timer = 2.0

    def _render_top_bar(self):
        """Save button, search, modified count."""
        if imgui.button("Save to XML"):
            self._do_save()
        imgui.same_line()
        if imgui.button("Reset All"):
            self.model.reset_all()

        if self.save_flash_timer > 0:
            imgui.same_line()
            imgui.text_colored(imgui.ImVec4(0.2, 1.0, 0.2, 1.0), "Saved!")

        imgui.same_line()
        mod_count = self.model.get_modified_count()
        if mod_count > 0:
            imgui.text_colored(imgui.ImVec4(1.0, 1.0, 0.0, 1.0), f"({mod_count} modified)")

        # Search box
        imgui.set_next_item_width(-1)
        changed, self.search_text = imgui.input_text_with_hint(
            "##search", "Search variables...", self.search_text
        )

        _, self.show_modified_only = imgui.checkbox("Show modified only", self.show_modified_only)
        _, self.show_stim_two = imgui.checkbox("Show Stim 2", self.show_stim_two)
        imgui.same_line()
        _, self.show_mask_two = imgui.checkbox("Show Mask 2", self.show_mask_two)

    def _render_grating_params(self):
        """Specialized sliders for grating preview parameters."""
        if not imgui.collapsing_header("Grating Preview Parameters", imgui.TreeNodeFlags_.default_open):
            return

        self._render_stim_section("Stimulus One", "stimOne")
        self._render_mask_section("Mask One", "maskOne")
        self._render_stim_section("Stimulus Two", "stimTwo")
        self._render_mask_section("Mask Two", "maskTwo")

    def _render_stim_section(self, label: str, prefix: str):
        """Render sliders for a stimulus (stimOne or stimTwo)."""
        if not imgui.tree_node(label):
            return

        tag_contrast = f"{prefix}GratingContrast"
        tag_direction = f"{prefix}GratingDirectionDeg"
        tag_sf = f"{prefix}GratingSpatialFreqCPD"
        tag_tf = f"{prefix}GratingTemporalFreqCPS"
        tag_phase = f"{prefix}GratingPhaseDeg"
        tag_diameter = f"{prefix}GratingDiameterDeg"
        tag_azimuth = f"{prefix}GratingAzimuthDeg"
        tag_elevation = f"{prefix}GratingElevationDeg"

        self._slider_float(tag_contrast, "Contrast", 0.0, 1.0)
        self._slider_float(tag_direction, "Direction (deg)", 0.0, 360.0)
        self._slider_float(tag_sf, "Spatial Freq (CPD)", 0.001, 1.0)
        self._slider_float(tag_tf, "Temporal Freq (CPS)", 0.0, 5.0)
        self._slider_float(tag_phase, "Phase (deg)", 0.0, 360.0)
        self._drag_float(tag_diameter, "Diameter (deg)", 1.0, 400.0)
        self._slider_float(tag_azimuth, "Azimuth (deg)", -100.0, 100.0)
        self._slider_float(tag_elevation, "Elevation (deg)", -100.0, 100.0)

        imgui.tree_pop()

    def _render_mask_section(self, label: str, prefix: str):
        """Render sliders for a mask — all 8 parameters (independent from parent)."""
        if not imgui.tree_node(label):
            return

        tag_contrast = f"{prefix}GratingContrast"
        tag_direction = f"{prefix}GratingDirectionDeg"
        tag_sf = f"{prefix}GratingSpatialFreqCPD"
        tag_tf = f"{prefix}GratingTemporalFreqCPS"
        tag_phase = f"{prefix}GratingPhaseDeg"
        tag_diameter = f"{prefix}GratingDiameterDeg"
        tag_azimuth = f"{prefix}GratingAzimuthDeg"
        tag_elevation = f"{prefix}GratingElevationDeg"

        self._slider_float(tag_contrast, "Contrast", 0.0, 1.0)
        self._slider_float(tag_direction, "Direction (deg)", 0.0, 360.0)
        self._slider_float(tag_sf, "Spatial Freq (CPD)", 0.001, 1.0)
        self._slider_float(tag_tf, "Temporal Freq (CPS)", 0.0, 5.0)
        self._slider_float(tag_phase, "Phase (deg)", 0.0, 360.0)
        self._drag_float(tag_diameter, "Diameter (deg)", 1.0, 400.0)
        self._slider_float(tag_azimuth, "Azimuth (deg)", -100.0, 100.0)
        self._slider_float(tag_elevation, "Elevation (deg)", -100.0, 100.0)

        imgui.tree_pop()

    def _slider_float(self, tag: str, label: str, v_min: float, v_max: float):
        """Render a slider + input field for a variable, with modified highlighting."""
        if tag not in self.model.variables:
            return
        mv = self.model.variables[tag]

        was_modified = mv.is_modified
        if was_modified:
            imgui.push_style_color(imgui.Col_.frame_bg, imgui.ImVec4(0.4, 0.4, 0.0, 0.5))

        val = float(mv.current_value)

        # Label text first, then measure remaining width for slider+input
        imgui.text(label)
        imgui.same_line()
        remaining = imgui.get_content_region_avail().x

        # Slider (hidden label via ## prefix)
        imgui.set_next_item_width(remaining * 0.65)
        changed1, new_val1 = imgui.slider_float(f"##{tag}", val, v_min, v_max)
        if changed1:
            self.model.set(tag, new_val1)
            val = new_val1

        # Right-click slider to reset
        if imgui.is_item_clicked(imgui.MouseButton_.right):
            self.model.reset(tag)

        # Input box next to slider
        imgui.same_line()
        imgui.set_next_item_width(remaining * 0.35)
        changed2, new_val2 = imgui.input_float(f"##input_{tag}", val, 0.0, 0.0, "%.4f")
        if changed2:
            self.model.set(tag, new_val2)

        if was_modified:
            imgui.pop_style_color()

    def _drag_float(self, tag: str, label: str, v_min: float, v_max: float):
        """Render a drag float + input field for a variable."""
        if tag not in self.model.variables:
            return
        mv = self.model.variables[tag]

        was_modified = mv.is_modified
        if was_modified:
            imgui.push_style_color(imgui.Col_.frame_bg, imgui.ImVec4(0.4, 0.4, 0.0, 0.5))

        val = float(mv.current_value)

        # Label text first, then measure remaining width for drag+input
        imgui.text(label)
        imgui.same_line()
        remaining = imgui.get_content_region_avail().x

        imgui.set_next_item_width(remaining * 0.65)
        changed1, new_val1 = imgui.drag_float(f"##{tag}", val, 1.0, v_min, v_max)
        if changed1:
            self.model.set(tag, new_val1)
            val = new_val1

        if imgui.is_item_clicked(imgui.MouseButton_.right):
            self.model.reset(tag)

        imgui.same_line()
        imgui.set_next_item_width(remaining * 0.35)
        changed2, new_val2 = imgui.input_float(f"##input_{tag}", val, 0.0, 0.0, "%.1f")
        if changed2:
            self.model.set(tag, new_val2)

        if was_modified:
            imgui.pop_style_color()

    def _apply_group_delta(self, suffix: str, diff: float):
        """Apply an incremental delta to all linked gratings for the given parameter suffix."""
        for prefix in _GROUP_PREFIXES:
            if not self._group_linked.get(prefix, False):
                continue
            tag = f"{prefix}{suffix}"
            if tag not in self.model.variables:
                continue
            cur = float(self.model.get(tag))
            self.model.set(tag, cur + diff)

    def _render_group_controls(self):
        """Collapsible section for adjusting multiple gratings simultaneously."""
        if not imgui.collapsing_header("Group Controls"):
            return

        # Linkage checkboxes
        imgui.text("Linked:")
        imgui.same_line()
        for prefix, label in [("stimOne", "Stim 1"), ("maskOne", "Mask 1"),
                               ("stimTwo", "Stim 2"), ("maskTwo", "Mask 2")]:
            _, self._group_linked[prefix] = imgui.checkbox(
                f"{label}##grp_link", self._group_linked[prefix])
            imgui.same_line()
        imgui.new_line()

        # Zero All button — reverses accumulated deltas
        if imgui.button("Zero All Deltas"):
            for suffix, *_ in _GROUP_PARAMS:
                old = self._group_deltas[suffix]
                if old != 0.0:
                    self._apply_group_delta(suffix, -old)
                    self._group_deltas[suffix] = 0.0

        imgui.separator()

        # Parameter rows
        for suffix, label, speed, fmt in _GROUP_PARAMS:
            enabled = self._group_param_enabled[suffix]
            _, self._group_param_enabled[suffix] = imgui.checkbox(
                f"##{suffix}_enabled", enabled)
            imgui.same_line()

            if not enabled:
                imgui.begin_disabled()

            imgui.text(label)
            imgui.same_line()
            remaining = imgui.get_content_region_avail().x - 50  # reserve space for Zero btn

            imgui.set_next_item_width(max(remaining, 60))
            old_delta = self._group_deltas[suffix]
            changed, new_delta = imgui.drag_float(
                f"##grp_{suffix}", old_delta, speed, 0.0, 0.0, fmt)
            if changed and enabled:
                diff = new_delta - old_delta
                self._apply_group_delta(suffix, diff)
                self._group_deltas[suffix] = new_delta

            if not enabled:
                imgui.end_disabled()

            imgui.same_line()
            if imgui.small_button(f"Zero##{suffix}"):
                old = self._group_deltas[suffix]
                if old != 0.0:
                    self._apply_group_delta(suffix, -old)
                    self._group_deltas[suffix] = 0.0

    # ------------------------------------------------------------------
    #  Vector diagram
    # ------------------------------------------------------------------

    def _get_grating_vector_data(self, prefix: str):
        """Return (direction_deg, tf_cps, sf_cpd) for a grating/mask prefix.

        Returns None if the direction variable is missing.
        """
        dir_tag = f"{prefix}GratingDirectionDeg"
        if dir_tag not in self.model.variables:
            return None
        direction = float(self.model.get(dir_tag))

        tf_tag = f"{prefix}GratingTemporalFreqCPS"
        tf = float(self.model.get(tf_tag)) if tf_tag in self.model.variables else 0.0

        sf_tag = f"{prefix}GratingSpatialFreqCPD"
        sf = float(self.model.get(sf_tag)) if sf_tag in self.model.variables else 0.0

        return direction, tf, sf

    def _draw_arrow(self, dl, cx: float, cy: float,
                    direction_deg: float, speed: float,
                    color_tuple: tuple, label: str,
                    radius: float = _VD_RADIUS, max_speed: float = _VD_MAX_SPEED,
                    clamp: bool = True):
        """Draw a single direction/speed arrow from centre."""
        col = imgui.IM_COL32(*color_tuple)

        if speed <= 0:
            # Dot at centre
            dl.add_circle_filled(imgui.ImVec2(cx, cy), 4 * radius / _VD_RADIUS, col)
            return

        length = (speed / max_speed) * radius
        length = max(length, 8 * radius / _VD_RADIUS)   # minimum visible length
        if clamp:
            length = min(length, radius)

        # Angle: 0 deg = right, 90 deg = up.  Negate for screen Y-down.
        angle = math.radians(-direction_deg)
        ex = cx + length * math.cos(angle)
        ey = cy + length * math.sin(angle)

        # Shaft
        scale = radius / _VD_RADIUS
        dl.add_line(imgui.ImVec2(cx, cy), imgui.ImVec2(ex, ey), col, 2.0 * scale)

        # Arrowhead (filled triangle)
        head_len = 8.0 * scale
        head_half = 4.0 * scale
        dx = ex - cx
        dy = ey - cy
        norm = math.hypot(dx, dy)
        if norm < 1e-6:
            return
        ux, uy = dx / norm, dy / norm  # unit along shaft
        px, py = -uy, ux               # perpendicular

        p1 = imgui.ImVec2(ex, ey)
        p2 = imgui.ImVec2(ex - ux * head_len + px * head_half,
                          ey - uy * head_len + py * head_half)
        p3 = imgui.ImVec2(ex - ux * head_len - px * head_half,
                          ey - uy * head_len - py * head_half)
        dl.add_triangle_filled(p1, p2, p3, col)

    def _draw_constraint_line(self, dl, cx: float, cy: float,
                              direction_deg: float, speed: float,
                              radius: float = _VD_RADIUS, max_speed: float = _VD_MAX_SPEED,
                              color_tuple: tuple = _VD_COLORS["CL"]):
        """Draw thin constraint line through arrow tip, perpendicular to direction."""
        if speed <= 0:
            return

        col = imgui.IM_COL32(*color_tuple)

        length = (speed / max_speed) * radius

        angle = math.radians(-direction_deg)
        # Arrow tip
        tx = cx + length * math.cos(angle)
        ty = cy + length * math.sin(angle)

        # Perpendicular direction (along grating orientation)
        px = -math.sin(angle)
        py = math.cos(angle)

        ext = radius + 15 * radius / _VD_RADIUS
        dl.add_line(imgui.ImVec2(tx - px * ext, ty - py * ext),
                    imgui.ImVec2(tx + px * ext, ty + py * ext), col, 1.0)

    @staticmethod
    def _signed_angle_diff(from_deg: float, to_deg: float) -> float:
        """Smallest signed angle from from_deg to to_deg, in [-180, 180]."""
        d = (to_deg - from_deg) % 360
        if d > 180:
            d -= 360
        return d

    @staticmethod
    def _compute_va(dir1: float, speed1: float,
                    dir2: float, speed2: float):
        """Compute Vector Average of two velocity vectors.

        Returns (direction_deg, speed) or None if magnitude < epsilon.
        """
        a1 = math.radians(dir1)
        a2 = math.radians(dir2)
        vx = (speed1 * math.cos(a1) + speed2 * math.cos(a2)) / 2.0
        vy = (speed1 * math.sin(a1) + speed2 * math.sin(a2)) / 2.0
        mag = math.hypot(vx, vy)
        if mag < 1e-6:
            return None
        d = math.degrees(math.atan2(vy, vx))
        return d % 360, mag

    @staticmethod
    def _compute_ioc(dir1: float, speed1: float,
                     dir2: float, speed2: float):
        """Compute Intersection of Constraints of two velocity vectors.

        Each grating defines a constraint line in velocity space:
          Point: P = speed*(cos(theta), sin(theta))
          Direction along line: (-sin(theta), cos(theta))
        Solve for intersection via Cramer's rule.
        Returns (direction_deg, speed) or None when det ~ 0 (parallel).
        """
        a1 = math.radians(dir1)
        a2 = math.radians(dir2)

        # Points on constraint lines (velocity tips)
        p1x = speed1 * math.cos(a1)
        p1y = speed1 * math.sin(a1)
        p2x = speed2 * math.cos(a2)
        p2y = speed2 * math.sin(a2)

        # Directions along constraint lines (perpendicular to velocity = along grating orientation)
        d1x = -math.sin(a1)
        d1y = math.cos(a1)
        d2x = -math.sin(a2)
        d2y = math.cos(a2)

        # Solve: P1 + t*D1 = P2 + s*D2
        # => t*D1 - s*D2 = P2 - P1
        # | d1x  -d2x | |t|   |p2x - p1x|
        # | d1y  -d2y | |s| = |p2y - p1y|
        det = d1x * (-d2y) - (-d2x) * d1y
        if abs(det) < 1e-6:
            return None

        dpx = p2x - p1x
        dpy = p2y - p1y

        t = (dpx * (-d2y) - (-d2x) * dpy) / det

        ix = p1x + t * d1x
        iy = p1y + t * d1y
        mag = math.hypot(ix, iy)
        if mag < 1e-6:
            return None
        d = math.degrees(math.atan2(iy, ix))
        return d % 360, mag

    def _draw_vector_diagram_core(self, dl, cx: float, cy: float,
                                   radius: float, max_speed: float,
                                   alpha_mult: float = 1.0, draw_legend: bool = True,
                                   clip_rect: tuple = None):
        """Draw the full vector diagram into an arbitrary draw list at (cx, cy).

        clip_rect: optional (x1, y1, x2, y2) to clip arrows/lines but not labels.
        """
        scale = radius / _VD_RADIUS

        def _a(color):
            """Scale alpha channel by alpha_mult."""
            return (*color[:3], int(color[3] * alpha_mult))

        faint = imgui.IM_COL32(*_a((255, 255, 255, 40)))
        label_col = imgui.IM_COL32(*_a((255, 255, 255, 100)))

        # Reference circles
        dl.add_circle(imgui.ImVec2(cx, cy), radius, faint, 64, 1.0)
        dl.add_circle(imgui.ImVec2(cx, cy), radius * 0.5, faint, 64, 1.0)

        # Crosshair axes
        ext_ax = radius + 10 * scale
        dl.add_line(imgui.ImVec2(cx - ext_ax, cy),
                    imgui.ImVec2(cx + ext_ax, cy), faint, 1.0)
        dl.add_line(imgui.ImVec2(cx, cy - ext_ax),
                    imgui.ImVec2(cx, cy + ext_ax), faint, 1.0)

        # Cardinal labels (0=right, 90=up, 180=left, 270=down)
        off = radius + 14 * scale
        dl.add_text(imgui.ImVec2(cx + off, cy - 6), label_col, "0")
        dl.add_text(imgui.ImVec2(cx - 8, cy - off - 12), label_col, "90")
        dl.add_text(imgui.ImVec2(cx - off - 20, cy - 6), label_col, "180")
        dl.add_text(imgui.ImVec2(cx - 10, cy + off), label_col, "270")

        # Clip arrows/constraint lines so they don't escape the diagram
        if clip_rect:
            dl.push_clip_rect(imgui.ImVec2(clip_rect[0], clip_rect[1]),
                              imgui.ImVec2(clip_rect[2], clip_rect[3]), True)

        # Determine which arrows to draw
        entries: list[tuple[str, str, tuple]] = []
        entries.append(("S1", "stimOne", _VD_COLORS["S1"]))
        entries.append(("M1", "maskOne", _VD_COLORS["M1"]))
        if self.show_stim_two:
            entries.append(("S2", "stimTwo", _VD_COLORS["S2"]))
        if self.show_mask_two:
            entries.append(("M2", "maskTwo", _VD_COLORS["M2"]))

        # Collect direction/speed data for all entries
        arrow_data: dict[str, tuple[float, float]] = {}
        for label, prefix, color in entries:
            data = self._get_grating_vector_data(prefix)
            if data is None:
                continue
            direction, tf, sf = data
            speed = tf / sf if sf > 1e-6 else 0.0
            arrow_data[label] = (direction, speed)
            self._draw_arrow(dl, cx, cy, direction, speed, _a(color), label,
                             radius=radius, max_speed=max_speed)

        # VA/IOC for Pair 1 (S1 + M1)
        va_ioc_legend: list[tuple[str, tuple]] = []
        if "S1" in arrow_data and "M1" in arrow_data:
            d1, sp1 = arrow_data["S1"]
            d2, sp2 = arrow_data["M1"]
            if sp1 > 0 and sp2 > 0:
                self._draw_constraint_line(dl, cx, cy, d1, sp1, radius=radius, max_speed=max_speed, color_tuple=_a(_VD_COLORS["CL"]))
                self._draw_constraint_line(dl, cx, cy, d2, sp2, radius=radius, max_speed=max_speed, color_tuple=_a(_VD_COLORS["CL"]))
                va = self._compute_va(d1, sp1, d2, sp2)
                if va:
                    self._draw_arrow(dl, cx, cy, va[0], va[1], _a(_VD_COLORS["VA1"]), "VA1",
                                     radius=radius, max_speed=max_speed)
                    va_ioc_legend.append(("VA1", _VD_COLORS["VA1"]))
                ioc = self._compute_ioc(d1, sp1, d2, sp2)
                if ioc:
                    self._draw_arrow(dl, cx, cy, ioc[0], ioc[1], _a(_VD_COLORS["IOC1"]), "IOC1",
                                     radius=radius, max_speed=max_speed, clamp=False)
                    va_ioc_legend.append(("IOC1", _VD_COLORS["IOC1"]))

        # VA/IOC for Pair 2 (S2 + M2)
        if self.show_stim_two and self.show_mask_two:
            if "S2" in arrow_data and "M2" in arrow_data:
                d1, sp1 = arrow_data["S2"]
                d2, sp2 = arrow_data["M2"]
                if sp1 > 0 and sp2 > 0:
                    self._draw_constraint_line(dl, cx, cy, d1, sp1, radius=radius, max_speed=max_speed)
                    self._draw_constraint_line(dl, cx, cy, d2, sp2, radius=radius, max_speed=max_speed)
                    va = self._compute_va(d1, sp1, d2, sp2)
                    if va:
                        self._draw_arrow(dl, cx, cy, va[0], va[1], _a(_VD_COLORS["VA2"]), "VA2",
                                         radius=radius, max_speed=max_speed)
                        va_ioc_legend.append(("VA2", _VD_COLORS["VA2"]))
                    ioc = self._compute_ioc(d1, sp1, d2, sp2)
                    if ioc:
                        self._draw_arrow(dl, cx, cy, ioc[0], ioc[1], _a(_VD_COLORS["IOC2"]), "IOC2",
                                         radius=radius, max_speed=max_speed, clamp=False)
                        va_ioc_legend.append(("IOC2", _VD_COLORS["IOC2"]))

        # Pop arrow clip before drawing legend text
        if clip_rect:
            dl.pop_clip_rect()

        # Legend
        if draw_legend:
            legend_y = cy + radius + 28 * scale
            lx = cx - radius
            sq = 10 * scale
            for label, _, color in entries:
                col = imgui.IM_COL32(*_a(color))
                dl.add_rect_filled(imgui.ImVec2(lx, legend_y),
                                   imgui.ImVec2(lx + sq, legend_y + sq), col)
                dl.add_text(imgui.ImVec2(lx + sq + 4, legend_y - 2), label_col, label)
                lx += 50 * scale

            dl.add_text(imgui.ImVec2(lx + 10 * scale, legend_y - 2), label_col,
                        f"Ring = {_VD_MAX_SPEED:.0f} deg/s")

            if va_ioc_legend:
                legend_y2 = legend_y + 16 * scale
                lx2 = cx - radius
                for label, color in va_ioc_legend:
                    col = imgui.IM_COL32(*_a(color))
                    dl.add_rect_filled(imgui.ImVec2(lx2, legend_y2),
                                       imgui.ImVec2(lx2 + sq, legend_y2 + sq), col)
                    dl.add_text(imgui.ImVec2(lx2 + sq + 4, legend_y2 - 2), label_col, label)
                    lx2 += 60 * scale

    def _render_vector_diagram(self):
        """Collapsible section showing direction/speed vector arrows with VA/IOC."""
        if not imgui.collapsing_header("Direction / Speed Vectors",
                                       imgui.TreeNodeFlags_.default_open):
            return

        avail_w = imgui.get_content_region_avail().x
        cursor = imgui.get_cursor_screen_pos()
        cx = cursor.x + avail_w * 0.5
        cy = cursor.y + _VD_SIZE * 0.4
        imgui.dummy(imgui.ImVec2(avail_w, _VD_SIZE))
        dl = imgui.get_window_draw_list()

        clip = (cursor.x, cursor.y, cursor.x + avail_w, cursor.y + _VD_SIZE)
        self._draw_vector_diagram_core(dl, cx, cy, _VD_RADIUS, _VD_MAX_SPEED, 1.0, True,
                                       clip_rect=clip)

        # Angle readouts (M1 → VA, M1 → IOC)
        s1_data = self._get_grating_vector_data("stimOne")
        m1_data = self._get_grating_vector_data("maskOne")
        if s1_data and m1_data:
            d1, tf1, sf1 = s1_data
            d2, tf2, sf2 = m1_data
            sp1 = tf1 / sf1 if sf1 > 1e-6 else 0.0
            sp2 = tf2 / sf2 if sf2 > 1e-6 else 0.0
            if sp1 > 0 and sp2 > 0:
                va = self._compute_va(d1, sp1, d2, sp2)
                ioc = self._compute_ioc(d1, sp1, d2, sp2)
                if va:
                    angle_va = self._signed_angle_diff(d2, va[0])
                    imgui.text(f"M1 -> VA:  {angle_va:+.1f} deg")
                if ioc:
                    angle_ioc = self._signed_angle_diff(d2, ioc[0])
                    imgui.text(f"M1 -> IOC: {angle_ioc:+.1f} deg")

        imgui.dummy(imgui.ImVec2(avail_w, 4))

        # Visual overlay controls
        imgui.text("Viewport Overlay")
        imgui.same_line()
        _, self.overlay_enabled = imgui.checkbox("Enable##overlay", self.overlay_enabled)
        if self.overlay_enabled:
            imgui.same_line()
            imgui.set_next_item_width(imgui.get_content_region_avail().x)
            _, self.overlay_opacity = imgui.slider_float(
                "##overlay_opacity", self.overlay_opacity, 0.05, 1.0, "%.2f")

    def render_overlay(self, viewport_w: float, viewport_h: float):
        """Draw the vector diagram as a semi-transparent overlay on the grating viewport."""
        if not self.overlay_enabled:
            return
        shorter = min(viewport_w, viewport_h)
        overlay_radius = shorter * 0.22
        cx = viewport_w * 0.5
        cy = viewport_h * 0.5

        # Offset overlay to track stimOne azimuth/elevation (deg -> pixels)
        az_tag = "stimOneGratingAzimuthDeg"
        el_tag = "stimOneGratingElevationDeg"
        if az_tag in self.model.variables:
            cx += (float(self.model.get(az_tag)) / 100.0) * (viewport_w * 0.5)
        if el_tag in self.model.variables:
            cy -= (float(self.model.get(el_tag)) / 100.0) * (viewport_h * 0.5)

        dl = imgui.get_foreground_draw_list()
        dl.push_clip_rect(imgui.ImVec2(0, 0),
                          imgui.ImVec2(viewport_w, viewport_h), True)
        self._draw_vector_diagram_core(dl, cx, cy, overlay_radius, _VD_MAX_SPEED,
                                       self.overlay_opacity, True)
        dl.pop_clip_rect()

    def _render_all_variables(self):
        """Tree of ALL variables organized by group, with auto-selected widgets."""
        header_open = imgui.collapsing_header("All Variables")
        if not header_open:
            return

        search = self.search_text.lower()

        # Render grouped variables
        for group_name, tags in self.model.groups.items():
            filtered = [t for t in tags if self._passes_filter(t, search)]
            if not filtered:
                continue

            if imgui.tree_node(f"{group_name} ({len(filtered)})"):
                for tag in filtered:
                    self._render_variable_widget(tag)
                imgui.tree_pop()

        # Render ungrouped variables by folder
        if imgui.tree_node("Ungrouped"):
            for folder_name, tags in self.model.folder_vars.items():
                ungrouped = [t for t in tags
                            if self.model.variables[t].groups == "" and self._passes_filter(t, search)]
                if not ungrouped:
                    continue
                if imgui.tree_node(f"{folder_name} ({len(ungrouped)})"):
                    for tag in ungrouped:
                        self._render_variable_widget(tag)
                    imgui.tree_pop()
            imgui.tree_pop()

    def _passes_filter(self, tag: str, search: str) -> bool:
        """Check if a variable passes the current search/modified filter."""
        mv = self.model.variables[tag]
        if self.show_modified_only and not mv.is_modified:
            return False
        if search and search not in tag.lower() and search not in mv.groups.lower():
            return False
        return True

    def _render_variable_widget(self, tag: str):
        """Render the appropriate widget for a variable based on its type."""
        mv = self.model.variables[tag]

        was_modified = mv.is_modified
        if was_modified:
            imgui.push_style_color(imgui.Col_.text, imgui.ImVec4(1.0, 1.0, 0.0, 1.0))

        vt = mv.var_type.lower()

        if vt == "boolean":
            val = bool(mv.current_value)
            changed, new_val = imgui.checkbox(f"{tag}", val)
            if changed:
                self.model.set(tag, new_val)

        elif vt == "float":
            val = float(mv.current_value)
            imgui.set_next_item_width(150)
            changed, new_val = imgui.input_float(f"{tag}", val, 0.0, 0.0, "%.4f")
            if changed:
                self.model.set(tag, new_val)

        elif vt == "integer":
            val = int(mv.current_value)
            imgui.set_next_item_width(150)
            changed, new_val = imgui.input_int(f"{tag}", val)
            if changed:
                self.model.set(tag, new_val)

        else:  # string, list, selection
            val = str(mv.current_value)
            imgui.set_next_item_width(200)
            changed, new_val = imgui.input_text(f"{tag}", val)
            if changed:
                self.model.set(tag, new_val)

        if was_modified:
            imgui.pop_style_color()

        # Right-click context menu
        if imgui.is_item_clicked(imgui.MouseButton_.right):
            self.model.reset(tag)
