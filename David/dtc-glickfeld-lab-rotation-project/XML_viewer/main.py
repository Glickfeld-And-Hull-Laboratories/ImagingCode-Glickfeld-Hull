"""MonkeyML Grating Viewer — main entry point."""

import os
import sys
import time

from imgui_bundle import imgui, hello_imgui
import moderngl

from xml_model import ExperimentModel
from grating_renderer import GratingRenderer, GratingInstance
from imgui_panel import ImguiPanel


class App:
    def __init__(self):
        # Find XML file
        self.base_dir = os.path.dirname(os.path.abspath(__file__))
        xml_path = os.path.join(self.base_dir, 'twoStim_2P_Frames.xml')

        # Load experiment model
        self.model = ExperimentModel(xml_path)

        # Grating instances
        self.stim_one = GratingInstance("stimOne")
        self.mask_one = GratingInstance("maskOne")
        self.stim_two = GratingInstance("stimTwo")
        self.mask_two = GratingInstance("maskTwo")
        self.gratings = [self.stim_one, self.mask_one, self.stim_two, self.mask_two]

        # Will be initialized after OpenGL context exists
        self.renderer = None
        self.ctx = None
        self.panel = ImguiPanel(self.model)

        # Timing
        self.last_time = time.perf_counter()
        self.playing = True
        self.fps_counter = 0
        self.fps_time = 0.0
        self.fps_display = 0

    def post_init(self):
        """Called after hello_imgui creates the OpenGL context."""
        self.ctx = moderngl.create_context(standalone=False)
        self.renderer = GratingRenderer(self.ctx)

    def sync_model_to_gratings(self):
        """Copy current XML parameter values into grating instances."""
        m = self.model

        # Stimulus One
        self.stim_one.contrast = float(m.get('stimOneGratingContrast') or 1.0)
        self.stim_one.direction_deg = float(m.get('stimOneGratingDirectionDeg') or 0.0)
        self.stim_one.spatial_freq = float(m.get('stimOneGratingSpatialFreqCPD') or 0.05)
        self.stim_one.temporal_freq = float(m.get('stimOneGratingTemporalFreqCPS') or 0.0)
        self.stim_one.phase_deg = float(m.get('stimOneGratingPhaseDeg') or 0.0)
        self.stim_one.diameter_deg = float(m.get('stimOneGratingDiameterDeg') or 200.0)
        self.stim_one.azimuth_deg = float(m.get('stimOneGratingAzimuthDeg') or 0.0)
        self.stim_one.elevation_deg = float(m.get('stimOneGratingElevationDeg') or 0.0)

        # Mask One: fully independent parameters
        self.mask_one.contrast = float(m.get('maskOneGratingContrast') or 0.0)
        self.mask_one.direction_deg = float(m.get('maskOneGratingDirectionDeg') or 0.0)
        self.mask_one.phase_deg = float(m.get('maskOneGratingPhaseDeg') or 0.0)
        self.mask_one.spatial_freq = float(m.get('maskOneGratingSpatialFreqCPD') or 0.05)
        self.mask_one.temporal_freq = float(m.get('maskOneGratingTemporalFreqCPS') or 0.0)
        self.mask_one.diameter_deg = float(m.get('maskOneGratingDiameterDeg') or 200.0)
        self.mask_one.azimuth_deg = float(m.get('maskOneGratingAzimuthDeg') or 0.0)
        self.mask_one.elevation_deg = float(m.get('maskOneGratingElevationDeg') or 0.0)

        # Stimulus Two
        self.stim_two.contrast = float(m.get('stimTwoGratingContrast') or 1.0)
        self.stim_two.direction_deg = float(m.get('stimTwoGratingDirectionDeg') or 0.0)
        self.stim_two.spatial_freq = float(m.get('stimTwoGratingSpatialFreqCPD') or 0.05)
        self.stim_two.temporal_freq = float(m.get('stimTwoGratingTemporalFreqCPS') or 0.0)
        self.stim_two.phase_deg = float(m.get('stimTwoGratingPhaseDeg') or 0.0)
        self.stim_two.diameter_deg = float(m.get('stimTwoGratingDiameterDeg') or 200.0)
        self.stim_two.azimuth_deg = float(m.get('stimTwoGratingAzimuthDeg') or 0.0)
        self.stim_two.elevation_deg = float(m.get('stimTwoGratingElevationDeg') or 0.0)

        # Mask Two: fully independent parameters
        self.mask_two.contrast = float(m.get('maskTwoGratingContrast') or 0.0)
        self.mask_two.direction_deg = float(m.get('maskTwoGratingDirectionDeg') or 0.0)
        self.mask_two.phase_deg = float(m.get('maskTwoGratingPhaseDeg') or 0.0)
        self.mask_two.spatial_freq = float(m.get('maskTwoGratingSpatialFreqCPD') or 0.05)
        self.mask_two.temporal_freq = float(m.get('maskTwoGratingTemporalFreqCPS') or 0.0)
        self.mask_two.diameter_deg = float(m.get('maskTwoGratingDiameterDeg') or 200.0)
        self.mask_two.azimuth_deg = float(m.get('maskTwoGratingAzimuthDeg') or 0.0)
        self.mask_two.elevation_deg = float(m.get('maskTwoGratingElevationDeg') or 0.0)

    def render_background(self):
        """custom_background callback: OpenGL grating rendering."""
        if self.renderer is None:
            return

        now = time.perf_counter()
        dt = now - self.last_time
        self.last_time = now

        # FPS counter
        self.fps_counter += 1
        self.fps_time += dt
        if self.fps_time >= 1.0:
            self.fps_display = self.fps_counter
            self.fps_counter = 0
            self.fps_time = 0.0

        # Sync parameters
        self.sync_model_to_gratings()

        # Update phase animation
        if self.playing:
            for g in self.gratings:
                g.update_phase(dt)

        # Get viewport size (leave room for imgui panel on the right)
        io = imgui.get_io()
        display_w = int(io.display_size.x)
        display_h = int(io.display_size.y)

        # Retina scale: framebuffer pixels vs logical pixels
        fb_scale_x = io.display_framebuffer_scale.x if io.display_framebuffer_scale.x > 0 else 1.0
        fb_scale_y = io.display_framebuffer_scale.y if io.display_framebuffer_scale.y > 0 else 1.0

        fb_w = int(display_w * fb_scale_x)
        fb_h = int(display_h * fb_scale_y)
        panel_fb_w = int(self.panel.panel_width * fb_scale_x)

        vp_w = max(1, fb_w - panel_fb_w)
        vp_h = max(1, fb_h)

        # All gratings render in the same centered viewport (overlaid)
        self.ctx.viewport = (0, 0, vp_w, vp_h)
        self.ctx.clear(0.5, 0.5, 0.5, 1.0)

        # Collect active gratings
        active = []
        if self.stim_one.contrast > 0:
            active.append(self.stim_one)
        if self.mask_one.contrast > 0:
            active.append(self.mask_one)
        if self.panel.show_stim_two and self.stim_two.contrast > 0:
            active.append(self.stim_two)
        if self.panel.show_mask_two and self.mask_two.contrast > 0:
            active.append(self.mask_two)

        # Single-pass additive compositing — no blending needed
        if active:
            self.renderer.render_gratings(active, vp_w, vp_h)

    def render_gui(self):
        """show_gui callback: imgui parameter panel."""
        now = time.perf_counter()
        dt = now - self.last_time  # approximate dt for UI

        self.panel.render(max(dt, 1.0 / 120.0))

        # Vector diagram overlay on the grating viewport
        io = imgui.get_io()
        viewport_w = io.display_size.x - self.panel.panel_width
        viewport_h = io.display_size.y
        self.panel.render_overlay(viewport_w, viewport_h)

        # Play/pause button and FPS in a small overlay
        imgui.set_next_window_pos(imgui.ImVec2(10, 10))
        imgui.set_next_window_size(imgui.ImVec2(200, 0))
        flags = (imgui.WindowFlags_.no_title_bar | imgui.WindowFlags_.no_resize |
                 imgui.WindowFlags_.always_auto_resize | imgui.WindowFlags_.no_saved_settings)
        imgui.begin("##overlay", True, flags)

        if imgui.button("Pause" if self.playing else "Play"):
            self.playing = not self.playing
        imgui.same_line()
        if imgui.button("Reset Phase"):
            for g in self.gratings:
                g.reset_phase()
        imgui.same_line()
        imgui.text(f"FPS: {self.fps_display}")

        imgui.end()


def main():
    app = App()

    runner_params = hello_imgui.RunnerParams()
    runner_params.app_window_params.window_geometry.size = (1600, 900)
    runner_params.app_window_params.window_title = "MonkeyML Grating Viewer"

    # Request OpenGL 4.1 Core profile
    runner_params.renderer_backend_options.open_gl_options = hello_imgui.OpenGlOptions()
    runner_params.renderer_backend_options.open_gl_options.major_version = 4
    runner_params.renderer_backend_options.open_gl_options.minor_version = 1
    runner_params.renderer_backend_options.open_gl_options.use_core_profile = True
    runner_params.renderer_backend_options.open_gl_options.use_forward_compat = True

    # Set callbacks
    runner_params.callbacks.post_init = app.post_init
    runner_params.callbacks.custom_background = app.render_background
    runner_params.callbacks.show_gui = app.render_gui

    # Use default backend (should pick OpenGL3 on macOS)
    runner_params.renderer_backend_type = hello_imgui.RendererBackendType.open_gl3

    # Disable idle throttling so animation renders continuously
    runner_params.fps_idling.enable_idling = False

    hello_imgui.run(runner_params)


if __name__ == '__main__':
    main()
