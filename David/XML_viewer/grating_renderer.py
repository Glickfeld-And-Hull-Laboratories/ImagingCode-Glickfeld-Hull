"""OpenGL grating rendering using moderngl."""

import math
import os
import moderngl
import numpy as np


class GratingInstance:
    """Per-grating CPU state — stores parameters and computes shader uniforms."""

    def __init__(self, name: str = ""):
        self.name = name
        # Visual parameters (degrees)
        self.direction_deg = 0.0
        self.spatial_freq = 0.05   # cycles per degree
        self.temporal_freq = 0.0   # cycles per second
        self.phase_deg = 0.0       # starting phase in degrees
        self.diameter_deg = 200.0  # size in degrees
        self.azimuth_deg = 0.0     # x position
        self.elevation_deg = 0.0   # y position
        self.contrast = 1.0        # alpha multiplier

        # Gaussian mask parameters (from XML stimulus definition)
        self.std_dev = 0.3
        self.mean = 0.1

        # Internal phase accumulator (radians)
        self._phase_rad = 0.0
        self._time = 0.0

    def compute_grating_coords(self):
        """
        Port of MWorks' grating coordinate math.
        Returns (bl, br, tl, tr) corner coordinates for the shader.
        """
        direction_rad = math.radians(self.direction_deg)
        phase_rad = math.radians(self.phase_deg) + self._phase_rad

        f = math.cos(direction_rad)
        g = math.sin(direction_rad)
        d = ((f + g) - 1.0) / 2.0

        bl = -d
        br = bl + f
        tl = bl + g
        tr = 1.0 + d

        # Scale by cycle_proportion (SF * diameter since x_size == y_size)
        cycle_prop = self.spatial_freq * self.diameter_deg
        bl *= cycle_prop
        br *= cycle_prop
        tl *= cycle_prop
        tr *= cycle_prop

        # Offset by phase
        phase_prop = phase_rad / (2.0 * math.pi)
        bl += phase_prop
        br += phase_prop
        tl += phase_prop
        tr += phase_prop

        return bl, br, tl, tr

    def update_phase(self, dt: float):
        """Advance phase based on temporal frequency. dt in seconds."""
        if self.temporal_freq != 0:
            # Phase advances by TF * dt cycles, i.e. TF * dt * 2pi radians
            self._phase_rad -= 2.0 * math.pi * self.temporal_freq * dt

    def reset_phase(self):
        """Reset accumulated phase to zero."""
        self._phase_rad = 0.0


class GratingRenderer:
    """Manages the moderngl shader program and renders gratings."""

    def __init__(self, ctx: moderngl.Context):
        self.ctx = ctx

        # Load shaders
        shader_dir = os.path.join(os.path.dirname(__file__), 'shaders')
        with open(os.path.join(shader_dir, 'grating.vert'), 'r') as f:
            vert_src = f.read()
        with open(os.path.join(shader_dir, 'grating.frag'), 'r') as f:
            frag_src = f.read()

        self.prog = ctx.program(vertex_shader=vert_src, fragment_shader=frag_src)

        # Create unit quad: 2 triangles covering (-1,-1) to (1,1)
        # Each vertex: (x, y, u, v)
        vertices = np.array([
            # Triangle 1
            -1.0, -1.0, 0.0, 0.0,
             1.0, -1.0, 1.0, 0.0,
             1.0,  1.0, 1.0, 1.0,
            # Triangle 2
            -1.0, -1.0, 0.0, 0.0,
             1.0,  1.0, 1.0, 1.0,
            -1.0,  1.0, 0.0, 1.0,
        ], dtype='f4')

        self.vbo = ctx.buffer(vertices.tobytes())
        self.vao = ctx.vertex_array(
            self.prog,
            [(self.vbo, '2f 2f', 'in_position', 'in_uv')],
        )

    def _set_grating_uniforms(self, grating: GratingInstance,
                              viewport_width: int, viewport_height: int,
                              reference_diameter: float = 200.0):
        """Compute and set all uniforms for a grating (except u_dev_sign)."""
        scale = grating.diameter_deg / reference_diameter

        aspect = viewport_width / max(viewport_height, 1)
        size_x = scale
        size_y = scale
        if aspect > 1.0:
            size_x /= aspect
        else:
            size_y *= aspect

        pos_x = (grating.azimuth_deg / (reference_diameter * 0.5))
        pos_y = (grating.elevation_deg / (reference_diameter * 0.5))
        if aspect > 1.0:
            pos_x /= aspect
        else:
            pos_y *= aspect

        bl, br, tl, tr = grating.compute_grating_coords()

        self.prog['u_position'].value = (pos_x, pos_y)
        self.prog['u_size'].value = (size_x, size_y)
        self.prog['u_gc_bl'].value = bl
        self.prog['u_gc_br'].value = br
        self.prog['u_gc_tl'].value = tl
        self.prog['u_gc_tr'].value = tr
        self.prog['u_contrast'].value = grating.contrast
        self.prog['u_std_dev'].value = grating.std_dev
        self.prog['u_mean'].value = grating.mean

    def render_grating(self, grating: GratingInstance, viewport_width: int, viewport_height: int,
                       reference_diameter: float = 200.0, num_gratings: int = 1):
        """Render a grating with two-pass additive deviation blending.

        Pass 1 (FUNC_ADD): adds bright deviations from mid-gray.
        Pass 2 (FUNC_REVERSE_SUBTRACT): subtracts dark deviations.
        Net effect: framebuffer += (grating - 0.5) * contrast * weight * gaussian_mask
        where weight = 1/num_gratings to prevent clipping when overlaying.
        """
        self._set_grating_uniforms(grating, viewport_width, viewport_height, reference_diameter)
        self.prog['u_weight'].value = 1.0 / max(num_gratings, 1)

        # Bright pass: add positive deviations
        self.ctx.blend_equation = moderngl.FUNC_ADD
        self.ctx.blend_func = (moderngl.SRC_ALPHA, moderngl.ONE)
        self.prog['u_dev_sign'].value = 1.0
        self.vao.render(moderngl.TRIANGLES)

        # Dark pass: subtract negative deviations
        self.ctx.blend_equation = moderngl.FUNC_REVERSE_SUBTRACT
        self.ctx.blend_func = (moderngl.SRC_ALPHA, moderngl.ONE)
        self.prog['u_dev_sign'].value = -1.0
        self.vao.render(moderngl.TRIANGLES)

        # Restore default blend state
        self.ctx.blend_equation = moderngl.FUNC_ADD

    def release(self):
        """Clean up OpenGL resources."""
        self.vao.release()
        self.vbo.release()
        self.prog.release()
