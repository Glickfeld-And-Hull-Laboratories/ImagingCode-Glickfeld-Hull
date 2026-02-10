#version 410 core

layout(location = 0) in vec2 in_position;  // unit quad: (-1,-1) to (1,1)
layout(location = 1) in vec2 in_uv;        // 0..1 UVs

// Grating quad placement in NDC
uniform vec2 u_position;   // NDC center of the grating
uniform vec2 u_size;        // NDC half-extents

// Grating texture coordinate corners (computed on CPU from direction, SF, phase)
uniform float u_gc_bl;  // bottom-left
uniform float u_gc_br;  // bottom-right
uniform float u_gc_tl;  // top-left
uniform float u_gc_tr;  // top-right

out float v_grating_coord;
out vec2 v_mask_coord;      // -1..1 for Gaussian mask

void main() {
    // Position the quad in NDC
    vec2 pos = in_position * u_size + u_position;
    gl_Position = vec4(pos, 0.0, 1.0);

    // Interpolate grating coordinate from 4 corners using UV
    float bottom = mix(u_gc_bl, u_gc_br, in_uv.x);
    float top    = mix(u_gc_tl, u_gc_tr, in_uv.x);
    v_grating_coord = mix(bottom, top, in_uv.y);

    // Mask coordinate: -1..1 from center
    v_mask_coord = in_uv * 2.0 - 1.0;
}
