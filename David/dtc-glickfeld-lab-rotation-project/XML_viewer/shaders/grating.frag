#version 410 core

in vec2 v_ndc;

const int MAX_GRATINGS = 4;

uniform int   u_num_gratings;
uniform vec2  u_position[MAX_GRATINGS];   // NDC center
uniform vec2  u_size[MAX_GRATINGS];       // NDC half-extents
uniform float u_gc_bl[MAX_GRATINGS];      // grating coord corners
uniform float u_gc_br[MAX_GRATINGS];
uniform float u_gc_tl[MAX_GRATINGS];
uniform float u_gc_tr[MAX_GRATINGS];
uniform float u_contrast[MAX_GRATINGS];
uniform float u_std_dev[MAX_GRATINGS];
uniform float u_mean[MAX_GRATINGS];

out vec4 frag_color;

void main() {
    float result = 0.5;

    for (int i = 0; i < u_num_gratings; i++) {
        // Map NDC back to this grating's local (-1..1) space
        vec2 local = (v_ndc - u_position[i]) / u_size[i];

        // Convert to (0..1) UV for bilinear grating-coord interpolation
        vec2 uv = (local + 1.0) * 0.5;

        // Skip fragments outside this grating's quad
        if (uv.x < 0.0 || uv.x > 1.0 || uv.y < 0.0 || uv.y > 1.0) continue;

        // Bilinear interpolation of grating coordinate from 4 corners
        float bottom = mix(u_gc_bl[i], u_gc_br[i], uv.x);
        float top    = mix(u_gc_tl[i], u_gc_tr[i], uv.x);
        float grating = 0.5 * (1.0 + cos(6.283185307 * mix(bottom, top, uv.y)));

        // Gaussian mask in local space
        float dist = length(local);
        float mask = exp(-0.5 * pow((dist - u_mean[i]) / u_std_dev[i], 2.0));

        // Additive: accumulate signed deviation (no 1/N weight)
        result += (grating - 0.5) * u_contrast[i] * mask;
    }

    frag_color = vec4(vec3(clamp(result, 0.0, 1.0)), 1.0);
}
