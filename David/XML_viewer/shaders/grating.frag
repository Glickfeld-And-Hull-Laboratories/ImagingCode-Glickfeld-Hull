#version 410 core

in float v_grating_coord;
in vec2 v_mask_coord;

uniform float u_contrast;    // grating contrast (0..1)
uniform float u_std_dev;     // Gaussian mask std_dev (default 0.3)
uniform float u_mean;        // Gaussian mask mean (default 0.1)
uniform float u_dev_sign;    // +1.0 for bright pass, -1.0 for dark pass
uniform float u_weight;      // 1.0/N where N = number of active gratings

out vec4 frag_color;

void main() {
    // Sinusoidal grating: 0.5 * (1 + cos(2*pi*coord))
    float grating = 0.5 * (1.0 + cos(6.283185307 * v_grating_coord));

    // Gaussian mask
    float dist = length(v_mask_coord);
    float mask_val = exp(-0.5 * pow((dist - u_mean) / u_std_dev, 2.0));

    // Deviation from mid-gray, scaled by contrast, weight, and pass sign
    // Weight = 1/N prevents clipping when multiple gratings overlap
    // Bright pass (u_dev_sign=+1): outputs max(deviation, 0) — added to framebuffer
    // Dark pass  (u_dev_sign=-1): outputs max(-deviation, 0) — subtracted from framebuffer
    float deviation = (grating - 0.5) * u_contrast * u_weight * u_dev_sign;
    float val = max(deviation, 0.0);

    frag_color = vec4(val, val, val, mask_val);
}
