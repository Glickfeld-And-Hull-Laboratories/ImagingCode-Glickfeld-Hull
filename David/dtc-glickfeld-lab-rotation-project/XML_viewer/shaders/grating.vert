#version 410 core

layout(location = 0) in vec2 in_position;  // fullscreen quad: (-1,-1) to (1,1)

out vec2 v_ndc;

void main() {
    gl_Position = vec4(in_position, 0.0, 1.0);
    v_ndc = in_position;
}
