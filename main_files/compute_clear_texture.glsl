#[compute]
#version 450
#extension GL_EXT_shader_atomic_float:enable


layout(local_size_x = 4, local_size_y = 4) in;


layout(r32f, set = 0, binding = 3) uniform image2D current_image;
layout(r32f, set = 0, binding = 1) uniform image2D blend_image;

layout(push_constant) uniform Params {
    int tex_size_x;
    int tex_size_y;
} params;

layout(set = 0, binding = 2, std430) buffer FS {
    vec3[] debug;
} fbuffer;



float lerp(float v0, float v1, float t) {
  return v0 + t * (v1 - v0);
}


void main() {
    uvec2 id2 = uvec2(gl_GlobalInvocationID.x, gl_GlobalInvocationID.y);
    ivec2 uv = ivec2(id2);
    vec4 point = imageLoad(current_image, uv);
    imageStore(blend_image, uv, point);
    float val = point.r;
    float new_val = lerp(val, 0.0, 0.1);
    imageStore(current_image, uv, vec4(new_val, 0.0, 0.0, 0.0)); 
}
