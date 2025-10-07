#[compute]
#version 450
#extension GL_EXT_shader_atomic_float:enable


layout(local_size_x = 4, local_size_y = 4, local_size_z = 4) in;

struct Particle {
    vec3 x;  // position
    float pad0;   

    float mass;
    float pad1; 
    float pad2;  
    float pad3; 

    vec3 v;   
    float pad4;  

    vec3 c1;
    float pad5;
    vec3 c2; 
    float pad6; 
    vec3 c3; 
    float pad7; 
};

struct Cell {
    vec3 v;       
    float pad0;   

    float mass; 
    float pad1; 
    float pad2; 
    float pad3;  
};

layout(set = 0, binding = 0, std430) buffer Particles {
    Particle particles[];
} pbuffer;

layout(set = 0, binding = 1, std430) buffer Grid {
    //float dummy[];
    Cell grid[];
} gridbuffer;

layout(set = 0, binding = 2, std430) buffer FS {
    vec3 debug;
} fbuffer;

layout(push_constant) uniform Params {
    int grid_res;
    int num_particles;
    float dt;
    float gravity;
    float elastic_mu;
    float elastic_lambda;
    int grid_res_x;
    int grid_res_y;
    int grid_res_z;
} params;



void main() {
    uvec3 id3 = gl_GlobalInvocationID;
    
    uint index = id3.x * params.grid_res_y * params.grid_res_z + 
                 id3.y * params.grid_res_z  + 
                 id3.z;

    if (index >= gridbuffer.grid.length()) return;

    //uint i = gl_GlobalInvocationID.x;
    gridbuffer.grid[index].mass = 0.0;
    gridbuffer.grid[index].v.x = 0.0;
    gridbuffer.grid[index].v.y = 0.0;
    gridbuffer.grid[index].v.z = 0.0;


}
