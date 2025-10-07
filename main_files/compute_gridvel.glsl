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
    Cell grid[];
} gridbuffer;

layout(set = 0, binding = 2, std430) buffer FS {
    vec3[] debug;
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



void grid_velocity_update_n(uint index, uvec3 id3) {

    if (gridbuffer.grid[index].mass > 0) { 
        float mass = gridbuffer.grid[index].mass;

        vec3 float_v = vec3(
            gridbuffer.grid[index].v.x / mass, 
            gridbuffer.grid[index].v.y / mass, 
            gridbuffer.grid[index].v.z / mass
        );
        gridbuffer.grid[index].v.x = float_v.x;
        gridbuffer.grid[index].v.y = float_v.y + (params.gravity * params.dt);
        gridbuffer.grid[index].v.z = float_v.z;

        int x = int(id3.x);
        int y = int(id3.y);
        int z = int(id3.z);

        if (x < 2 || x > int(ceil(params.grid_res_x) - 3)) { gridbuffer.grid[index].v.x = 0; } 
        if (y < 2 || y > int(ceil(params.grid_res_y) - 3)) { gridbuffer.grid[index].v.y = 0; }
        if (z < 2 || z > int(ceil(params.grid_res_z) - 3)) { gridbuffer.grid[index].v.z = 0; }

    }
}

void main() {
    uvec3 id3 = gl_GlobalInvocationID;
    uint index = id3.x * params.grid_res_y * params.grid_res_z +
                 id3.y * params.grid_res_z +
                 id3.z;
    if (index >= gridbuffer.grid.length()) return;
    grid_velocity_update_n(index, id3);
}
