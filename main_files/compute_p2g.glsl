#[compute]
#version 450
#extension GL_EXT_shader_atomic_float:enable


layout(local_size_x = 4) in;

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


layout(set=0, binding = 0, std430) buffer Particles {
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




void p2g() {
    uint i = gl_GlobalInvocationID.x;
    if (i >= pbuffer.particles.length()) return;

    Particle p = pbuffer.particles[i];


    vec3 center_cell_x = floor(pbuffer.particles[i].x); // x as in positional vector, not x axis
    vec3 cell_diff = p.x - (vec3(center_cell_x) + vec3(0.5));

    vec3 w[3];
    w[0] = 0.5 * ((vec3(0.5, 0.5, 0.5) - cell_diff) * (vec3(0.5, 0.5, 0.5) - cell_diff));
    w[1] =  vec3(0.75, 0.75, 0.75) - (cell_diff * cell_diff);
    w[2] = 0.5 * (vec3(0.5, 0.5, 0.5) + cell_diff) * (vec3(0.5, 0.5, 0.5) + cell_diff);



    for (int gx = 0; gx < 3; ++gx)
    for (int gy = 0; gy < 3; ++gy)
    for (int gz = 0; gz < 3; ++gz) {
        float weight = w[gx].x * w[gy].y * w[gz].z;
        ivec3 cell_x = ivec3(center_cell_x.x + gx - 1, center_cell_x.y + gy - 1, center_cell_x.z + gz - 1);
        vec3 cell_dist = (vec3(cell_x) + vec3(0.5)) - p.x;
        int cell_index = cell_x.x * params.grid_res_y * params.grid_res_z +
                         cell_x.y * params.grid_res_z +
                         cell_x.z;
        mat3 C = mat3(p.c1,
                      p.c2,
                      p.c3);

        vec3 Q = C * cell_dist;
        float weighted_mass_contribution = weight * p.mass;
        vec3 momentum_contribution = weighted_mass_contribution * (p.v + Q);

        atomicAdd(gridbuffer.grid[cell_index].mass, weighted_mass_contribution);
        atomicAdd(gridbuffer.grid[cell_index].v[0], momentum_contribution.x);
        atomicAdd(gridbuffer.grid[cell_index].v[1], momentum_contribution.y);
        atomicAdd(gridbuffer.grid[cell_index].v[2], momentum_contribution.z);

    }
}


void main() {
    uint i = gl_GlobalInvocationID.x;
    if (i >= pbuffer.particles.length()) return;
    p2g();
}
