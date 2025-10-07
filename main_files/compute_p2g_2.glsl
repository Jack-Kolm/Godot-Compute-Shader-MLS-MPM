#[compute]
#version 450
#extension GL_EXT_shader_atomic_float:enable


layout(local_size_x = 4) in;
// we now have 2 P2G phases as we need to ensure we have scattered particle masses to the grid once,
// in order to get our density estimate at each frame

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



void p2g2() {
    uint i = gl_GlobalInvocationID.x;
    if (i >= pbuffer.particles.length()) return;

    Particle p = pbuffer.particles[i];


    vec3 center_cell_x = floor(p.x);
    vec3 cell_diff = p.x - (vec3(center_cell_x) + vec3(0.5));

    vec3 w[3] = {vec3(1.0, 1.0, 1.0), vec3(1.0, 1.0, 1.0), vec3(1.0, 1.0, 1.0)};
    w[0] = 0.5 * (vec3(0.5, 0.5, 0.5) - cell_diff) * (vec3(0.5, 0.5, 0.5) - cell_diff);
    w[1] =  vec3(0.75, 0.75, 0.75) - (cell_diff * cell_diff);
    w[2] = 0.5 * (vec3(0.5, 0.5, 0.5) + cell_diff) * (vec3(0.5, 0.5, 0.5) + cell_diff);

    float density = 0;

    for (int gx = 0; gx < 3; ++gx)
    for (int gy = 0; gy < 3; ++gy)
    for (int gz = 0; gz < 3; ++gz) {
        float weight = w[gx].x * w[gy].y * w[gz].z;
        ivec3 cell_x = ivec3(center_cell_x.x + gx - 1, center_cell_x.y + gy - 1, center_cell_x.z + gz - 1);
        int cell_index = cell_x.x * params.grid_res_y * params.grid_res_z +
                         cell_x.y * params.grid_res_z +
                         cell_x.z;
        density += gridbuffer.grid[cell_index].mass * weight;
    }
    float elastic_lambda = 4.0; //note to self: defined here instead of using params.elastic_lambda; //originally 4.0
    float elastic_mu = 0.2; //note to self: defined here instead of using params.elastic_mu; // originally 0.1 or 0.2
    float eos_stiffness = 1.0; // stiffness
    float volume = 1.0 / density; 
    float pressure = max(-0.0, eos_stiffness * (pow(density/elastic_lambda, 5.0) - 1.0));
    mat3 stress = mat3(-pressure, 0, 0,
                       0, -pressure, 0,
                       0, 0.0, -pressure);
    mat3 C = mat3(pbuffer.particles[i].c1,
                  pbuffer.particles[i].c2,
                  pbuffer.particles[i].c3); 
    mat3 dudv = C;
    mat3 strain = dudv + transpose(dudv);
    stress = stress + (elastic_mu * strain);
    mat3 eq_16_term_0 = (-volume * 4 *  params.dt) * stress;

    
    for (int gx = 0; gx < 3; ++gx)
    for (int gy = 0; gy < 3; ++gy)
    for (int gz = 0; gz < 3; ++gz) {
        float weight = w[gx].x * w[gy].y * w[gz].z;
        ivec3 cell_x = ivec3(center_cell_x.x + gx - 1, center_cell_x.y + gy - 1, center_cell_x.z + gz - 1);
        vec3 cell_dist = (vec3(cell_x) - pbuffer.particles[i].x) + vec3(0.5);
        
        int cell_index = cell_x.x * params.grid_res_y * params.grid_res_z +
                         cell_x.y * params.grid_res_z +
                         cell_x.z;

        vec3 momentum = eq_16_term_0 * weight * cell_dist;
        atomicAdd(gridbuffer.grid[cell_index].v.x, momentum.x);
        atomicAdd(gridbuffer.grid[cell_index].v.y, momentum.y);
        atomicAdd(gridbuffer.grid[cell_index].v.z, momentum.z);

    }

}


void main() {
    uint i = gl_GlobalInvocationID.x;
    if (i >= pbuffer.particles.length()) return;
    p2g2();
}
