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

layout(set = 0, binding = 0, std430) buffer Particles {
    Particle particles[];
} pbuffer;

layout(set = 0, binding = 1, std430) buffer Grid {
    Cell grid[];
} gridbuffer;

layout(set = 0, binding = 2, std430) buffer FS {
    vec3[] debug;
} fbuffer;


layout(set = 0, binding = 4, std430) buffer Move {
    vec3[] move_pos;
} movebuffer;


layout(push_constant) uniform Params {
    int tex_size;
    int num_particles;
    float dt;
    float gravity;
    float elastic_mu;
    float elastic_lambda;
    int grid_res_x;
    int grid_res_y;
    int grid_res_z;
} params;

layout(r32f, set = 0, binding = 3) uniform image2D current_image;


void g2p() {
    uint i = gl_GlobalInvocationID.x;
    if (i >= pbuffer.particles.length()) return;

    Particle p = pbuffer.particles[i];

    vec3 new_particle_velocity = vec3(0.0);
    vec3 center_cell_x = floor(p.x);
    vec3 cell_diff = p.x - (vec3(center_cell_x) + vec3(0.5));


    vec3 w[3] = {vec3(1.0, 1.0, 1.0), vec3(1.0, 1.0, 1.0), vec3(1.0, 1.0, 1.0)};
    w[0] = 0.5 * (vec3(0.5, 0.5, 0.5) - cell_diff) * (vec3(0.5, 0.5, 0.5) - cell_diff);
    w[1] =  vec3(0.75, 0.75, 0.75) - (cell_diff * cell_diff);
    w[2] = 0.5 * (vec3(0.5, 0.5, 0.5) + cell_diff) * (vec3(0.5, 0.5, 0.5) + cell_diff);


    mat3 B = mat3(vec3(0.), vec3(0.), vec3(0.));


    for (int gx = 0; gx < 3; ++gx)
    for (int gy = 0; gy < 3; ++gy)
    for (int gz = 0; gz < 3; ++gz) {
        float weight = w[gx].x * w[gy].y * w[gz].z;
        ivec3 cell_x = ivec3(center_cell_x.x + gx - 1, center_cell_x.y + gy - 1, center_cell_x.z + gz - 1);
        int cell_index = cell_x.x * params.grid_res_y * params.grid_res_z +
                         cell_x.y * params.grid_res_z +
                         cell_x.z;

        vec3 cell_dist = (vec3(cell_x) - pbuffer.particles[i].x) + vec3(0.5);

        vec3 weighted_velocity = vec3(gridbuffer.grid[cell_index].v.x * weight, 
                                      gridbuffer.grid[cell_index].v.y * weight, 
                                      gridbuffer.grid[cell_index].v.z * weight);


        B += mat3(weighted_velocity * cell_dist.x,
                  weighted_velocity * cell_dist.y,
                  weighted_velocity * cell_dist.z);
        new_particle_velocity += weighted_velocity;
    }
    mat3 C = B * 4.0f;
    pbuffer.particles[i].c1 = C[0];
    pbuffer.particles[i].c2 = C[1];
    pbuffer.particles[i].c3 = C[2];

    pbuffer.particles[i].v = new_particle_velocity;
    pbuffer.particles[i].x += new_particle_velocity * vec3(params.dt);
    pbuffer.particles[i].x = vec3(clamp(pbuffer.particles[i].x.x, 1.0, params.grid_res_x - 2.), 
                                  clamp(pbuffer.particles[i].x.y, 1.0, params.grid_res_y - 2.), 
                                  clamp(pbuffer.particles[i].x.z, 1.0, params.grid_res_z - 2.));
    vec3 wall_min = vec3(3.0);
    vec3 wall_max = vec3(float(params.grid_res_x - 4), float(params.grid_res_y - 4), float(params.grid_res_z - 4));
    float k = 3.0;
    

    vec3 x_n = pbuffer.particles[i].x + (pbuffer.particles[i].v);
    float wall_scale = 1.0;
    float wall_scale_x = 1.0;

    if (x_n.x < wall_min.x) {
        pbuffer.particles[i].v.x += (wall_min.x - x_n.x) * wall_scale;
    }
    if (x_n.x > wall_max.x) {
        pbuffer.particles[i].v.x += (wall_max.x - x_n.x) * wall_scale;
    }
    if (x_n.y < wall_min.y) {
        pbuffer.particles[i].v.y += (wall_min.y - x_n.y)* wall_scale;
        }
    if (x_n.y > wall_max.y) {
        pbuffer.particles[i].v.y += (wall_max.y - x_n.y)* wall_scale;
    }
    if (x_n.z < wall_min.z) {
        pbuffer.particles[i].v.z += (wall_min.z - x_n.z)* wall_scale;
        }
    if (x_n.z > wall_max.z) {
        pbuffer.particles[i].v.z += (wall_max.z - x_n.z)* wall_scale;
    }


    // Done with G2P, doing texture stuff
    p = pbuffer.particles[i];

    float box_min_x, box_min_y, box_min_z = 0;
    float box_max_x = params.grid_res_x;
    float box_max_y = params.grid_res_y;
    float box_max_z = params.grid_res_z;
    vec2 tex_size = vec2(params.tex_size, params.tex_size);
    vec2 normalised = vec2(
        (p.x.x) / (box_max_x),
        (p.x.z) / (box_max_z)
    );

    ivec2 uv = ivec2(
        int(clamp(normalised.x * tex_size.x, 0.0, tex_size.x)),
        int(clamp(normalised.y * tex_size.y, 0.0, tex_size.y))
    );
    vec4 current_value_at_uv = imageLoad(current_image, uv);
    if (p.x.y > current_value_at_uv.x){
        float particle_h_val = p.x.y;
        imageStore(current_image, uv, vec4(particle_h_val * 1.0, 0.0, 0.0, 0.0));

    }
    //fbuffer.debug[i] = vec3(uv.x, uv.y, 0); // for checking the debug buffer :)

}

void test(){
    for(int idx = 0; idx < (params.tex_size); idx++){
        for(int jdx = 0; jdx < params.tex_size; jdx++){
            ivec2 loop_uv = ivec2(idx, jdx);

            if (idx < (params.tex_size / 2)){
                imageStore(current_image, loop_uv, vec4(255.0, 0.0, 0.0, 0.0));

            }
            else{
                imageStore(current_image, loop_uv, vec4(0.0, 0.0, 0.0, 0.0));

            }
        }
    }
}

void main() {
    uint i = gl_GlobalInvocationID.x;
    if (i >= pbuffer.particles.length()) return;

    g2p();
    //Particle p = pbuffer.particles[i];

    // Map particle position x,z → texture coords
    //ivec2 uv = ivec2(
    //    clamp(int(p.x.x / params.grid_res_x * 8192), 0, 8192 - 1),
    //    clamp(int(p.x.z / params.grid_res_z * 8192), 0, 8192 - 1)
    //);
    //fbuffer.debug[i] = vec3(params.tex_size, 0, 0);
    
    //imageStore(current_image, uv, vec4(255.0, 0.0, 0.0, 0.0)); // RGBA red
    //imageStore(current_image, ivec2(100, 200), vec4(50.0, 30.0, 10.0, 0.5)); // RGBA red

// Example: store Y as depth in red channel
    //imageStore(current_image, uv, vec4(p.x.y, 0.0, 0.0, 1.0));
    // Use particle.y as height contribution
    //current_image[100] = 5.0;
    //ivec2 coord = ivec2(100, 200);
    //imageStore(current_image, uv, vec4(50.0, 30.0, 10.0, 0.5)); // RGBA red
    fbuffer.debug[0] = vec3(movebuffer.move_pos[0].x, movebuffer.move_pos[0].z, 0.0);
}
