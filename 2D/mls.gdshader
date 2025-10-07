shader_type spatial;


const int n_balls = 8;
uniform vec3 particle_pos[8];

uniform float separation_limit = 2.0; // max distance to still connect

void vertex() {
    // Get particle index & cluster ID from vertex color
    uint particle_index = uint(COLOR.r * 65535.0 + COLOR.g * 65535.0 * 256.0);
    float cluster_id = COLOR.b;

    vec3 my_pos = particle_pos[particle_index];

    // Optional: smooth toward neighbor
    uint neighbor_index = particle_index + 1u; // naive neighbor
    if (int(neighbor_index) < n_balls) {
        vec3 neighbor_pos = particle_pos[neighbor_index];

        if (COLOR.b == cluster_id) { // only same cluster
            float d = distance(my_pos, neighbor_pos);
            if (d < separation_limit) {
                my_pos = mix(my_pos, neighbor_pos, 0.2);
            }
        }
    }

    VERTEX += (my_pos - VERTEX); // move vertex to particle
}
