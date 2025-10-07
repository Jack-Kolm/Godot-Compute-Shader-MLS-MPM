shader_type spatial;
render_mode blend_mix, unshaded, cull_back, world_vertex_coords; //, diffuse_burley, specular_schlick_ggx, world_vertex_coords;

uniform sampler2D tex;
uniform sampler2D blend_tex;

uniform int tex_size = 512;


float get_tex_val(sampler2D t, vec2 uv) {
	//return clamp(texture(t, uv).r, 0.0, 2000000.0);
	return texture(t, uv).r;
}


float smoothGaussian(vec2 uv, sampler2D t, float radius) {
    float sum = 0.0;
    float weightSum = 0.0;
    int r = int(radius);
    for (int x = -r; x <= r; x++) {
        for (int y = -r; y <= r; y++) {
            vec2 offset = vec2(float(x), float(y)) / float(tex_size);
            float dist2 = float(x*x + y*y);
            float w = exp(-dist2 / (2.0 * radius * radius)); // Gaussian weight
			float val = get_tex_val(t, uv+offset) * w;
            sum += val;//texture(t, uv + offset).r * w;
            weightSum += w;
        }
    }
	float final = sum / weightSum;
	return final;
}
varying vec2 my_uv;
varying vec3 vertex_point;
const float lower_limit = 0.05;
const float upper_limit = 1.0 - lower_limit;

void vertex() {
    //vec2 uv = (VERTEX.xz / 1.0) * 0.5 + vec2(0.5);
	//vec2 uv = VERTEX.xz + vec2(0.0);  // maps -0.5..0.5 → 0..1
	vec2 scale = vec2(0.9);    // zoom factor
	vec2 offset = vec2(0.05);  // move inside the 0–1 range
	vec2 uv = VERTEX.xz * vec2(1.0);//vec2(0.1, 0.1);
	uv = uv; //* scale + offset;
	float test_val = smoothGaussian(uv, tex, 28.0) * 1.0;
	float prev_test_val = smoothGaussian(uv, blend_tex, 28.0) * 1.0;
	float mixed = mix(test_val, prev_test_val, 0.5);
	float local_height_scale = 0.02;
	float local_height_offset = 0.0;
	float val = texture(tex, uv).r;
    float height = smoothstep(0.01, 0.5, test_val);
	
    if (VERTEX.y > 0.0) { // assuming the top face has positive y
		//if ((VERTEX.x > 0.1) && (VERTEX.x < 0.9)&& (VERTEX.z < 0.9)&& (VERTEX.z > 0.1)){
		VERTEX.y += mixed * local_height_scale;
		
		//if ((VERTEX.x < lower_limit) || (VERTEX.x > upper_limit) || (VERTEX.z < lower_limit) || (VERTEX.z > upper_limit)){
		//	VERTEX.y = 0.0; //VERTEX.y; //smoothGaussian(uv, tex, 28.0) * 1.0;
		//}
		//}
		//else{
		//	VERTEX = vec3(0.0, VERTEX.y, 0.0);
		//}
	}
	vertex_point = VERTEX;
	my_uv = uv;
}



void fragment() {
    //vec2 uv = (VERTEX.xz / 1.0) * 0.5 + vec2(0.5);
    vec2 uxv = UV + vec2(0.5); // but remember: BoxMesh UVs are per-face islands!
	//vec2 uv = VERTEX.xz + vec2(0.0);
	// vec2 uv = VERTEX.xz + vec2(0.5);
    float height_val = texture(tex, my_uv).r;
	height_val = smoothGaussian(my_uv, tex, 16.0) * 0.08;
	//float test_val = smoothGaussian(my_uv, tex, 12.0) * 1.0;

	//ALBEDO = vec3(test_val, 0.0, 0.0);
	//ALPHA = 0.5;


    // Water colors
    vec3 shallowestColor = vec3(0.1, 0.5, 0.55);

    vec3 shallowerColor = vec3(0.02, 0.3, 0.4);
    vec3 shallowColor = vec3(0.0, 0.2, 0.3);
    vec3 deepColor = vec3(0.0, 0.05, 0.2);

    // Approximate "depth" using world height
    //float depthFactor = clamp((VERTEX.y + 1.0) * 0.5, 0.0, 1.0);
    float depthFactor = -VERTEX.y;
	vec3 world_pos = (INV_VIEW_MATRIX * vec4(VERTEX, 1.0)).xyz;
	
	vec3 obj_pos = (inverse(MODEL_MATRIX) * vec4(world_pos, 1.0)).xyz;
    // Fresnel
    float fresnel = pow(1.0 - dot(normalize(NORMAL), normalize(CAMERA_POSITION_WORLD - VERTEX)), 3.0);

    // Foam: use heightmap texture
    //float foam = smoothstep(0.6, 1.0, height_val + fresnel * 0.5);


			//if ((VERTEX.x < lower_limit) || (VERTEX.x > upper_limit) || (VERTEX.z < lower_limit) || (VERTEX.z > upper_limit)){
// Limit foam to only near top surface (world Y close to water level)
    float foam = smoothstep(0.2, 0.6, height_val);

	float foam_surface = smoothstep(0.0, 0.2, world_pos.y);

	// Final foam factor (texture + fade with height)
	//foam_surface = smoothstep(0.3, -0.2, world_pos.y);
	//float foam = foam_t * foam_surface;
	vec3 foamer = vec3(0.0, world_pos.y, 0.0);
	float transparency = smoothstep(0.1, 0.0, world_pos.y);
    // Mix shallow vs deep
    vec3 waterColor = mix(deepColor, shallowColor, vec3(0.0, world_pos.y, 0.0));
    float aaalpha = mix(1.0, 0.5, world_pos.y);


	if (vertex_point.y > 1.0){//((world_pos.y > 1.0) && (cross(NORMAL, vec3(0.0, 1.0, 0.0)) == vec3(0.0, 0.0, 0.0))){
		waterColor = mix(waterColor, shallowerColor, vertex_point.y - 1.0);
		aaalpha = mix(0.5, 0.3, vertex_point.y - 1.0);
	}
	//if ((vertex_point.x <= lower_limit) || (vertex_point.x >= upper_limit) || (vertex_point.z <= lower_limit) || (vertex_point.z >= upper_limit)){
	    //foam = 0.0;//smoothstep(0.2, 0.6, 0.0);
	//	vec3 allo = mix(waterColor, shallowestColor, foam * 0.7);
	//	ALBEDO = mix(waterColor, allo, vertex_point.y + 0.65);
	//}
	//else{
	    vec3 n = normalize(NORMAL);                    // surface normal
    float upness = clamp(dot(n, vec3(0.0, 1.0, 0.0)), 0.0, 1.0); // 1 on top face, 0 on vertical
    float sideness = 1.0 - upness;                 // 1 on perfect vertical sides

    // ---- side fade (use per-face V coordinate if possible) ----
    // my_uv is the UV used for sampling the heightmap in your shader
    // on BoxMesh each face is an island so my_uv.y usually maps vertical for that face
    float face_v = my_uv.y;

    // if your side UVs are flipped, invert: face_v = 1.0 - my_uv.y;

    // control where the foam starts down the side (tweak)
    float side_fade_start = 0.65;   // start fading below this V (0..1)
    float side_fade = smoothstep(side_fade_start, 1.0, face_v);

    // reduce how strong side foam can be compared to top face
    float side_strength = 0.8;

    // combine top-face foam (texture-driven) and side foam (texture * vertical fade)
    float foam_top  = foam * upness;                               // top surface foam
    float foam_side = foam * sideness * side_fade * side_strength;  // side foam near top edge

    foam = clamp(foam_top + foam_side, 0.0, 1.0);

    // optional: bias foam by fresnel / rim so edges pop a little
    //fresnel = pow(1.0 - dot(n, normalize(CAMERA_POSITION_WORLD - VERTEX)), 2.0);
    foam *= mix(1.0, fresnel, 0.4); // 0.4 = how much fresnel affects foam (tweak)

    // ---- final color mix ----
    ALBEDO = mix(waterColor, shallowestColor, foam * 0.7);
    ALPHA  = aaalpha;
	//ALBEDO = mix(waterColor, shallowestColor, foam * 0.7);

	//}
	//ALPHA = aaalpha;

    // Final albedo with foam overlay
	//if (!(cross(NORMAL, vec3(0.0, 1.0, 0.0)) == vec3(0.0, 0.0, 0.0))){

    // Transparency
	//if ((world_pos.x < 0.1) || (world_pos.x > 0.9) || (world_pos.z < 0.1) || (world_pos.z > 0.9)){
	//	ALPHA = 0.0;
	//}
    //ALPHA = mix(0.9, 0.6, world_pos.y) + fresnel * 0.2;

    // Specular shine
    //SPECULAR = 0.9;
    //ROUGHNESS = 0.05 + foam * 0.4;
    //ALBEDO = vec3(height_val, 0, 0);
    // Get the world normal and position

    //uv *= uv_scale;
	//CAMERA_POSITION_WORLD

}
