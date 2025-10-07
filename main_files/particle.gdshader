
shader_type particles;

uniform vec3 particles[3199];

void start() {
	// Called when a particle is spawned.
}

void process() {
	vec3 offset = vec3(-7, -4, -22);
	vec3 position = particles[INDEX];
	TRANSFORM[3].xyz = (position + offset) * 0.01;
	
	//SCAL
	//vec3 current_pos = TRANSFORM[3].xyz;
	//vec3 target_pos = particles[INDEX];
	//TRANSFORM[3].xyz = mix(current_pos, target_pos, DELTA);
	// Called every frame on existing particles (according to the Fixed FPS property).
}
