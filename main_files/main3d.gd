#@tool
extends Node3D
const RUN_P2G = true
const RUN_GRIDVEL = true
const RUN_G2P = true
const floor_height = 50

#const grid_res : int = 64

#const grid_resos : int = 64

#const num_cells :int = grid_res * grid_res * grid_res

const gravity : float = -0.9 # og 0.2 or 0.9, between is good

const elastic_lambda : float = 10.5 #0.8
const elastic_mu : float = 10.0
const rest_density = 4.0
const dynamic_viscosity = 0.2
# equation of state
const eos_stiffness = 3.0
const eos_power = 4
#var num_particles : int
var all_num_particles : int

var all_particles : Array = []
var grid : Array = []
var Fs : Array = []
var weights: Array[Vector3] = [
	Vector3.ZERO,
	Vector3.ZERO,
	Vector3.ZERO
]

var manual_dt = 0.2 # reasonably between 1.0 and 0.05, higher means more aggressive and faster simulation
var iterations : int = int(1.0/manual_dt)

var thread1 : Thread


var rd: RenderingDevice

var shader_is_prepped = false
var parameters : PackedByteArray
var compute_parameters : Dictionary


func _ready():

	await init_3()
	await create_visual_particles_basic()
	RenderingServer.call_on_render_thread(shader_prep.bind())
	var material: ShaderMaterial = $MeshInstance3D.material_override
	if material:
		rdtexture = material.get_shader_parameter("tex")
		rdtex_blend = material.get_shader_parameter("blend_tex")
		material.set_shader_parameter("tex_size", tex_size.x)

func _process(delta):
	if shader_is_prepped:
		
		RenderingServer.call_on_render_thread(compute_shader_iteration.bind(delta))
		if rdtexture:
			rdtexture.texture_rd_rid = rdtexture_rid
			rdtex_blend.texture_rd_rid = rdtex_blend_rid

var visual_particles = []
func create_visual_particles_basic():
	var count = 0
	for i in range(all_num_particles):
		if count == 64:
			var p = Particle3D.new()  # or preload/instantiate a scene, or create a struct/dictionary
			p.x = all_particles[i].x
			p.v = all_particles[i].v
			p.C = all_particles[i].C
			p.mass = all_particles[i].mass
			p.index = i
			visual_particles.append(p)
			count = 0
			#self.add_child(p)
		all_particles[i].queue_free()
		count += 1
	print(len(visual_particles))


var grid_res_v : Vector3i = Vector3i(64, 64, 64)
var num_cells :int = grid_res_v.x * grid_res_v.y * grid_res_v.z
func make_alot_of_particles():
	var inner_box_size : Vector3 = Vector3(48, 48, 48)

	var grid_center = grid_res_v / 2.0
	var inner_box_min = grid_center - inner_box_size / 2.0
	var inner_box_max = grid_center + inner_box_size / 2.0
	var particles = 102400 * 2
	var particles_per_axis = int(round(pow(particles, 1.0/3.0))) # ~12.7 -> 13
	var particles_list = []
	var step = inner_box_size / float(particles_per_axis - 1)

	for x in range(particles_per_axis):
		for y in range(particles_per_axis):
			for z in range(particles_per_axis):
				var pos = inner_box_min + Vector3(x, y, z) * step
				particles_list.append(pos)
				if particles_list.size() >= particles:
					break
			if particles_list.size() >= particles:
				break
		if particles_list.size() >= particles:
			break
			
	return particles_list


func init_3():
	var temp_positions = make_alot_of_particles()
	all_num_particles = len(temp_positions)
	var rng = RandomNumberGenerator.new()

	print(all_num_particles)
	for i in all_num_particles:
		var p = Particle3D.new()  # or preload/instantiate a scene, or create a struct/dictionary
		p.x = temp_positions[i];

		var speed_x = rng.randf_range(0.0, 7.0)
		var speed_z = rng.randf_range(-1.0, 1.0)
		p.v = Vector3(0.0, 0.0, 0.0)
		p.mass = 1.0
		all_particles.append(p)
	for i in num_cells:
		var c = Cell3D.new()
		c.v = Vector3(0.0, 0.0, 0.0)
		grid.append(c)
	return


var rdtexture : Texture2DRD
var rdtexture_rid : RID
var rdtex_blend : Texture2DRD
var rdtex_blend_rid : RID
var rdtexture_rids: Array[RID] = [RID(), RID(), RID()]
var tex_size = Vector2i(1024, 1024)
var clear_tex_params

func create_shader(file_path="res://main_files/compute_p2g.glsl"):
	var shader_file := load(file_path)
	var shader_spirv : RDShaderSPIRV = shader_file.get_spirv()
	return rd.shader_create_from_spirv(shader_spirv)
	

func shader_prep():
	rd = RenderingServer.get_rendering_device()
	for i in num_cells: # prep the cells
		var cell = grid[i]
		cell.mass = 0.0
		cell.v = Vector3(0, 0, 0)
		grid[i] = cell # needed?
	var p2g_shader = create_shader("res://main_files/compute_p2g.glsl")
	var p2g2_shader = create_shader("res://main_files/compute_p2g_2.glsl")

	var gridvel_shader = create_shader("res://main_files/compute_gridvel.glsl")
	var g2p_shader = create_shader("res://main_files/compute_g2p.glsl")
	var zero_shader = create_shader("res://main_files/compute_zerogrid.glsl")

	var particles_constants := PackedFloat32Array() #PackedVector3Array(particles_positions()) #
	particles_constants = pack_particles_std430(all_particles)

	var particle_bytes := particles_constants.to_byte_array()

	var particle_buffer = rd.storage_buffer_create(particle_bytes.size(), particle_bytes)
	var particle_uniform := RDUniform.new()
	particle_uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_STORAGE_BUFFER
	particle_uniform.binding = 0 # this needs to match the binding in our shader file
	particle_uniform.add_id(particle_buffer)
	
	var grid_constants := pack_cells_std430(grid) 
	var grid_bytes := grid_constants.to_byte_array()
	var grid_buffer = rd.storage_buffer_create(grid_bytes.size(), grid_bytes)
	var grid_uniform := RDUniform.new()
	grid_uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_STORAGE_BUFFER
	grid_uniform.binding = 1 
	grid_uniform.add_id(grid_buffer)
	
	# Buffer for Fs[] # repurposed as a debug buffer
	var f_constants := pack_Fs_std430(Fs)
	var f_bytes := f_constants.to_byte_array()
	var f_buffer = rd.storage_buffer_create(f_bytes.size(), f_bytes)
	var f_uniform := RDUniform.new()
	f_uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_STORAGE_BUFFER
	f_uniform.binding = 2
	f_uniform.add_id(f_buffer)
	
	# Move buffer # i never fully implemented this, please ignore!
	var move_constants := pack_move_std430()
	var move_bytes := move_constants.to_byte_array()
	var move_buffer = rd.storage_buffer_create(move_bytes.size(), move_bytes)
	var move_uniform := RDUniform.new()
	move_uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_STORAGE_BUFFER
	move_uniform.binding = 4
	move_uniform.add_id(move_buffer)
	
	# Instancing the texture uniform
	var tf: RDTextureFormat = RDTextureFormat.new()
	tf.format = RenderingDevice.DATA_FORMAT_R32_SFLOAT
	tf.texture_type = RenderingDevice.TEXTURE_TYPE_2D
	tf.width = tex_size.x
	tf.height = tex_size.y
	tf.depth = 1
	tf.array_layers = 1
	tf.mipmaps = 1
	tf.usage_bits = (
			RenderingDevice.TEXTURE_USAGE_SAMPLING_BIT |
			RenderingDevice.TEXTURE_USAGE_COLOR_ATTACHMENT_BIT |
			RenderingDevice.TEXTURE_USAGE_STORAGE_BIT |
			RenderingDevice.TEXTURE_USAGE_CAN_UPDATE_BIT |
			RenderingDevice.TEXTURE_USAGE_CAN_COPY_TO_BIT |
			RenderingDevice.TEXTURE_USAGE_CAN_COPY_FROM_BIT
	)
	rdtexture_rid = rd.texture_create(tf, RDTextureView.new(), [])
	rdtex_blend_rid = rd.texture_create(tf, RDTextureView.new(), [])
	rd.texture_clear(rdtexture_rid, Color(0, 0, 0, 0), 0, 1, 0, 1)
	rd.texture_clear(rdtex_blend_rid, Color(0, 0, 0, 0), 0, 1, 0, 1)

	var texture_uniform := RDUniform.new()
	texture_uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_IMAGE
	texture_uniform.binding = 3
	texture_uniform.add_id(rdtexture_rid)
	var blend_texture_uniform := RDUniform.new()
	blend_texture_uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_IMAGE
	blend_texture_uniform.binding = 1
	blend_texture_uniform.add_id(rdtex_blend_rid)


	
	# Create sets
	# these could maybe be reused, but this may be clearer
	var uniform_set_p2g = rd.uniform_set_create([particle_uniform, grid_uniform, f_uniform], p2g_shader, 0) # the last parameter (the 0) needs to match the set in our shader file
	var uniform_set_p2g2 = rd.uniform_set_create([particle_uniform, grid_uniform, f_uniform], p2g2_shader, 0) 

	var uniform_set_gridvel = rd.uniform_set_create([particle_uniform, grid_uniform, f_uniform], gridvel_shader, 0) 
	var uniform_set_g2p = rd.uniform_set_create([particle_uniform, grid_uniform, f_uniform, texture_uniform, move_uniform], g2p_shader, 0) 
	var uniform_set_zero = rd.uniform_set_create([particle_uniform, grid_uniform, f_uniform], zero_shader, 0) 

	# Create pipelines
	var pipeline_p2g = rd.compute_pipeline_create(p2g_shader)
	var pipeline_p2g2 = rd.compute_pipeline_create(p2g2_shader)

	var pipeline_gridvel = rd.compute_pipeline_create(gridvel_shader)
	var pipeline_g2p = rd.compute_pipeline_create(g2p_shader)
	var pipeline_zero = rd.compute_pipeline_create(zero_shader)
	
	var clear_tex_shader = create_shader("res://main_files/compute_clear_texture.glsl") 	# additional pipeline for clearing the texture properly
	var uniform_set_clear_tex = rd.uniform_set_create([f_uniform, texture_uniform, blend_texture_uniform], clear_tex_shader, 0)
	var pipeline_clear_tex = rd.compute_pipeline_create(clear_tex_shader)


	var params_array := PackedFloat32Array([
		0.0,
		manual_dt,
		gravity,
		elastic_mu,
		elastic_lambda,
		grid_res_v.x, #grid_res_x,
		grid_res_v.y, #grid_res_y,
		grid_res_v.z  #grid_res_z,
	])
	parameters = PackedByteArray()
	parameters.resize(48)
	parameters.encode_s32(0, tex_size.x)
	parameters.encode_s32(4, all_num_particles)

	parameters.encode_float(8, manual_dt)
	parameters.encode_float(12, gravity)
	parameters.encode_float(16, elastic_mu)
	parameters.encode_float(20, elastic_lambda)
	parameters.encode_s32(24, int(grid_res_v.x))
	parameters.encode_s32(28, int(grid_res_v.y))
	parameters.encode_s32(32, int(grid_res_v.z))

	clear_tex_params = PackedByteArray()
	clear_tex_params.resize(16)
	clear_tex_params.encode_s32(0, tex_size.x)
	clear_tex_params.encode_s32(4, tex_size.y)



	var gvu_shader_local_size = Vector3(4, 4, 4) # same as size in the shader
	var grid_workgroups = Vector3i(grid_res_v.x / gvu_shader_local_size.x,
									grid_res_v.y / gvu_shader_local_size.y,
									grid_res_v.z / gvu_shader_local_size.z)
	
	var particle_workgroups = int(all_num_particles / 4)


	compute_parameters = {
		"pipeline_p2g": pipeline_p2g,
		"pipeline_p2g2": pipeline_p2g2,
		"pipeline_gridvel": pipeline_gridvel,
		"pipeline_g2p": pipeline_g2p,
		"pipeline_clear_tex": pipeline_clear_tex,
		"uniform_set_p2g": uniform_set_p2g,
		"uniform_set_p2g2": uniform_set_p2g2,
		"uniform_set_gridvel": uniform_set_gridvel,
		"uniform_set_g2p": uniform_set_g2p,
		"uniform_set_clear_tex": uniform_set_clear_tex,
		"grid_buffer": grid_buffer,
		"particle_buffer": particle_buffer,
		"f_buffer": f_buffer,
		"move_buffer": move_buffer,
		"grid_bytes": grid_bytes,
		"pipeline_zero": pipeline_zero,
		"uniform_set_zero": uniform_set_zero,
		"particle_workgroups": particle_workgroups,
		"grid_workgroups": grid_workgroups
	}
	shader_is_prepped = true



func compute_shader_iteration(dt):
	"""
	(clear texture ->) zero the grid -> P2G -> grid velocity update -> P2G -> repeat
	"""

	## clear tex
	var clear_tex_workgroup_size = tex_size.x / 4
	var compute_list := rd.compute_list_begin()
	rd.compute_list_bind_compute_pipeline(compute_list, compute_parameters.pipeline_clear_tex)
	rd.compute_list_bind_uniform_set(compute_list, compute_parameters.uniform_set_clear_tex, 0)
	rd.compute_list_set_push_constant(compute_list, clear_tex_params, clear_tex_params.size())
	rd.compute_list_dispatch(compute_list, clear_tex_workgroup_size, clear_tex_workgroup_size, 1)
	rd.compute_list_end()

	## zero out grid
	compute_list = rd.compute_list_begin()
	rd.compute_list_bind_compute_pipeline(compute_list, compute_parameters.pipeline_zero)
	rd.compute_list_bind_uniform_set(compute_list, compute_parameters.uniform_set_zero, 0)
	rd.compute_list_set_push_constant(compute_list, parameters, parameters.size())
	rd.compute_list_dispatch(compute_list, compute_parameters.grid_workgroups.x, compute_parameters.grid_workgroups.y, compute_parameters.grid_workgroups.z)
	rd.compute_list_end()

	## p2g
	compute_list = rd.compute_list_begin()
	rd.compute_list_bind_compute_pipeline(compute_list, compute_parameters.pipeline_p2g)
	rd.compute_list_bind_uniform_set(compute_list, compute_parameters.uniform_set_p2g, 0)
	rd.compute_list_set_push_constant(compute_list, parameters, parameters.size())
	rd.compute_list_dispatch(compute_list, compute_parameters.particle_workgroups, 1, 1)
	rd.compute_list_end()
	
	## p2g - 2
	compute_list = rd.compute_list_begin()
	rd.compute_list_bind_compute_pipeline(compute_list, compute_parameters.pipeline_p2g2)
	rd.compute_list_bind_uniform_set(compute_list, compute_parameters.uniform_set_p2g2, 0)
	rd.compute_list_set_push_constant(compute_list, parameters, parameters.size())
	rd.compute_list_dispatch(compute_list, compute_parameters.particle_workgroups, 1, 1)
	rd.compute_list_end()

	## update grid velocity
	compute_list = rd.compute_list_begin()
	rd.compute_list_bind_compute_pipeline(compute_list, compute_parameters.pipeline_gridvel)
	rd.compute_list_bind_uniform_set(compute_list, compute_parameters.uniform_set_gridvel, 0)
	rd.compute_list_set_push_constant(compute_list, parameters, parameters.size())
	rd.compute_list_dispatch(compute_list, compute_parameters.grid_workgroups.x, compute_parameters.grid_workgroups.y, compute_parameters.grid_workgroups.z)
	rd.compute_list_end()

		## g2p
	compute_list = rd.compute_list_begin()
	rd.compute_list_bind_compute_pipeline(compute_list, compute_parameters.pipeline_g2p)
	rd.compute_list_bind_uniform_set(compute_list, compute_parameters.uniform_set_g2p, 0)
	rd.compute_list_set_push_constant(compute_list, parameters, parameters.size())
	rd.compute_list_dispatch(compute_list, compute_parameters.particle_workgroups, 1, 1)
	rd.compute_list_end()
	
	## Iteration complete, get updated particles from buffer and sync with CPU particles (if needed)
	#var particle_output = rd.buffer_get_data(compute_parameters.particle_buffer).to_float32_array()
	#var updated_particles = unpack_particles_std430(particle_output)

	## Debug
	#debug_particles_std430(particle_output)
	#var debug_output = rd.buffer_get_data(compute_parameters.f_buffer).to_float32_array()
	#RenderingServer.call_on_render_thread(unpack_Fs_std430.bind(debug_output))

	## Move buffer
	#var data = repack_character().to_byte_array()
	#rd.buffer_update(compute_parameters.move_buffer, 0.0, data.size(), data)
	#var debug_output = rd.buffer_get_data(compute_parameters.move_buffer).to_float32_array()
	#RenderingServer.call_on_render_thread(unpack_move_std430.bind(debug_output))

	return




func get_texel_value(tex_rid: RID, width: int, height: int, x: int, y: int) -> float:
	var data: PackedByteArray = rd.texture_get_data(tex_rid, 0)
	var offset = (y * width + x) * 4
	var res = data.decode_float(offset)
	print(x, " ", y, " ", " that is ", res)
	return res


func debug_print_texture(tex_rid: RID, w: int, h: int, max_samples: int = 10):
	var data: PackedByteArray = rd.texture_get_data(tex_rid, 0)
	var idx = 0
	for y in range(h):
		for x in range(w):
			var f = data.decode_float(idx * 4)
			if idx < max_samples:
				print("(", x, ",", y, ") = ", f)
			idx += 1


func pack_particles_std430(all_particles_: Array) -> PackedFloat32Array:
	var result := PackedFloat32Array()
	for p in all_particles_:
		# x: vec3 (with padding)
		result.push_back(p.x.x)
		result.push_back(p.x.y)
		result.push_back(p.x.z)
		result.push_back(0.0)

		result.push_back(p.mass)
		result.push_back(0.0) 
		result.push_back(0.0)
		result.push_back(0.0)

		# v: vec3 (with padding)
		result.push_back(p.v.x)
		result.push_back(p.v.y)
		result.push_back(p.v.z)
		result.push_back(0.0) 

		# C: mat3 as 3 padded vec3 columns
		result.push_back(0.0)
		result.push_back(0.0)
		result.push_back(0.0)
		result.push_back(0.0)
		result.push_back(0.0)
		result.push_back(0.0)
		result.push_back(0.0)
		result.push_back(0.0)
		result.push_back(0.0)
		result.push_back(0.0)
		result.push_back(0.0)
		result.push_back(0.0)
	return result


func pack_cells_std430(cells: Array) -> PackedFloat32Array:
	var result := PackedFloat32Array()

	for c in cells:
		result.push_back(c.v.x)
		result.push_back(c.v.y)
		result.push_back(c.v.z)
		result.push_back(0.0)

		result.push_back(c.mass)
		result.push_back(0.0) 
		result.push_back(0.0) 
		result.push_back(0.0)
	return result

func pack_Fs_std430(cells: Array) -> PackedFloat32Array:
	"""
	Turned into debug buffer
	"""
	var result := PackedFloat32Array()
	for i in range(all_num_particles): #	for f in Fs:
		result.push_back(-1.0)
		result.push_back(-1.0)
		result.push_back(-1.0)
		result.push_back(0.0)

	return result

func pack_move_std430() -> PackedFloat32Array:
	"""
		For player collision
	"""
	var result := PackedFloat32Array()
	for i in range(10):

		result.push_back(-1.0)
		result.push_back(-1.0)
		result.push_back(-1.0)
		result.push_back(0.0)
	return result

func repack_character() -> PackedFloat32Array:
	"""
		For player collision
	"""
	var result := PackedFloat32Array()
	var char = $CharacterBody3D
	result.push_back(char.position.x)
	result.push_back(char.position.y)
	result.push_back(char.position.z)
	result.push_back(0.0)
	return result


func unpack_Fs_std430(data: PackedFloat32Array) -> Array:
	"""
	Turned into debug buffer
	"""
	var result :=  []
	var stride := 4
	for i in range(0, data.size(), stride):
		var val := Vector3(data[i + 0], data[i + 1], data[i + 2])
		#if val.x < 50:
		#	print(val)
		result.append(val)
	#print(result)
	var particle_index_to_print = 1254
	#print(len(result))
	#print("debug vals at index ", int(particle_index_to_print), ": ", result[particle_index_to_print])
	if rdtex_blend_rid:
		print("FIRST")
		RenderingServer.call_on_render_thread(get_texel_value.bind(rdtex_blend_rid, tex_size.x, tex_size.y, int(result[particle_index_to_print].x), int(result[particle_index_to_print].y)))
		print("SECOND")
		RenderingServer.call_on_render_thread(get_texel_value.bind(rdtexture_rid, tex_size.x, tex_size.y, int(result[particle_index_to_print].x), int(result[particle_index_to_print].y)))
		
	return result


func unpack_move_std430(data: PackedFloat32Array) -> Array:
	"""
	Turned into debug buffer
	"""
	var result :=  []
	var stride := 4
	var val := Vector3(data[0], data[1], data[2])
	#if val.x < 50:
	#	print(val)
	result.append(val)
	print(result)
	return result

func unpack_particles_positions(data: PackedFloat32Array) -> Array:
	var particles_positions := []
	var stride := 24  # floats per particle

	for i in range(0, data.size(), stride):
		var pos := Vector3(data[i + 0], data[i + 1], data[i + 2])

		particles_positions.append(pos)
		
	return particles_positions


func debug_particles_std430(data: PackedFloat32Array):
	var particles := []
	var stride := 24  # floats per particle
	var count = 0

	for i in range(0, data.size(), stride):
		var pos := Vector3(data[i + 0], data[i + 1], data[i + 2])
		var vel := Vector3(data[i + 8], data[i + 9], data[i + 10])
		var c1 := Vector3(data[i + 12], data[i + 13], data[i + 14])
		var c2 := Vector3(data[i + 16], data[i + 17], data[i + 18])
		var c3 := Vector3(data[i + 20], data[i + 21], data[i + 22])

		var p := {
			"x": pos,
			"v": vel,
			"c1": c1,
			"c2": c2,
			"c3": c3,
		}
		if p.v.x == NAN:
			print("Vel x is nan!")
		if p.v.y == NAN:
			print("Vel y is nan!")
		if p.v.z == NAN:
			print("Vel z is nan!")
		particles.append(p)
	print("---------\n", particles[50].x)
	print("c1: ", particles[50].c1)
	print("c2: ", particles[50].c2)
	print("c3: ", particles[50].c3)
	print("\n---------")
	print(visual_particles[5].x)




func unpack_particles_std430(data: PackedFloat32Array) -> Array:
	var particles : PackedVector3Array = []
	var radii = []
	var stride := 24  # floats per particle
	var count = 0
	for v_i in len(visual_particles):
		var i = visual_particles[v_i].index * stride
	#for p_i in all_num_particles:
	#	var i = p_i * stride

		var pos := Vector3(data[i + 0], data[i + 1], data[i + 2])
		particles.append(pos)
		if (pos.z <= 1):
			count += 1
	var gpup : GPUParticles3D = $GPUParticles3D
	gpup.process_material.set_shader_parameter("particles", particles)

	return particles


func unpack_cells_std430(packed: PackedFloat32Array):
	var cells := []
	var stride := 8  # 4 floats for vec3+pad, 4 more for mass+pad
	var count = 0
	for i in range(0, packed.size(), stride):
		var c_v = Vector3(packed[i + 0], packed[i + 1], packed[i + 2])
		var c_mass = packed[i + 4]
		grid[count].v = c_v
		grid[count].mass = c_mass
		count += 1



func float32_array_to_vec3_array(padded: PackedFloat32Array) -> Array:
	var result := []
	for i in range(0, padded.size(), 4):
		result.append(Vector3(padded[i], padded[i + 1], padded[i + 2]))
	return result
