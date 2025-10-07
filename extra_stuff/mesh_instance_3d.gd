extends MeshInstance3D

@export var motion_radius: float = 0.3
@export var min_radius: float = 0.05
@export var max_radius: float = 0.15
@export var speed: float = 1.0
@export var blend_r: float = 0.2
const BALL_COUNT := 174

var phases_pos: Array = []
var phases_dir: Array = []
var phases_rad: Array = []
var mat: ShaderMaterial

var particle_positions = []

func _ready():
	# Make sure we have a unique material instance
	if get_surface_override_material(0) == null:
		set_surface_override_material(0, mesh.surface_get_material(0).duplicate())
	mat = get_surface_override_material(0)

	# Random phases
	for i in range(BALL_COUNT):
		phases_pos.append(Vector3(randf() * TAU, randf() * TAU, randf() * TAU))
		phases_dir.append(Vector3(randf_range(0.5, 1.5), randf_range(0.5, 1.5), randf_range(0.5, 1.5)))
		phases_rad.append(randf() * TAU)

	mat.set_shader_parameter("blend_r", blend_r)

func _tprocess(delta: float):
	var time = Time.get_ticks_msec() / 1000.0 * speed

	var positions := PackedVector3Array()
	var radii := PackedFloat32Array()
	
	for i in range(BALL_COUNT):
		var phase = phases_pos[i] + phases_dir[i] * time
		var pos = Vector3(sin(phase.x), sin(phase.y * 1.2), sin(phase.z * 0.8)) * motion_radius
		positions.append(pos)

		var r_phase = phases_rad[i] + time
		var final = lerp(min_radius, max_radius, 0.5 + 0.5 * sin(r_phase))
		final = 0.0

		radii.append(final)
	mat.set_shader_parameter("ball_pos", positions)
	mat.set_shader_parameter("radi", radii)

func iteration(dt, i, positions, radii):
	#var time = Time.get_ticks_msec() / 1000.0 * speed

	#var positions := PackedVector3Array()
	#var radii := PackedFloat32Array()
	
	#var phase = phases_pos[i] + phases_dir[i] * time
	#var pos = Vector3(sin(phase.x), sin(phase.y * 1.2), sin(phase.z * 0.8)) * motion_radius
	#positions.append(pos)

	#var r_phase = phases_rad[i] + 5.0
	#radii.append(lerp(min_radius, max_radius, 0.5 + 0.5 * sin(r_phase)))
	#particle_positions[i] = po
	mat.set_shader_parameter("ball_pos", positions)
	#mat.set_shader_parameter("radi", radii)
