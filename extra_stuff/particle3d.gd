extends Node3D
class_name Particle3D

var x : Vector3 = Vector3.ZERO
var v : Vector3 = Vector3.ZERO
var mass : float = 0.0
var C : Transform3D = Transform3D(Vector3.ZERO, Vector3.ZERO, Vector3.ZERO, Vector3.ZERO)

var volume_0 : float = 0.0
var index = -1
var cube_model = preload("res://extra_stuff/cube.glb")
var cube_scene = preload("res://extra_stuff/cube_scene.tscn")#preload("res://cube.glb")

func _ready():
	#s = Sprite2D.new()
	var cube = cube_model.instantiate()
	self.add_child(cube)
	cube.scale *= 5.0
	
func _physics_process(delta):
	#self.global_position.x = x.x
	#self.global_position.y = x.y
	#self.global_position.z = x.z
	#lerp_own_position(delta)
	set_own_position()
	#self.global_position.x += 1.0 * delta
	#// Scaling vertices by our base size param (configurable in the material) and the mass of the particle
   # var localPosition = v.vertex.xyz * (_Size * data.w);
	#var worldPosition = data.xyz + localPosition;
func set_own_position():
	#self.global_position.x = x.x
	#self.global_position.y = x.y
	#self.global_position.z = x.z
	#var lpos = self.global_position
	#var wpos = lpos + x
	#self.global_position = wpos
	self.position = x #* 4 - Vector3(10, 10, 10)
func lerp_own_position(delta):
	
	self.global_position = lerp(self.global_position, x, delta)

	#s.position = self.position
	#self.rotation.y += 0.01 * delta
