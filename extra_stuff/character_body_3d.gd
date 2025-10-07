extends CharacterBody3D

@export var speed: float = 5.0
@export var gravity: float = -9.8

var jump_speed = 5
var mouse_sensitivity = 0.008

func _physics_process(delta):
	#velocity.y += -gravity * delta
	var input = Input.get_vector("left", "right", "forward", "back")
	var movement_dir = transform.basis * Vector3(input.x, 0, input.y)
	velocity.x = movement_dir.x * speed
	velocity.z = movement_dir.z * speed
	var input2 = Input.get_axis("down", "up")
	velocity.y =  input2 * speed /  2


	move_and_slide()


func _input(event):
	if event is InputEventMouseMotion: #and Input.mouse_mode == Input.MOUSE_MODE_CAPTURED:
		#rotate_y(-event.relative.x * mouse_sensitivity)
		#$Camera3D.rotate_x(-event.relative.y * mouse_sensitivity)
		#$Camera3D.rotation.x = clampf($Camera3D.rotation.x, -deg_to_rad(70), deg_to_rad(70))
		pass
