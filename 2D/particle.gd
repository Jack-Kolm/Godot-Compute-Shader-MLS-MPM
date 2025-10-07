extends Node2D
class_name Particle

var x : Vector2 = Vector2.ZERO
var v : Vector2 = Vector2.ZERO
var mass : float = 0.0
var C : Transform2D = Transform2D(Vector2.ZERO, Vector2.ZERO, Vector2.ZERO)
var s : Sprite2D
var volume_0 : float = 0.0
func _ready():
	s = Sprite2D.new()
	s.texture = ResourceLoader.load("res://square.png")
	#var scene = ResourceLoader.load("res://2D/PBody.tscn")
	#var s = scene.instantiate()
	s.scale *= 0.001
	self.add_child(s)

	
func _physics_process(delta):
	self.position.x = x.x
	self.position.y = x.y
	#s.position = self.position
	#self.rotation.y += 0.01 * delta
