extends Node3D


const floor_height = 50

const grid_res : int = 16
const num_cells :int = grid_res * grid_res * grid_res

const gravity : float = -0.3

var all_num_particles : int

var all_particles : Array = []
var grid : Array = []
var Fs : Array = []


var manual_dt = 0.2
var iterations : int = int(1.0/manual_dt)

var thread1 : Thread
var thread2 : Thread
var thread3 : Thread
var mutex : Mutex

func _ready():
	await init()
	#shader_prep()
	manual_dt = 0.2
	thread1 = Thread.new()
	thread1.start(thread_loop.bind(manual_dt, all_particles, all_num_particles))

func _physics_process(delta):
	#loop_iteration()
	pass



var debugbuffer = []

func init():

	var spacing = 1.0
	var box = Vector3(8, 8, 8)
	var s = Vector3(grid_res/2.0, grid_res/2.0, grid_res/2.0)
	
	var temp_positions: Array = []
	var index_i = s.x - box.x / 2.0
	while index_i < s.x + box.x / 2.0:
		var index_j = s.y - box.y / 2.0
		while index_j < s.y + box.y / 2.0:
			var index_k = s.z - box.z / 2.0
			while index_k < s.z + box.z / 2.0:
				var pos = Vector3(index_i, index_j, index_k)
				temp_positions.append(pos)
				index_k += spacing
			index_j += spacing
		index_i += spacing
	all_num_particles = len(temp_positions)
	for i in all_num_particles:
		var p = Particle3D.new()  # or preload/instantiate a scene, or create a struct/dictionary
		p.x = temp_positions[i];
		#p.v = Vector3(randf() - 0.5, (randf() - 0.5 + 2.75), randf() - 0.5) * 0.5 # random initial velocity 
		p.v = Vector3(0.0, 0.0, 0.0)
		# p.C = 0 # not yet implemented
		p.mass = 1.0
		all_particles.append(p)
		self.add_child(p)
		
		# Fs
		var fi = Transform3D.IDENTITY
		Fs.append(fi)
		debugbuffer.append(Vector3.ZERO)
	for i in num_cells:
		var c = Cell3D.new()
		grid.append(c)
	print(all_num_particles)
	#p2g(manual_dt, all_particles, all_num_particles)
	#pre_calculation(all_particles, all_num_particles)
	return

func pre_calculation(particles, num_particles):
	for i in num_particles:
		var p = particles[i]

		var cell_idx: Vector3i = Vector3i(floor(p.x.x), floor(p.x.y), floor(p.x.z)) # I hope this is correct
		var cell_diff: Vector3 = (p.x - Vector3(cell_idx)) - Vector3(0.5, 0.5, 0.5)#p.x - cell_center
		var weights: Array[Vector3] = [
			Vector3.ZERO,
			Vector3.ZERO,
			Vector3.ZERO
		]
		weights[0] = 0.5 * Math.pow3(Vector3(0.5, 0.5, 0.5) - cell_diff) # I hope this is correct
		weights[1] = Vector3(0.75, 0.75, 0.75) - Math.pow3(cell_diff)
		weights[2] = 0.5 * Math.pow3(Vector3(0.5, 0.5, 0.5) + cell_diff)
		var density : float = 0.0
		
		# for all surrounding 9 cells
		for gx : int  in range(3):
			for gy : int  in range(3):
				for gz : int  in range(3):
					var weight : float = weights[gx].x * weights[gy].y * weights[gz].z
					var cell_index : int = int(cell_idx.x + (gx -1)) * grid_res + int(cell_idx.y + (gy - 1)) + int(cell_idx.z + (gz -1)) * grid_res * grid_res
					density += grid[cell_index].mass * weight
		var volume : float = p.mass / density
		p.volume_0 = volume
		particles[i] = p


func thread_loop(dt, particles, num_particles):
	while true:
		loop(dt, particles, num_particles)
const elastic_lambda : float = 0.01
const elastic_mu : float = 10.0
var debug = Vector3(-1.0, -1.0, -1.0)

func loop(dt, particles, num_particles):
	for i in num_cells:
		var cell = grid[i]
		cell.mass = 0.0
		cell.v = Vector3(0, 0, 0)
		grid[i] = cell # needed?
	p2g(dt, particles, num_particles)
	#p2g_2(dt, particles, num_particles)
	grid_velocity_update(dt, particles, num_particles)
	g2p(dt, particles, num_particles)
	#print(debug)
	print(debugbuffer[0])

# P2G
func p2g(dt, particles, num_particles):
	
	for i in num_particles:

		var p = particles[i]
		var cell_idx: Vector3i = Vector3i(floor(p.x.x), floor(p.x.y), floor(p.x.z)) # I hope this is correct
		var cell_diff: Vector3 = p.x - (Vector3(cell_idx) + Vector3(0.5, 0.5, 0.5))
		# quadratic interpolation weights
		var weights: Array[Vector3] = [
			Vector3.ZERO,
			Vector3.ZERO,
			Vector3.ZERO
		]
		weights[0] = 0.5 * Math.pow3(Vector3(0.5, 0.5, 0.5) - cell_diff) # I also hope this is correct
		weights[1] = Vector3(0.75, 0.75, 0.75) - Math.pow3(cell_diff)
		weights[2] = 0.5 * Math.pow3(Vector3(0.5, 0.5, 0.5) + cell_diff)

		#print("w", weights[0])
		# For all surrounding 9 cells
		for gx : int  in range(3):
			for gy : int  in range(3):
				for gz : int  in range(3):

					var weight : float = weights[gx].x * weights[gy].y * weights[gz].z
					var cell_x : Vector3i = Vector3i(cell_idx.x + gx - 1, cell_idx.y + gy - 1, cell_idx.z + gz - 1)
					var cell_dist : Vector3 = (Vector3(cell_x) - p.x) + Vector3(0.5, 0.5, 0.5)
					var Q : Vector3 = p.C * cell_dist
					debugbuffer[i] = weight

					var mass_contrib : float = weight * p.mass
					
					# converting 3d index to 1d
					var cell_index : int = int(cell_x.x) * grid_res * grid_res + int(cell_x.y)  * grid_res + int(cell_x.z)
					var cell : Cell3D = grid[cell_index]
					# scatter mass and momentum to grid
					cell.mass += mass_contrib
					#Q = Vector3(0, 0, 0)
					cell.v += mass_contrib * (p.v + Q)
					
					
					grid[cell_index] = cell
						
func p2g_2(dt, particles, num_particles):
	
	for i in num_particles:

		var p = particles[i]

		
		# quadratic interpolation weights
		var cell_idx: Vector3i = Vector3i(floor(p.x.x), floor(p.x.y), floor(p.x.z)) # I hope this is correct
		#var cell_diff: Vector3 = (p.x - Vector3(cell_idx)) - Vector3(0.5, 0.5, 0.5)#p.x - cell_center
		var cell_diff: Vector3 = p.x - (Vector3(cell_idx) + Vector3(0.5, 0.5, 0.5))

		var weights: Array[Vector3] = [
			Vector3.ZERO,
			Vector3.ZERO,
			Vector3.ZERO
		]
		weights[0] = 0.5 * Math.pow3(Vector3(0.5, 0.5, 0.5) - cell_diff) # I hope this is correct
		weights[1] = Vector3(0.75, 0.75, 0.75) - Math.pow3(cell_diff)
		weights[2] = 0.5 * Math.pow3(Vector3(0.5, 0.5, 0.5) + cell_diff)

		var density = 0
		# For all surrounding 9 cells
		for gx : int  in range(3):
			for gy : int  in range(3):
				for gz : int  in range(3):
					var weight : float = weights[gx].x * weights[gy].y
					var cell_x : Vector3i = Vector3i(cell_idx.x + gx - 1, cell_idx.y + gy - 1, cell_idx.z + gz - 1)
					var cell_index : int = int(cell_x.x) * grid_res * grid_res + int(cell_x.y)  * grid_res + int(cell_x.z)
					var cell_dist : Vector3 = (Vector3(cell_x) - p.x) + Vector3(0.5, 0.5, 0.5)
					density += grid[cell_index].mass * weight
		var rest_density = 4.0
		var dynamic_viscosity = 0.1
		var eos_stiffness = 3.0
		var eos_power = 5.0
		var volume = p.mass / density
		var pressure = max(-0.0, eos_stiffness * (((density/rest_density)**eos_power)-1))
		var stress : Transform3D = Transform3D(Vector3(-pressure, 0.0, 0.0), Vector3(0.0, -pressure, 0.0), Vector3(0.0, 0.0, -pressure), Vector3.ZERO)
		var dudv = p.C
		var strain = Math.add_t3d(dudv, Math.transpose_t3d(dudv))
		stress = Math.add_t3d(stress, strain * dynamic_viscosity)

		var eq_16_term_0 =  stress * (-volume * 4 * dt) # ash ash # 

		for gx : int  in range(3):
			for gy : int  in range(3):
				for gz : int  in range(3):

					var weight = weights[gx].x * weights[gy].y * weights[gz].z
					var cell_x : Vector3i = Vector3i(cell_idx.x + gx - 1, cell_idx.y + gy - 1, cell_idx.z + gz - 1)
					var cell_index : int = int(cell_x.x) * grid_res * grid_res + int(cell_x.y)  * grid_res + int(cell_x.z)
					var cell_dist : Vector3 = (Vector3(cell_x) - p.x) + Vector3(0.5, 0.5, 0.5)
					var cell = grid[cell_index]
					var momentum = eq_16_term_0 * weight * cell_dist
					cell.v += momentum
					grid[cell_index] = cell

func grid_velocity_update(dt, particles, num_particles):
	for i in range(num_cells):
		var cell = grid[i]
		if cell.mass > 0:
			cell.v = cell.v / cell.mass
			cell.v += dt * Vector3(0.0, gravity, 0.0)
			var x = i / (grid_res * grid_res)#(i / (grid_res * grid_res)) % grid_res
			var y = (i / grid_res) % grid_res#(i / grid_res) % grid_res
			var z = i % grid_res
			if x < 2 or x > grid_res-3:
				cell.v.x = 0
			if y < 2 or y > grid_res-3:
				cell.v.y = 0
			if z < 2 or z > grid_res-3:
				cell.v.z = 0
			grid[i] = cell
		grid[i] = cell

func g2p(dt, particles, num_particles):
	for i in num_particles:
		var p = particles[i] # MUTEX
		# reset particle velocity. we calculate it from scratch each iteration
		p.v = Vector3(0, 0, 0)
		var cell_idx : Vector3i = Vector3i(p.x)
		var cell_diff: Vector3 = p.x - (Vector3(cell_idx) + Vector3(0.5, 0.5, 0.5))
		# quadratic interpolation weights
		var weights: Array[Vector3] = [
			Vector3.ZERO,
			Vector3.ZERO,
			Vector3.ZERO
		]
		weights[0] = 0.5 * Math.pow3(Vector3(0.5, 0.5, 0.5) - cell_diff) # I hope this is correct # MUTEX
		weights[1] = Vector3(0.75, 0.75, 0.75) - Math.pow3(cell_diff)# MUTEX
		weights[2] = 0.5 * Math.pow3(Vector3(0.5, 0.5, 0.5) + cell_diff) # MUTEX

		# constructing affine matrix from APIC / MLS-MPM.
		var B : Transform3D = Transform3D(Vector3.ZERO, Vector3.ZERO, Vector3.ZERO, Vector3.ZERO)
		for gx : int in range(3):
			for gy : int in range(3):
				for gz : int in range(3):
					var weight = weights[gx].x * weights[gy].y * weights[gz].z
					var cell_x : Vector3i = Vector3i(cell_idx.x + gx - 1, cell_idx.y + gy - 1, cell_idx.z + gz - 1)
					var cell_index : int = int(cell_x.x) * grid_res * grid_res + int(cell_x.y)  * grid_res + int(cell_x.z)
					var cell_dist : Vector3 = (Vector3(cell_x) - p.x) + Vector3(0.5, 0.5, 0.5)
					var weighted_velocity = grid[cell_index].v * weight #+ Vector2(weight, weight)
					var term : Transform3D = Transform3D(weighted_velocity * cell_dist.x, weighted_velocity * cell_dist.y, weighted_velocity * cell_dist.z, Vector3.ZERO)
					B = Math.add_t3d(B, term)
					p.v += weighted_velocity
		p.C = B * 4
		
		var clamped_v = Vector3(
			clamp(p.v.x, -1.0, 1.0),
			clamp(p.v.y, -1.0, 1.0),
			clamp(p.v.z, -1.0, 1.0)
		)
		var dv = p.v * dt
		var new_x = p.x + dv
		p.x = new_x
		
		# safety clamp to ensure particles don't exit simulation domain # this doesn't seem to transfer that well into godot?
		var clamped_x = Vector3(
			clamp(p.x.x, 1.0, grid_res - 2),
			clamp(p.x.y, 1.0, grid_res - 2),
			clamp(p.x.z, 1.0, grid_res - 2)
		)
		p.x = clamped_x
		var x_n = p.x + p.v
		var wall_min = 3
		var wall_max = grid_res - 4
		if (x_n.x < wall_min): p.v.x += wall_min - x_n.x
		if (x_n.x > wall_max): p.v.x += wall_max - x_n.x
		if (x_n.y < wall_min): p.v.y += wall_min - x_n.y
		if (x_n.y > wall_max): p.v.y += wall_max - x_n.x
		if (x_n.z < wall_min): p.v.z += wall_min - x_n.z
		if (x_n.z > wall_max): p.v.z += wall_max - x_n.z
		
		particles[i] = p
