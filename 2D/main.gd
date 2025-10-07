extends Node2D


const floor_height = 50

const grid_res : int = 64
const num_cells :int = grid_res * grid_res

const gravity : float = 0.3

var num_particles : int

var particles : Array = []
var grid : Array = []
var Fs : Array = []
var weights: Array[Vector2] = [
	Vector2.ZERO,
	Vector2.ZERO,
	Vector2.ZERO
]

var manual_dt = 0.2
var iterations : int = int(1.0/manual_dt)

var thread1 : Thread
var thread2 : Thread
var thread3 : Thread
var mutex : Mutex

func _ready():
	init()
	thread1 = Thread.new()
	#thread2 = Thread.new()
	#thread3 = Thread.new()
	#mutex = Mutex.new()
	#thread.start(_thread_function.bind("Wafflecopter"))
	thread1.start(loop.bind(manual_dt))
	#thread2.start(loop.bind(1.0))

func _physics_process(delta):
	var dt = delta * 50
	#print(dt)
	#iterations = int(1.0/dt)
	#thread1.start(p2g.bind(dt))

	#loop(dt)


func init():
	"""
	// 1.  initialise your grid - fill your grid array with (grid_res * grid_res) cells.
	
	// 2. create a bunch of particles. set their positions somewhere in your simulation domain.
	// initialise their deformation gradients to the identity matrix, as they're in their undeformed state.

	// 3. optionally precompute state variables e.g. particle initial volume, if your model calls for it
	"""
	var spacing = 1.0
	var box = Vector2(32, 32)
	var s = Vector2(grid_res/2.0, grid_res/2.0)
	var temp_positions: Array = []
	var index_i = s.x - box.x / 2.0
	while index_i < s.x + box.x / 2.0:
		var index_j = s.y - box.y / 2.0
		while index_j < s.y + box.y / 2.0:
			var pos = Vector2(index_i, index_j)
			temp_positions.append(pos)
			index_j += spacing
		index_i += spacing
	num_particles = len(temp_positions)
	for i in num_particles:
		var p = Particle.new()  # or preload/instantiate a scene, or create a struct/dictionary
		p.x = temp_positions[i];
		# Random initial velocity (similar to Random.value in Unity)
		p.v = Vector2(randf() - 0.5, -(randf() - 0.5 + 2.75)) * 0.5
		# p.C = 0 # not yet implemented
		p.v = Vector2(0, 0)
		p.mass = 1.0
		#p.C = 0;
		particles.append(p)
		self.add_child(p)
		
		# Fs
		var fi = Transform2D.IDENTITY
		Fs.append(fi)
		
	for i in num_cells:
		var c = Cell.new()
		grid.append(c)
	#p2g(manual_dt)
	#pre_calculation()

func pre_calculation():
	for i in num_particles:
		var p = particles[i]
		#print(p)
		#var stress : Transform2D = Transform2D(Vector2.ZERO, Vector2.ZERO, Vector2.ZERO)
		#var F : Transform2D = Fs[i]
		#var J = Math.determinant_t2d(F)
		#var volume = p.volume_0 * J;
		# quadratic interpolation weights
		var cell_idx: Vector2i = Vector2i(floor(p.x.x), floor(p.x.y)) # I hope this is correct
		#var cell_center: Vector2 = cell_idx + Vector2(0.5, 0.5)
		var cell_diff: Vector2 = (p.x - Vector2(cell_idx)) - Vector2(0.5, 0.5)#p.x - cell_center
		
		#mutex.lock()
		weights[0] = 0.5 * Math.pow2(Vector2(0.5, 0.5) - cell_diff) # I hope this is correct
		weights[1] = Vector2(0.75, 0.75) - Math.pow2(cell_diff)
		weights[2] = 0.5 * Math.pow2(Vector2(0.5, 0.5) + cell_diff)
		#mutex.unlock()
		var density : float = 0.0
		# For all surrounding 9 cells
		for gx : int  in range(3):
			for gy : int  in range(3):
				var weight : float = weights[gx].x * weights[gy].y
			   # int cell_index = ((int)cell_idx.x + (gx - 1)) * grid_res + ((int)cell_idx.y + gy - 1)
				var cell_index : int = int(cell_idx.x + (gx -1)) * grid_res + int(cell_idx.y + (gy - 1))
				density += grid[cell_index].mass * weight
		
		var volume : float = p.mass / density
		p.volume_0 = volume
		particles[i] = p



const elastic_lambda : float = 4.0
const elastic_mu : float = 0.01
func loop(dt):
	while true:
		for i in num_cells:
			#print(grid)
			var cell = grid[i]
			cell.mass = 0
			cell.v = Vector2(0, 0)
			grid[i] = cell # needed?
		#if (!thread1.is_alive() and !thread1.is_started()) or !thread1.is():
		#	print("yes?")
		#if !thread2.is_alive() and !thread2.is_started():
		#thread2.start(grid_velocity_update.bind(dt))
		#if !thread3.is_alive() and !thread3.is_started():
		#thread3.start(g2p.bind(dt))
		p2g(dt)
		p2g_2(dt)
		grid_velocity_update(dt)
		g2p(dt)
	# P2G
func p2g(dt):
	
	for i in num_particles:

		var p = particles[i]
		#print(p)
		"""
		var stress : Transform2D = Transform2D(Vector2.ZERO, Vector2.ZERO, Vector2.ZERO)
		# deformation gradient
		var F : Transform2D = Fs[i]
		var J = Math.determinant_t2d(F)
		
		#mpm course page 46!?!?
		var volume = p.volume_0 * J;
		
		# matrices for neo-hookean model
		var F_T = Math.transpose_t2d(F)
		var F_inv_T = Math.inverse_t2d(F_T)
		var F_minus_F_inv_T = Math.subtract_t2d(F, F_inv_T)
		#print(P_term_0)

		var P_term_0 = Math.float_mul_t2d(elastic_mu, F_minus_F_inv_T) # elastic_mu * (F_minus_F_inv_T)
		var P_term_1 = Math.float_mul_t2d(elastic_lambda * log(J), F_inv_T)
		var P = Math.add_t2d(P_term_0, P_term_1)
		stress = Math.float_mul_t2d((1.0 / J), (P * F_T))
		
		"""		
		#var eq_16_term_0 = Math.float_mul_t2d(-volume * 4, stress) * dt
		
		# quadratic interpolation weights
		var cell_idx: Vector2i = Vector2i(floor(p.x.x), floor(p.x.y)) # I hope this is correct
		#var cell_center: Vector2 = cell_idx + Vector2(0.5, 0.5)
		var cell_diff: Vector2 = (p.x - Vector2(cell_idx)) - Vector2(0.5, 0.5)#p.x - cell_center
		#mutex.lock()
		weights[0] = 0.5 * Math.pow2(Vector2(0.5, 0.5) - cell_diff) # I hope this is correct
		weights[1] = Vector2(0.75, 0.75) - Math.pow2(cell_diff)
		weights[2] = 0.5 * Math.pow2(Vector2(0.5, 0.5) + cell_diff)
		#mutex.unlock()

		# For all surrounding 9 cells
		for gx : int  in range(3):
			for gy : int  in range(3):
				var weight : float = weights[gx].x * weights[gy].y
				var cell_x : Vector2i = Vector2i(cell_idx.x + gx - 1, cell_idx.y + gy - 1)
				var cell_index : int = int(cell_x.x) * grid_res + int(cell_x.y)

				var cell_dist : Vector2 = (Vector2(cell_x) - p.x) + Vector2(0.5, 0.5)
				var Q : Vector2 = p.C * cell_dist
				
				var mass_contrib : float = weight * p.mass
				
				# converting 2d index to 1d
				
				if cell_index < len(grid) and cell_index > 0: # so we're not out of bounds!
					#mutex.lock()
					var cell : Cell = grid[cell_index]
					# scatter mass to grid
					cell.mass += mass_contrib
					cell.v += mass_contrib * (p.v + Q)
					
					# neo-hookean
					#var momentum : Vector2 = (eq_16_term_0 * weight) * cell_dist
					#cell.v += momentum
					
					grid[cell_index] = cell
					

func p2g_2(dt):
	
	for i in num_particles:

		var p = particles[i]
		#print(p)

		
		# quadratic interpolation weights
		var cell_idx: Vector2i = Vector2i(floor(p.x.x), floor(p.x.y)) # I hope this is correct
		#var cell_center: Vector2 = cell_idx + Vector2(0.5, 0.5)
		var cell_diff: Vector2 = (p.x - Vector2(cell_idx)) - Vector2(0.5, 0.5)#p.x - cell_center
		#mutex.lock()
		weights[0] = 0.5 * Math.pow2(Vector2(0.5, 0.5) - cell_diff) # I hope this is correct
		weights[1] = Vector2(0.75, 0.75) - Math.pow2(cell_diff)
		weights[2] = 0.5 * Math.pow2(Vector2(0.5, 0.5) + cell_diff)
		#mutex.unlock()
		var density = 0
		# For all surrounding 9 cells
		for gx : int  in range(3):
			for gy : int  in range(3):
				var weight : float = weights[gx].x * weights[gy].y
				var cell_x : Vector2i = Vector2i(cell_idx.x + gx - 1, cell_idx.y + gy - 1)
				var cell_dist : Vector2 = (Vector2(cell_x) - p.x) + Vector2(0.5, 0.5)
				var cell_index : int = int(cell_x.x) * grid_res + int(cell_x.y)
				density += grid[cell_index].mass * weight
		var rest_density = 4.0
		var dynamic_viscosity = 0.1
		var eos_stiffness = 3.0
		var eos_power = 5.0
		var volume = p.mass / density
		var pressure = max(-0.1, eos_stiffness * (((density/rest_density)**eos_power)-1))
		var stress : Transform2D = Transform2D(Vector2(-pressure, 0.0), Vector2(0.0, -pressure), Vector2.ZERO)
		var dudv = p.C
		var strain = dudv
		var trace = strain.y.x + strain.x.y
		strain.x.y = trace
		strain.y.x = trace
		# trace??
		# viscosity_term = dynamic_viscosity * strain
		# Transform2D doesn't support scalar multiplication directly, so multiply columns
		var viscosity_term := Transform2D(
			strain.x * dynamic_viscosity,
			strain.y * dynamic_viscosity,
			Vector2.ZERO
		)

		# stress += viscosity_term
		stress.x += viscosity_term.x
		stress.y += viscosity_term.y
		var eq_16_term_0 = Math.float_mul_t2d(-volume * 4 * dt, stress)
		for gx : int  in range(3):
			for gy : int  in range(3):
				var weight : float = weights[gx].x * weights[gy].y
				var cell_x : Vector2i = Vector2i(cell_idx.x + gx - 1, cell_idx.y + gy - 1)
				var cell_dist : Vector2 = (Vector2(cell_x) - p.x) + Vector2(0.5, 0.5)
				var cell_index : int = int(cell_x.x) * grid_res + int(cell_x.y)
				var cell = grid[cell_index]
				var momentum = (eq_16_term_0 * weight) * cell_dist
				cell.v += momentum
				grid[cell_index] = cell
func grid_velocity_update(dt):
	for i in range(num_cells):
		var cell = grid[i]
		if cell.mass > 0:
			cell.v = cell.v / cell.mass
			cell.v += dt * Vector2(0.0, -gravity)
			var x : int = i / grid_res
			var y = i % grid_res;
			if x < 2 or x > grid_res-3:
				cell.v.x = 0
			if y < 2 or y > grid_res-3:
				cell.v.y = 0
			#mutex.lock()
			grid[i] = cell
			#mutex.unlock()
		#mutex.lock()
		grid[i] = cell
		#mutex.unlock()
func g2p(dt):
	var max_x = -INF
	var min_x = INF
	var max_y = -INF
	var min_y = INF

	#print("aaa")
	for i in num_particles:
		var p = particles[i] # MUTEX
		# reset particle velocity. we calculate it from scratch each step using the grid
		p.v = Vector2(0, 0)
		# quadratic interpolation weights
		var cell_idx : Vector2i = Vector2i(p.x)
		var cell_diff : Vector2 = (p.x - Vector2(cell_idx)) - Vector2(0.5, 0.5) # MUTEX
		#mutex.lock()

		weights[0] = 0.5 * Math.pow2(Vector2(0.5, 0.5) - cell_diff) # I hope this is correct # MUTEX
		weights[1] = Vector2(0.75, 0.75) - Math.pow2(cell_diff)# MUTEX
		weights[2] = 0.5 * Math.pow2(Vector2(0.5, 0.5) + cell_diff) # MUTEX
		#mutex.unlock()

		# constructing affine per-particle momentum matrix from APIC / MLS-MPM.
		# see APIC paper (https://web.archive.org/web/20190427165435/https://www.math.ucla.edu/~jteran/papers/JSSTS15.pdf), page 6
		# below equation 11 for clarification. this is calculating C = B * (D^-1) for APIC equation 8,
		# where B is calculated in the inner loop at (D^-1) = 4 is a constant when using quadratic interpolation functions
		var B : Transform2D = Transform2D(Vector2.ZERO, Vector2.ZERO, Vector2.ZERO)
		for gx : int in range(3):
			for gy : int in range(3):
				var weight = weights[gx].x * weights[gy].y
				var cell_x : Vector2i = Vector2i(cell_idx.x + gx - 1, cell_idx.y + gy - 1)
				var cell_index : int = int(cell_x.x) * grid_res + int(cell_x.y)
				var dist : Vector2 = (Vector2(cell_x) - p.x) + Vector2(0.5, 0.5)
				#if cell_index < len(grid) and cell_index > 0:

				var weighted_velocity = grid[cell_index].v * weight #+ Vector2(weight, weight)
				
				# APIC paper equation 10, constructing inner term for B
				var term : Transform2D = Transform2D(weighted_velocity * dist.x, weighted_velocity * dist.y, Vector2.ZERO)
				#var terma = Transform2D(
				#	Vector2(weight * grid[cell_index].v.x * dist.x, weight * grid[cell_index].v.x * dist.y),
				#	Vector2(weight * grid[cell_index].v.y * dist.x, weight * grid[cell_index].v.y * dist.y),
				#	Vector2(0, 0)
				#)
				B = Math.add_t2d(B, term)
				p.v += weighted_velocity
		p.C = B * 4
		
		
		# This is test
		
		#var clamped_v = Vector2(
		#	clamp(p.v.x, -0.5, 0.5),
		#	clamp(p.v.y, -0.5, 0.5)
		#)
		#p.v = clamped_v
		# advect particles
		#var dv = p.v * dt
		#var new_x = p.x + dv
		#print(new_x.y)
		#if new_x.y > floor_height:
		#	print(new_x.y)
		#	p.x.x += dv.x
		#else:
		p.x += p.v * dt;
		#print(p.v)
		# safety clamp to ensure particles don't exit simulation domain # this doesn't seem to transfer that well into godot
		#p.x = clamp(p.x, Vector2(1, 1), Vector2(grid_res - 2, grid_res - 2))
		#print(grid_res - 2)
		var clamped_x = Vector2(
			clamp(p.x.x, 1.0, grid_res - 2),
			clamp(p.x.y, 1.0, grid_res - 2)
		)
		p.x = clamped_x
		#if Input.is_mouse_button_pressed(MOUSE_BUTTON_LEFT):
		#	#print("MOUSE")
		#	var mouse_pos = get_viewport().get_mouse_position()
		#	var camera_2d = $Camera2D
		#	mouse_pos = camera_2d.get_screen_to_world(mouse_pos)
		#	var mouse_radius = 1.0
		#	var dist = p.x - mouse_pos
		#	if dist.length_squared() < mouse_radius * mouse_radius:
		#		var norm_factor = dist.length() / mouse_radius
		#		norm_factor = pow(sqrt(norm_factor), 8)
		#		var force = dist.normalized() * norm_factor * 0.5
		#		p.v += force
		#mutex.lock()
		#mutex.unlock()
		#p.global_position.x = p.x.x
		#p.global_position.y = p.x.y
		
		#var Fp_new = Transform2D.IDENTITY
		#Fp_new = Math.add_t2d(Fp_new, Math.float_mul_t2d(dt, p.C))
		#Fs[i] = Fp_new * Fs[i]
		var x_n = p.x + p.v
		var wall_min = 3
		var wall_max = grid_res - 4
		if (x_n.x < wall_min): p.v.x += wall_min - x_n.x
		if (x_n.x > wall_max): p.v.x += wall_max - x_n.x
		if (x_n.y < wall_min): p.v.y += wall_min - x_n.y
		if (x_n.y > wall_max): p.v.y += wall_max - x_n.x
		particles[i] = p

		#if p.x.x < min_x: min_x = p.x.x
		#if p.x.x > max_x: max_x = p.x.x
		#if p.x.y < min_y: min_y = p.x.y
		#if p.x.y > max_y: max_y = p.x.y
	#print(min_y, " ", max_y) # x: 1.58, 63.2
	#print(particles[8].x, particles[8].position, particles[8].global_position)
