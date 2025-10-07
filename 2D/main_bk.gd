extends Node3D



const grid_res : int = 64
const num_cells :int = grid_res * grid_res

const gravity : float = -0.05

var num_particles : int

var particles : Array = []
var grid

var weights: Array[Vector2] = [
	Vector2.ZERO,
	Vector2.ZERO,
	Vector2.ZERO
]
var dt : float = 1.0
var iterations : int = int(1.0/dt)

func _physics_process(delta):
	dt = delta
	iterations = int(1.0/dt)
	loop(dt)

func init():
	"""
	// 1.  initialise your grid - fill your grid array with (grid_res * grid_res) cells.
	
	// 2. create a bunch of particles. set their positions somewhere in your simulation domain.
	// initialise their deformation gradients to the identity matrix, as they're in their undeformed state.

	// 3. optionally precompute state variables e.g. particle initial volume, if your model calls for it
	"""
	var spacing = 1.0
	var box = Vector2(16, 16)
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
		p.v = Vector2(randf() - 0.5, randf() - 0.5 + 2.75) * 0.5
		# p.C = 0 # not yet implemented
		p.mass = 1.0
		#p.C = 0;
		particles.append(p)
		self.add_child(p)
	for i in num_cells:
		var c = Cell.new()
		grid.append(c)
		


func loop(dt):
	for i in num_cells:
		var cell = grid[i]
		cell.mass = 0
		cell.v = 0
		grid[i] = cell # needed?
	
	# P2G
	for i in num_particles:
		var p = particles[i]
		# quadratic interpolation weights
		var cell_idx: Vector2i = Vector2i(floor(p.x.x), floor(p.x.y)) # I hope this is correct
		#var cell_center: Vector2 = cell_idx + Vector2(0.5, 0.5)
		var cell_diff: Vector2 = (p.x - cell_idx) - Vector2(0.5, 0.5)#p.x - cell_center
		weights[0] = 0.5 * Math.pow2(Vector2(0.5, 0.5) - cell_diff) # I hope this is correct
		weights[1] = Vector2(0.75, 0.75) - Math.pow2(cell_diff)
		weights[2] = 0.5 * Math.pow2(Vector2(0.5, 0.5) + cell_diff)
		# For all surrounding 9 cells
		for gx : int  in range(3):
			for gy : int  in range(3):
				var weight : float = weights[gx].x * weights[gy].y
				var cell_x : Vector2i = Vector2i(cell_idx.x + gx -1, cell_idx.y + gy -1)
				var cell_dist : Vector2 = (cell_x - p.x) + Vector2(0.5, 0.5)
				var Q : float = p.C * cell_dist
				
				var mass_contrib : float = weight * p.mass
				
				# converting 2d index to 1d
				var cell_index : int = int(cell_x.x) * grid_res + int(cell_x.y)
				var cell : Cell = grid[cell_index]
				
				# scatter mass to grid
				cell.mass += mass_contrib
				cell.v += mass_contrib * (p.v + Q)
				grid[cell_index] = cell
	
	# grid velocity update
	for i in range(num_cells):
		var cell = grid[i]
		if cell.mass > 0:
			cell.v = cell.v / cell.mass
			cell.v += dt * Vector2(0.0, gravity)
			var x : int = i / grid_res
			var y = i % grid_res;
			if x < 2 or x > grid_res-3:
				cell.v.x = 0
			if y < 2 or y > grid_res-3:
				cell.v.y = 0
			grid[i] = cell
		grid[i] = cell
	
	# G2P
	for i in num_particles:
		var p = particles[i]
		# reset particle velocity. we calculate it from scratch each step using the grid
		p.v = 0
		# quadratic interpolation weights
		var cell_idx : Vector2i = Vector2i(p.x)
		var cell_diff : Vector2 = (p.x - cell_idx) - Vector2(0.5, 0.5)
		weights[0] = 0.5 * Math.pow2(Vector2(0.5, 0.5) - cell_diff) # I hope this is correct
		weights[1] = Vector2(0.75, 0.75) - Math.pow2(cell_diff)
		weights[2] = 0.5 * Math.pow2(Vector2(0.5, 0.5) + cell_diff)
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
				var dist : Vector2 = (cell_x - p.x) + Vector2(0.5, 0.5)
				var weighted_velocity = grid[cell_index].v + weight
				
				# APIC paper equation 10, constructing inner term for B
				var term = Transform2D(weighted_velocity * dist.x, weighted_velocity * dist.y, Vector2.ZERO)
				B += term
				p.v += weighted_velocity
		p.C = B * 4
		
		# advect particles
		p.x += p.v * dt;
		
		# safety clamp to ensure particles don't exit simulation domain
		p.x = clamp(p.x, 1, grid_res - 2)
		particles[i] = p


	# Assuming you have a Particle class with properties like:
# F: Transform2D (or 2x2 matrix), C: Transform2D, v: Vector2, x: Vector2, etc.
var elastic_lambda = 1.0
var elastic_mu = 1.0
func loop_(dt):
	for p in particles:
		var F = p.F  # deformation gradient as Transform2D

		# MPM course page 13 - "Kinematics Theory"
		var J = determinant_2x2(F)
		var volume = p.volume_0 * J

		# Neo-Hookean terms (eq. 48)
		var F_T = F.transposed()
		var F_inv_T = inverse_2x2(F_T)
		var F_minus_F_inv_T = subtract_2x2(F, F_inv_T)

		var P_term_0 = elastic_mu * F_minus_F_inv_T
		var P_term_1 = elastic_lambda * log(J) * F_inv_T
		var P = add_2x2(P_term_0, P_term_1)

		# Cauchy stress = (1 / det(F)) * P * F^T
		var stress = (1.0 / J) * mul_2x2(P, F_T)

		# Eq. 16 term from MLS-MPM
		var eq_16_term_0 = -volume * 4.0 * dt * stress

		# Loop over neighboring grid cells
		for cell in p.neighbourhood:
			var cell_dist = (cell.x - p.x) + Vector2(0.5, 0.5)
			var Q = mul_2x2(p.C, cell_dist)
			var weight = bspline_weight(p.x, cell.x)
			# Equation (172)
			var weighted_mass = weight * p.mass
			cell.mass += weighted_mass
			cell.v += weighted_mass * (p.v + Q)

			# Fused momentum + force (MLS-MPM)
			var momentum = mul_2x2(eq_16_term_0, cell_dist) * cell.weight
			cell.v += momentum
