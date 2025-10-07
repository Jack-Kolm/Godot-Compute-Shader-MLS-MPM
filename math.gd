extends Node


func pow2(v: Vector2) -> Vector2:
	return Vector2(v.x * v.x, v.y * v.y)

func pow3(v: Vector3) -> Vector3:
	return Vector3(v.x * v.x, v.y * v.y, v.z * v.z)

## Mat2

func add_t2d(a: Transform2D, b: Transform2D) -> Transform2D:
	return Transform2D(
		a.x + b.x,
		a.y + b.y,
		a.origin + b.origin
	)

func subtract_t2d(a: Transform2D, b: Transform2D) -> Transform2D:
	var x = a.x - b.x
	var y = a.y - b.y
	var origin = a.origin - b.origin
	return Transform2D(x, y, origin)


func determinant_t2d(F: Transform2D) -> float:
	return F.x.x * F.y.y - F.y.x * F.x.y

# Transpose for Transform2D
func transpose_t2d(t: Transform2D) -> Transform2D:
	var x = Vector2(t.x.x, t.y.x)
	var y = Vector2(t.x.y, t.y.y)
	var origin = Vector2.ZERO  # Reset origin when transposing
	return Transform2D(x, y, origin)

# Inverse for Transform2D (more precise than using .affine_inverse())
func inverse_t2d(t: Transform2D) -> Transform2D:
	var det = t.x.x * t.y.y - t.x.y * t.y.x
	if abs(det) < 1e-8:
		push_error("Transform2D is not invertible (determinant is zero)")
		return Transform2D()
	
	var inv_det = 1.0 / det
	var x = Vector2( t.y.y * inv_det, -t.x.y * inv_det)
	var y = Vector2(-t.y.x * inv_det,  t.x.x * inv_det)
	var origin = -Vector2(x.dot(t.origin), y.dot(t.origin))
	
	return Transform2D(x, y, origin)

func float_mul_t2d(fp : float, t : Transform2D):
	return Transform2D(
		t.x * fp,
		t.y * fp,
		t.origin
	)

## 3D matrices

func add_t3d(a: Transform3D, b: Transform3D) -> Transform3D:
	var basis = Basis(
		a.basis.x + b.basis.x,
		a.basis.y + b.basis.y,
		a.basis.z + b.basis.z
	)
	var origin = a.origin + b.origin
	return Transform3D(basis, origin)


func determinant_t3d(F: Transform3D) -> float:
	return F.basis.determinant()


func subtract_t3d(a: Transform3D, b: Transform3D) -> Transform3D:
	var basis = Basis(
		a.basis.x - b.basis.x,
		a.basis.y - b.basis.y,
		a.basis.z - b.basis.z
	)
	var origin = a.origin - b.origin
	return Transform3D(basis, origin)


# Transpose for Transform3D
func transpose_t3d(t: Transform3D) -> Transform3D:
	var basis = Basis(
		Vector3(t.basis.x.x, t.basis.y.x, t.basis.z.x),
		Vector3(t.basis.x.y, t.basis.y.y, t.basis.z.y),
		Vector3(t.basis.x.z, t.basis.y.z, t.basis.z.z)
	)
	# Reset origin when transposing
	return Transform3D(basis, Vector3.ZERO)

# Inverse for Transform3D (manual implementation)
func inverse_t3d(t: Transform3D) -> Transform3D:
	var basis_inv = t.basis.inverse()
	var origin_inv = basis_inv * -t.origin
	return Transform3D(basis_inv, origin_inv)
	

func float_mul_t3d_(fp: float, t: Transform3D) -> Transform3D:
	var scaled_basis = Basis(
		t.basis.x * fp,
		t.basis.y * fp,
		t.basis.z * fp
	)
	return Transform3D(scaled_basis, t.origin)
	
func float_mul_t3d(fp: float, t: Transform3D) -> Transform3D:
	return Transform3D(t.basis * fp, t.origin)
