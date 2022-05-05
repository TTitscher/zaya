import taichi as ti

ti.init(kernel_profiler=True, arch=ti.cuda)

RES_X = 1000
RES_Y = RES_X // 2

MAX_STEPS = 100
MAX_DIST = 1000.0
EPS = 0.01

pixels = ti.Vector.field(3, dtype=ti.f32, shape=(RES_X, RES_Y))

camera_pos = ti.Vector.field(3, dtype=float, shape=())
camera_rot = ti.Vector.field(2, dtype=float, shape=())
camera_pos[None].y = 1.


@ti.func
def get_distance(p):
    sphere= ti.Vector((0.0, 1.0, 6.0, 1.))
    sphere2 = ti.Vector((1.0, 2.0, 6.0, 1.))
    sphere_distance = (p - sphere.xyz).norm() - sphere.w
    sphere_distance2 = (p - sphere2.xyz).norm() - sphere.w

    op = ti.min(sphere_distance, sphere_distance2)

    # plane at (0,0,0) pointing in y direction
    plane_distance = p.y
    return ti.min(op, plane_distance)


@ti.func
def ray_march(ray_origin, ray_direction):
    """
    Finds the distance `d` in
        closest_surface_coordinate = `ray_origin` + d * `ray_direction`
    via ray marching. See wiki or so.
    """
    d = 0.0
    for i in range(MAX_STEPS):
        p = ray_origin + d * ray_direction
        distance = get_distance(p)
        d += distance

        if d > MAX_DIST or distance < EPS:
            break

    return d


@ti.func
def get_normal(p):
    d = get_distance(p)

    dx = get_distance(p - ti.Vector((EPS, 0.0, 0.0)))
    dy = get_distance(p - ti.Vector((0.0, EPS, 0.0)))
    dz = get_distance(p - ti.Vector((0.0, 0.0, EPS)))

    return (d - ti.Vector((dx, dy, dz))).normalized()


@ti.func
def get_light(p, t, light_pos):

    """
    Lighing is brightest when the light vector (current point to light source)
    is colinear with the surface normal.
    """
    light_vector = (light_pos - p).normalized()
    surface_normal = get_normal(p)
    diffuse_light = surface_normal.dot(light_vector)

    """
    Shadow is surprisingly easy: If the distance from the point to the light
    source is smaller than the distance to the light source, there is something
    blocking it.

    However, p already intersects with the scene, so the distance will be 0.
    So we have to offset it a little bit (here in normal direction) to avoid
    this first intersection.
    """
    d = ray_march(p + surface_normal * 2.0 * EPS, light_vector)
    d_ref = (p - light_pos).norm()

    if d < d_ref:
        diffuse_light *= 0.5

    return ti.max(diffuse_light, 0.0)


@ti.kernel
def paint(t: float):
    for i, j in pixels:  # Parallelized over all pixels

        x = (i - 0.5 * RES_X) / RES_Y  # / RES_Y is not a typo. Makes dx == dy
        y = (j - 0.5 * RES_Y) / RES_Y

        ray_origin = camera_pos[None]
        ray_direction = ti.Vector((x, y, 1.0)).normalized()

        Rx = ti.Matrix.rotation2d(camera_rot[None].x)
        Ry = ti.Matrix.rotation2d(camera_rot[None].y)
        ray_direction.zy = Ry @ ray_direction.zy
        ray_direction.zx = Rx @ ray_direction.zx

        # find the intersection point `p` in in distance `d` along the ray
        d = ray_march(ray_origin, ray_direction)
        p = ray_origin + ray_direction * d

        light_pos = ti.Vector((0.0, 20.0, 6.0))
        light_pos.x += ti.sin(t) * 30
        light_pos.z += ti.cos(t) * 30

        light_pos2 = ti.Vector((0.0, 20.0, 6.0))
        light_pos2.x -= ti.sin(-t) * 30
        light_pos2.z -= ti.cos(t) * 30

        diffuse_light = get_light(p, t, light_pos)
        diffuse_light2 = get_light(p, t, light_pos2)

        pixels[i, j] = diffuse_light, 0.0, diffuse_light2



gui = ti.GUI("test stuff", res=(RES_X, RES_Y))

t = 0
dt = 0.02
CAM_SPEED = 0.1
CAM_ROT_SPEED= 0.1

while True:
    if gui.get_event(ti.GUI.PRESS):
        if gui.event.key in [ti.GUI.ESCAPE, ti.GUI.EXIT]:
            break
    if gui.is_pressed("w"):
        camera_pos[None].z += CAM_SPEED
    if gui.is_pressed("s"):
        camera_pos[None].z -= CAM_SPEED
    if gui.is_pressed("a"):
        camera_pos[None].x -= CAM_SPEED
    if gui.is_pressed("d"):
        camera_pos[None].x += CAM_SPEED
   
    if gui.is_pressed(ti.GUI.UP):
        camera_rot[None].y += CAM_ROT_SPEED
    if gui.is_pressed(ti.GUI.DOWN):
        camera_rot[None].y -= CAM_ROT_SPEED
    if gui.is_pressed(ti.GUI.LEFT):
        camera_rot[None].x -= CAM_ROT_SPEED
    if gui.is_pressed(ti.GUI.RIGHT):
        camera_rot[None].x += CAM_ROT_SPEED

    paint(t)
    gui.set_image(pixels)
    gui.show()
    t += dt
