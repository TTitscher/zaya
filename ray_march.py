import taichi as ti
import numpy as np

ti.init(arch=ti.gpu)

RES_X = 1000
RES_Y = RES_X // 2

MAX_STEPS = 200
MAX_DIST = 1000.0
EPS = 0.01

pixels = ti.Vector.field(3, dtype=ti.f32, shape=(RES_X, RES_Y))

time = ti.Vector.field(1, dtype=float, shape=())

cam_pos = ti.Vector.field(3, dtype=float, shape=())
cam_rot = ti.Vector.field(2, dtype=float, shape=())
cam_pos[None].y = 3.0


# N = 10
# spheres_np = np.random.random((N, 4)).astype(float)
# spheres = ti.Vector.field(4, dtype=ti.f32, shape=N)
# spheres.from_numpy(spheres_np)


@ti.func
def sd_sphere(p, center, radius):
    return (p - center).norm() - radius


@ti.func
def sd_box(p, center):
    p = ti.abs(p) - center
    return ti.max(p, 0.).norm() + ti.min(ti.max(p.x, max(p.y, p.z)), 0)

@ti.func
def sd_gyroid(p, scale):
    p *= scale
    return ti.sin(p).dot(ti.cos(p.zxy)) / scale


@ti.func
def get_distance(p):
    op = MAX_DIST
    # spheres[:].xyz
    t = time[None].x

    d1 = sd_gyroid(p, 5.)
    d2 = sd_gyroid(p + [0, t, 0], 5.)

    # scale = 4.;
    # d = sd_box((p/scale - [0, 1, 2]), [1,1,1] ) * scale
    d = sd_sphere(p % 10, [5.0, 5.0, 5.0], 4.0 + 0.5 * ti.sin(p.y + 5 * t))

    d = ti.max(d, d1, d2)
    
    op = ti.min(d, op)

    # plane at (0,0,0) pointing in y direction
    plane_distance = (
        p.y
    )  # - ti.exp(-0.01*p.xz.norm())*(0.4*ti.sin(p.x-t) + 0.2*ti.cos(p.z))
    return ti.min(op, plane_distance)


@ti.func
def ray_march(ray_origin, ray_direction):
    """
    Finds the distance `d` in
        closest_surface_coordinate = `ray_origin` + d * `ray_direction`
    via ray marching. See wiki or so.
    """
    d = 0.0
    d_min = MAX_DIST
    for i in range(MAX_STEPS):
        p = ray_origin + d * ray_direction
        distance = get_distance(p)
        d_min = ti.min(d_min, distance)
        d += distance

        if d > MAX_DIST or abs(distance) < EPS:
            break

    return d, d_min


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

    if d[0] < d_ref:
        diffuse_light *= 0.5

    return ti.max(diffuse_light, 0.0)


@ti.kernel
def paint():
    for i, j in pixels:  # Parallelized over all pixels

        x = (i - 0.5 * RES_X) / RES_Y  # / RES_Y is not a typo. Makes dx == dy
        y = (j - 0.5 * RES_Y) / RES_Y

        ray_origin = cam_pos[None]
        ray_direction = ti.Vector((x, y, 1.0)).normalized()

        ray_direction.yz = ti.Matrix.rotation2d(-cam_rot[None].y) @ ray_direction.yz
        ray_direction.xz = ti.Matrix.rotation2d(-cam_rot[None].x) @ ray_direction.xz

        # find the intersection point `p` in in distance `d` along the ray
        d, d_min = ray_march(ray_origin, ray_direction)
        p = ray_origin + ray_direction * d

        t = time[None].x

        light_pos = ti.Vector((0.0, 10.0, 6.0))
        light_pos.x += ti.sin(0.2*t) * 30
        light_pos.z += ti.cos(0.2*t) * 30

        light_pos2 = ti.Vector((0.0, 10.0, 6.0))
        light_pos2.x -= ti.sin(-0.4*t) * 30
        light_pos2.z -= ti.cos(0.4*t) * 30

        diffuse_light = get_light(p, t, light_pos)
        diffuse_light2 = get_light(p, t, light_pos2)

        sunlight = ti.Vector((0.9922, 0.9842, 0.8275))
        blue = ti.Vector((0., 0.7, 1.0))
        
        pixels[i, j] = sunlight * diffuse_light
        pixels[i, j] += blue * diffuse_light2
        # pixels[i, j] += yellowish * diffuse_light2
        # pixels[i, j] = diffuse_light, 0.0, diffuse_light2


gui = ti.GUI("Move with WASD + mouse. Space/Shift = Up/Down ", res=(RES_X, RES_Y))

dt = 0.02
CAM_SPEED = 0.1
CAM_ROT_SPEED = 0.1

mx0, my0 = 0.0, 0.0

while not gui.get_event(ti.GUI.ESCAPE):

    if gui.is_pressed(ti.GUI.LMB):
        mx, my = gui.get_cursor_pos()
        cam_rot[None].x += (mx - mx0) * 3.1415
        cam_rot[None].y += (my - my0) * 6.2831
        mx0, my0 = gui.get_cursor_pos()
    else:
        mx0, my0 = gui.get_cursor_pos()

    front = ti.Matrix.rotation2d(cam_rot[None].x)

    if gui.is_pressed("w"):
        cam_pos[None].zx += CAM_SPEED * front[:, 0]
    if gui.is_pressed("s"):
        cam_pos[None].zx -= CAM_SPEED * front[:, 0]
    if gui.is_pressed("a"):
        cam_pos[None].zx -= CAM_SPEED * front[:, 1]
    if gui.is_pressed("d"):
        cam_pos[None].zx += CAM_SPEED * front[:, 1]
    if gui.is_pressed(ti.GUI.SPACE):
        cam_pos[None].y += CAM_SPEED
    if gui.is_pressed(ti.GUI.SHIFT):
        cam_pos[None].y -= CAM_SPEED

    paint()
    gui.set_image(pixels)
    gui.show()
    time[None] += dt
