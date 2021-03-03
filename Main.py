import matplotlib.pyplot as plt
import math
from matrix import *


re = 800  # big link
rf = 300  # small link
base_radius = 370  # radius of top platform
end_platform_radius = 80  # radius of end-effector

end_position = [-200, 400, -805]


def GetActiveJointAngle(x0, y0, z0):
    y1 = -base_radius
    y0 -= end_platform_radius

    a = (x0 * x0 + y0 * y0 + z0 * z0 + rf * rf - re * re - y1 * y1) / (2 * z0)
    b = (y1 - y0) / z0
    d = -(a + b * y1) * (a + b * y1) + rf * (b * b * rf + rf)

    if d < 0:
        print("The point doesnt exist")

    yj = (y1 - a * b - np.sqrt(d)) / (b * b + 1)
    zj = a + b * yj

    if yj > y1:
        theta = 180.0 * np.arctan(-zj / (y1 - yj)) / np.pi + 180
    else:
        theta = 180.0 * np.arctan(-zj / (y1 - yj)) / np.pi

    return theta


def IK(x0, y0, z0):
    cos120 = -0.5
    sin120 = np.sin(np.deg2rad(120))

    theta1 = GetActiveJointAngle(x0, y0, z0)
    theta2 = GetActiveJointAngle(x0 * cos120 + y0 * sin120, y0 * cos120 - x0 * sin120, z0)
    theta3 = GetActiveJointAngle(x0 * cos120 - y0 * sin120, y0 * cos120 + x0 * sin120, z0)

    return [theta1, theta2, theta3]


def FK(theta1, theta2, theta3):
    tan30 = np.tan(np.deg2rad(30))
    tan60 = np.tan(np.deg2rad(60))
    sin30 = np.sin(np.deg2rad(30))

    t = (base_radius - end_platform_radius) * tan30 / 2
    dtr = np.pi / 180

    theta1 *= dtr
    theta2 *= dtr
    theta3 *= dtr

    y1 = -(t + rf * np.cos(theta1))
    z1 = -rf * np.sin(theta1)

    y2 = (t + rf * np.cos(theta2)) * sin30
    x2 = y2 * tan60
    z2 = -rf * np.sin(theta2)

    y3 = (t + rf * np.cos(theta3)) * sin30
    x3 = -y3 * tan60
    z3 = -rf * np.sin(theta3)

    dnm = (y2 - y1) * x3 - (y3 - y1) * x2

    w1 = y1 * y1 + z1 * z1
    w2 = x2 * x2 + y2 * y2 + z2 * z2
    w3 = x3 * x3 + y3 * y3 + z3 * z3

    a1 = (z2 - z1) * (y3 - y1) - (z3 - z1) * (y2 - y1)
    b1 = -((w2 - w1) * (y3 - y1) - (w3 - w1) * (y2 - y1)) / 2

    a2 = -(z2 - z1) * x3 + (z3 - z1) * x2
    b2 = ((w2 - w1) * x3 - (w3 - w1) * x2) / 2

    a = a1 * a1 + a2 * a2 + dnm * dnm
    b = 2 * (a1 * b1 + a2 * (b2 - y1 * dnm) - z1 * dnm * dnm)
    c = (b2 - y1 * dnm) * (b2 - y1 * dnm) + b1 * b1 + dnm * dnm * (z1 * z1 - re * re)

    d = b * b - 4 * a * c

    if (d < 0):
        print("error")

    z0 = -0.5 * (b + np.sqrt(d)) / a
    x0 = (a1 * z0 + b1) / dnm
    y0 = (a2 * z0 + b2) / dnm
    return [x0, y0, z0]


def DrawCircle(ax, pos, radius):
    theta = np.linspace(0, 2 * np.pi, 200)

    x = radius * np.sin(theta) + pos[0]
    y = radius * np.cos(theta) + pos[1]

    ax.plot(x, y, pos[2])


def GetLegLinksPos(ee_pos, q0, q1):
    link_start_pos = np.linalg.multi_dot([Rz(q0), Ty(-base_radius)])[0:3, 3]
    elbow_pos = np.linalg.multi_dot([Rz(q0), Ty(-base_radius), Rx(q1), Ty(-rf)])[0:3, 3]
    wrist_pos = np.linalg.multi_dot([Tx(ee_pos[0]), Ty(ee_pos[1]), Tz(ee_pos[2]), Rz(q0), Ty(-end_platform_radius)])[0:3, 3]

    return link_start_pos, elbow_pos, wrist_pos


def DrawLeg(ax, base_center, leg_start_pos, elbow_pos, wrist_pos, end_pos):
    x = [base_center[0], leg_start_pos[0], elbow_pos[0], wrist_pos[0], end_pos[0]]
    y = [base_center[1], leg_start_pos[1], elbow_pos[1], wrist_pos[1], end_pos[1]]
    z = [base_center[2], leg_start_pos[2], elbow_pos[2], wrist_pos[2], end_pos[2]]

    ax.plot3D(x, y, z)


def Visualisation(ee_pos, q):
    ax = plt.axes(projection='3d')
    base_center = [0, 0, 0]

    # Draw base
    DrawCircle(ax, base_center, base_radius)

    # Draw end-effector
    DrawCircle(ax, ee_pos, end_platform_radius)

    # Sequence 1
    leg_start_pos, elbow_pos, wrist_pos = GetLegLinksPos(ee_pos, 0, np.deg2rad(q[0]))
    DrawLeg(ax, base_center, leg_start_pos, elbow_pos, wrist_pos, ee_pos)

    # Sequence 2
    leg_start_pos, elbow_pos, wrist_pos = GetLegLinksPos(ee_pos, np.deg2rad(120), np.deg2rad(q[1]))
    DrawLeg(ax, base_center, leg_start_pos, elbow_pos, wrist_pos, ee_pos)

    # Sequence 3
    leg_start_pos, elbow_pos, wrist_pos = GetLegLinksPos(ee_pos, np.deg2rad(-120), np.deg2rad(q[2]))
    DrawLeg(ax, base_center, leg_start_pos, elbow_pos, wrist_pos, ee_pos)

    print("RE length {}.".format(math.dist(elbow_pos, wrist_pos)))

    ax.set_xlim(-500, 500)
    ax.set_ylim(-500, 500)
    ax.set_zlim(-1000, 0)
    plt.show()


thetas = IK(end_position[0], end_position[1], end_position[2])
print(f'IK angles: {thetas}')
Visualisation(end_position, thetas)
