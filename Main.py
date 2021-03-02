import numpy as np
import matplotlib.pyplot as plt

def Rx(q):
    T = np.array([[1, 0, 0, 0],
                  [0, np.cos(q), -np.sin(q), 0],
                  [0, np.sin(q), np.cos(q), 0],
                  [0, 0, 0, 1]])
    return T


def Ry(q):
    T = np.array([[np.cos(q), 0, np.sin(q), 0],
                  [0, 1, 0, 0],
                  [-np.sin(q), 0, np.cos(q), 0],
                  [0, 0, 0, 1]])
    return T


def Rz(q):
    T = np.array([[np.cos(q), -np.sin(q), 0, 0],
                  [np.sin(q), np.cos(q), 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]])
    return T


def Tx(x):
    T = np.array([[1, 0, 0, x],
                  [0, 1, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]])
    return T


def Ty(y):
    T = np.array([[1, 0, 0, 0],
                  [0, 1, 0, y],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]])
    return T


def Tz(z):
    T = np.array([[1, 0, 0, 0],
                  [0, 1, 0, 0],
                  [0, 0, 1, z],
                  [0, 0, 0, 1]])
    return T


def dRx(q):
    T = np.array([[0, 0, 0, 0],
                  [0, -np.sin(q), -np.cos(q), 0],
                  [0, np.cos(q), -np.sin(q), 0],
                  [0, 0, 0, 0]])
    return T


def dRy(q):
    T = np.array([[-np.sin(q), 0, np.cos(q), 0],
                  [0, 0, 0, 0],
                  [-np.cos(q), 0, -np.sin(q), 0],
                  [0, 0, 0, 0]])
    return T


def dRz(q):
    T = np.array([[-np.sin(q), -np.cos(q), 0, 0],
                  [np.cos(q), -np.sin(q), 0, 0],
                  [0, 0, 0, 0],
                  [0, 0, 0, 0]])
    return T


def dTx(x):
    T = np.array([[0, 0, 0, 1],
                  [0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [0, 0, 0, 0]])
    return T


def dTy(y):
    T = np.array([[0, 0, 0, 0],
                  [0, 0, 0, 1],
                  [0, 0, 0, 0],
                  [0, 0, 0, 0]])
    return T


def dTz(z):
    T = np.array([[0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [0, 0, 0, 1],
                  [0, 0, 0, 0]])
    return T


re = 800  # big link
rf = 300  # small link
base_radius = 370  # radius of top platform
end_platform_radius = 80  # radius of end-effector


def calc_angle(x0, y0, z0):
    y1 = -0.5 * 0.57735 * base_radius
    y0 -= 0.5 * 0.57735 * end_platform_radius
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

    return theta, [0, -base_radius + yj, zj]


def IK(x0, y0, z0):
    cos120 = -0.5
    sin120 = np.sin(np.deg2rad(120))

    theta1, elbow_pos_1 = calc_angle(x0, y0, z0)
    theta2, elbow_pos_2 = calc_angle(x0 * cos120 + y0 * sin120, y0 * cos120 - x0 * sin120, z0)
    theta3, elbow_pos_3 = calc_angle(x0 * cos120 - y0 * sin120, y0 * cos120 + x0 * sin120, z0)

    return [theta1, theta2, theta3], [elbow_pos_1, elbow_pos_2, elbow_pos_3]


#print(IK(100, 300, -500))


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


#print(FK(15, -65, -40))


position = [0, 0, -900]
thetas, elbows_pos = IK(position[0], position[1], position[2])

#print(thetas)


def DrawCircle(ax, pos, radius):
    theta = np.linspace(0, 2 * np.pi, 200)

    x = radius * np.sin(theta) + pos[0]
    y = radius * np.cos(theta) + pos[1]

    ax.plot(x, y, pos[2])


def Visualisation(ee_pos, elbows_positions):
    ax = plt.axes(projection='3d')
    center = [0, 0, 0]

    # Draw base
    DrawCircle(ax, center, base_radius)

    # Draw end-effector
    DrawCircle(ax, ee_pos, end_platform_radius)

    elbow_point = elbows_positions[0]

    print(elbows_positions)

    # Sequence 1
    x = [center[0], center[0], elbow_point[0], ee_pos[0], ee_pos[0]]
    y = [center[1], center[1] - base_radius, elbow_point[1], ee_pos[1] - end_platform_radius, ee_pos[1]]
    z = [center[2], center[2], elbow_point[2], ee_pos[2], ee_pos[2]]

    ax.plot3D(x, y, z)

    # Sequence 2
    elbow_point = elbows_positions[1]

    cos120 = -0.5
    sin120 = np.sin(np.deg2rad(120))

    #theta2, elbow_pos_2 = calc_angle(x0 * cos120 + y0 * sin120, y0 * cos120 - x0 * sin120, z0)

    x = np.array([center[0], center[0], elbow_point[0], ee_pos[0], ee_pos[0]])
    y = np.array([center[1], center[1] - base_radius, elbow_point[1], ee_pos[1] - end_platform_radius, ee_pos[1]])
    z = [center[2], center[2], elbow_point[2], ee_pos[2], ee_pos[2]]

    x = x * cos120 + y * sin120
    y = y * cos120 - x * sin120

    ax.plot3D(x, y, z)

    ax.set_xlim(-500, 500)
    ax.set_ylim(-500, 500)
    ax.set_zlim(-1000, 0)
    plt.show()


Visualisation(position, elbows_pos)
