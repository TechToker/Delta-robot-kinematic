import numpy as np
import matplotlib.pyplot as plt
import math


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


rf = 300  # small link
re = 800  # big link

base_radius = 250  # radius of top platform
end_platform_radius = 100  # radius of end-effector

position = [0, -400, -600]


def FixedIK(x0, y0, z0):
    delta_radius = base_radius - end_platform_radius + y0

    x = np.sqrt(z0**2 + delta_radius**2)
    alpha = np.arctan2(abs(z0), delta_radius)
    beta = np.pi - alpha
    eta = np.arccos((rf**2 + x**2 - re**2) / (2 * rf * x))

    print(f'Eta: {eta}')
    print(f'Eta: {np.rad2deg(eta)}')

    q1 = eta - beta

    print(f'q1:{np.rad2deg(q1)}')

    q2 = np.pi - np.arccos((rf**2 + re**2 - x**2) / (2 * rf * re))

    #print(f'q2: {np.rad2deg(np.pi - np.arccos((rf**2 + re**2 - x**2) / (2 * rf * re)))}')

    return [q1, q2]


def calc_angle(x0, y0, z0):
    b_r = base_radius * 2 * np.sqrt(3)
    e_r = end_platform_radius * 2 * np.sqrt(3)

    y1 = -0.5 * np.tan(np.pi/6) * b_r
    y0 -= 0.5 * np.tan(np.pi/6) * e_r

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


    print(f'New Pos2: {[0, y1, zj]}')

    return theta #, [0, -base_radius + yj, zj]


def DrawPseudoCircle(ax, pos, radius):
    theta = np.linspace(0, 2 * np.pi, 200)

    y = radius * np.sin(theta) + pos[1]
    z = radius * np.cos(theta) + pos[2]

    zeros = np.zeros(200)

    ax.scatter(pos[0], pos[1], pos[2] + 30)
    ax.plot(zeros, y, z)


ax = plt.axes(projection='3d')

DrawPseudoCircle(ax, [0, - base_radius, 0], rf)
DrawPseudoCircle(ax, [position[0], position[1] - end_platform_radius, position[2]], re)

#anw = FixedIK(position[0], position[1], position[2])
#thetas = [anw[0], 0, 0]


def RobotFK(ax, q0, q1, q2):

    T = np.linalg.multi_dot([Rz(q0),
                             Ty(-base_radius)])

    pos1 = T[0:3, 3]

    T = np.linalg.multi_dot([Rz(q0),
                             Ty(-base_radius),
                             Rx(-q1),
                             Ty(-rf)])

    pos2 = T[0:3, 3]
    print(f'Correct Pos2: {pos2}')
    ax.scatter(pos2[0], pos2[1], pos2[2])

    T = np.linalg.multi_dot([Rz(q0),
                             Ty(-base_radius),
                             Rx(-q1),
                             Ty(-rf),
                             Rx(q2),
                             Ty(-re)])

    pos3 = T[0:3, 3]

    x = [pos1[0], pos2[0], pos3[0]]
    y = [pos1[1], pos2[1], pos3[1]]
    z = [pos1[2], pos2[2], pos3[2]]

    ax.plot3D(x, y, z)

    ax.set_xlim(-700, 700)
    ax.set_ylim(-700, 700)
    ax.set_zlim(-1400, 0)


def DrawCircle(ax, pos, radius):
    theta = np.linspace(0, 2 * np.pi, 200)

    x = radius * np.sin(theta) + pos[0]
    y = radius * np.cos(theta) + pos[1]

    ax.plot(x, y, pos[2])

# Draw base
DrawCircle(ax, [0, 0, 0], base_radius)

# Draw end-effector
DrawCircle(ax, position, end_platform_radius)

#RobotFK(ax, 0, anw[0], anw[1])
print("")
print("New IK")



def IK(x0, y0, z0):
    cos120 = -0.5
    sin120 = np.sin(np.deg2rad(120))

    theta1 = calc_angle(x0, y0, z0)
    theta2 = calc_angle(x0 * cos120 + y0 * sin120, y0 * cos120 - x0 * sin120, z0)
    theta3 = calc_angle(x0 * cos120 - y0 * sin120, y0 * cos120 + x0 * sin120, z0)

    return [theta1, theta2, theta3]


thetas = IK(position[0], position[1], position[2])
print(f'New q1: {thetas[0]}')


def VisualisationFK(ax, q0, q1):

    T = np.linalg.multi_dot([Rz(q0),
                             Ty(-base_radius)])

    pos1 = T[0:3, 3]

    T = np.linalg.multi_dot([Rz(q0),
                             Ty(-base_radius),
                             Rx(q1),
                             Ty(-rf)])

    pos2 = T[0:3, 3]

    # ax.scatter(pos1[0], pos1[1], pos1[2])
    # ax.scatter(pos2[0], pos2[1], pos2[2])

    return pos1, pos2



def VisualisationFK_ee(ax, ee_pos, q0):

    T = np.linalg.multi_dot([Tx(ee_pos[0]),
                             Ty(ee_pos[1]),
                             Tz(ee_pos[2])])

    pos0 = T[0:3, 3]

    T = np.linalg.multi_dot([Tx(ee_pos[0]),
                             Ty(ee_pos[1]),
                             Tz(ee_pos[2]),
                             Rz(q0),
                             Ty(-end_platform_radius)])

    pos1 = T[0:3, 3]

    # ax.scatter(pos1[0], pos1[1], pos1[2])
    # ax.scatter(pos0[0], pos0[1], pos0[2])

    return pos0, pos1


def PlotLeg(ax, center, pos0, pos1, posEnd, posEe):

    # Sequence 1
    x = [center[0], pos0[0], pos1[0], posEe[0], posEnd[0]]
    y = [center[1], pos0[1], pos1[1], posEe[1], posEnd[1]]
    z = [center[2], pos0[2], pos1[2], posEe[2], posEnd[2]]

    ax.plot3D(x, y, z)


def CalcDist(pos1, pos2):
    return np.sqrt((pos2[0] - pos1[0])**2 + (pos2[1] - pos1[1])**2 + (pos2[2] - pos1[2])**2)


def Visualisation(ax, ee_pos, q):
    #ax = plt.axes(projection='3d')
    center = [0, 0, 0]

    # Draw base
    DrawCircle(ax, center, base_radius)

    # Draw end-effector
    DrawCircle(ax, ee_pos, end_platform_radius)

    # Sequence 1
    # print(ax)
    # print(np.deg2rad(q[0]))

    pos0, posElbow = VisualisationFK(ax, 0, np.deg2rad(q[0]))
    posEndEffector, posWrist = VisualisationFK_ee(ax, ee_pos, np.deg2rad(0))

    PlotLeg(ax, center, pos0, posElbow, posEndEffector, posWrist)

    ax.scatter(posElbow[0], posElbow[1], posElbow[2])
    ax.scatter(posWrist[0], posWrist[1], posWrist[2])

    print(f'Pos el{posElbow}')
    print(f'Pos wr{posWrist}')

    dst = math.dist(posWrist, posElbow)
    print(f'Distance: {dst}')

    # # Sequence 2
    # pos0, pos1 = VisualisationFK(ax, np.deg2rad(120), np.deg2rad(q[1]))
    # posEnd, e_pos = VisualisationFK_ee(ax, ee_pos, np.deg2rad(120))
    #
    # PlotLeg(ax, center, pos0, pos1, posEnd, e_pos)
    # #
    # # dst1 = math.dist(posEnd, pos1)
    # # print(dst1)
    # #
    # # Sequence 3
    # pos0, pos1 = VisualisationFK(ax, np.deg2rad(-120), np.deg2rad(q[2]))
    # posEnd, e_pos = VisualisationFK_ee(ax, ee_pos, np.deg2rad(-120))
    #
    # PlotLeg(ax, center, pos0, pos1, posEnd, e_pos)
    #
    # dst2 = math.dist(posEnd, pos1)
    # print(dst2)

    ax.set_xlim(-500, 500)
    ax.set_ylim(-500, 500)
    ax.set_zlim(-1000, 0)
    plt.show()

Visualisation(ax, position, thetas)
