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
        theta = np.arctan(-zj / (y1 - yj)) + np.pi
    else:
        theta = np.arctan(-zj / (y1 - yj))

    return theta


def IK(x0, y0, z0):
    r1 = rf
    r2 = re

    cos120 = -0.5
    sin120 = np.sin(np.deg2rad(120))

    theta1 = GetActiveJointAngle(x0, y0, z0)
    theta2 = GetActiveJointAngle(x0 * cos120 + y0 * sin120, y0 * cos120 - x0 * sin120, z0)
    theta3 = GetActiveJointAngle(x0 * cos120 - y0 * sin120, y0 * cos120 + x0 * sin120, z0)

    alpha1 = np.arcsin(y0 / re)
    alpha2 = np.arcsin(y0 * cos120 - x0 * sin120 / re)
    alpha3 = np.arcsin(y0 * cos120 + x0 * sin120 / re)

    q1 = [theta1, theta2, theta3]
    q3 = [alpha1, alpha2, alpha3]

    beta1 = np.arccos((-(r1 * np.cos(theta1) - x0) * np.cos(theta1) + np.sqrt(
        -r1 ** 2 * np.cos(theta1) ** 2 + 2 * r1 * x0 * np.cos(theta1) + r2 ** 2 * np.cos(alpha1) ** 2 - x0 ** 2) * np.sin(theta1)) / (r2 * np.cos(alpha1)))

    beta2 = np.arccos((-(r1 * np.cos(theta2) - (x0 * cos120 + y0 * sin120)) * np.cos(theta2) + np.sqrt(
        -r1 ** 2 * np.cos(theta2) ** 2 + 2 * r1 * (x0 * cos120 + y0 * sin120) * np.cos(theta2) + r2 ** 2 * np.cos(alpha2) ** 2 - (x0 * cos120 + y0 * sin120) ** 2) * np.sin(theta2)) / (r2 * np.cos(alpha2)))

    beta3 = np.arccos((-(r1 * np.cos(theta3) - (x0 * cos120 - y0 * sin120)) * np.cos(theta3) + np.sqrt(
        -r1 ** 2 * np.cos(theta3) ** 2 + 2 * r1 * (x0 * cos120 - y0 * sin120) * np.cos(theta3) + r2 ** 2 * np.cos(alpha3) ** 2 - (x0 * cos120 - y0 * sin120) ** 2) * np.sin(theta3)) / (r2 * np.cos(alpha3)))

    q2 = [beta1, beta2, beta3]

    print(q1)
    print(q2)
    print(q3)

    return q1, q2, q3


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
#Visualisation(end_position, thetas)


####
space = [[-800, 800], [-800, 800], [-1100, -100]]
points_per_axis = 20

def JacobianTheta(q1, q2, q3):

    r1 = rf

    eq1 = -r1 * (-np.sin(q1[0]) * np.cos(q1[0] + q2[0]) - np.cos(q1[0]) * np.sin(q1[0]+q2[0]))
    eq2 = -r1 * (-np.sin(q1[1]) * np.cos(q1[1] + q2[1]) - np.cos(q1[1]) * np.sin(q1[1] + q2[1]))
    eq3 = -r1 * (-np.sin(q1[2]) * np.cos(q1[2] + q2[2]) - np.cos(q1[2]) * np.sin(q1[2] + q2[2]))

    J = np.diag([eq1, eq2, eq3])

    return J

def JacobianZ(q1, q2, q3):

    eq1 = np.hstack([-np.cos(q1[0] + q2[0]), -np.tan(q3[0]), np.sin(q1[0] + q2[0])])
    eq2 = np.hstack([-np.cos(q1[1] + q2[1]), -np.tan(q3[1]), np.sin(q1[1] + q2[1])])
    eq3 = np.hstack([-np.cos(q1[2] + q2[2]), -np.tan(q3[2]), np.sin(q1[2] + q2[2])])

    return np.vstack([eq1, eq2, eq3])


def CalculateDeflections():
    x_pos = []
    y_pos = []
    z_pos = []
    deflections = []

    x_linSpace = np.linspace(space[0][0], space[0][1], points_per_axis)
    y_linSpace = np.linspace(space[1][0], space[1][1], points_per_axis)
    z_linSpace = np.linspace(space[2][0], space[2][1], points_per_axis)

    for x in x_linSpace:
        for y in y_linSpace:
            for z in z_linSpace:

                end_position = [x, y, z]
                q1, q2, q3 = IK(end_position[0], end_position[1], end_position[2])
                J_theta = JacobianTheta(q1, q2, q3)
                J_z = JacobianZ(q1, q2, q3)
                J = np.dot(np.linalg.inv(J_theta), J_z)

                m = np.sqrt(np.linalg.det(np.dot(J, np.transpose(J))))

                print(m)

                x_pos.append(x)
                y_pos.append(y)
                z_pos.append(z)
                deflections.append(m)

    return np.array(x_pos), np.array(y_pos), np.array(z_pos), np.array(deflections)


def plotDeflectionMap(x_pos, y_pos, z_pos, deflection, colormap, s):
    plt.figure()
    ax = plt.axes(projection='3d')

    ax.set_xlim3d(space[0][0], space[0][1])
    ax.set_ylim3d(space[1][0], space[1][1])
    ax.set_zlim3d(space[2][0], space[2][1])

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.colorbar(ax.scatter3D(x_pos, y_pos, z_pos, c=deflection, cmap=colormap, s=s))
    plt.show()


color_map = plt.cm.get_cmap('viridis', 12)

xScatter, yScatter, zScatter, dScatter = CalculateDeflections()
plotDeflectionMap(xScatter, yScatter, zScatter, dScatter, color_map, 60)

r1 = 300
r2 =800
theta1 = 0.5
alpha1 = 0.5
x0 = 200
beta5 = np.arccos((-(r1 * np.cos(theta1) - x0) * np.cos(theta1) + np.sqrt(
    -r1 ** 2 * np.cos(theta1) ** 2 + 2 * r1 * x0 * np.cos(theta1) + r2 ** 2 * np.cos(alpha1) ** 2 - x0 ** 2) * np.sin(
    theta1)) / (r2 * np.cos(alpha1)))

print(beta5)
