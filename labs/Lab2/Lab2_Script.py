from cmath import sin
import modern_robotics as mr
import numpy as np

def PoE_2D(theta1, theta2):

    M = np.array([
        [1, 0, 0, (105 + 55.95 + 57.75)/1000],
        [0, 1, 0, (-98 + (-15))/1000],
        [0, 0, 1, (100 - 12.31)/1000],
        [0, 0, 0, 1]
    ])
    omega_1 = np.array([0,1,0])
    omega_2 = np.array([0,1,0])
    q_1 = np.array([0,-98,100])/1000
    q_2 = np.array([105,-113,100])/1000
    v_1 = - np.cross(omega_1,q_1) 
    v_2 = - np.cross(omega_2,q_2) 
    
    S1 = mr.MatrixExp6(mr.VecTose3(np.concatenate((omega_1,v_1),axis=0))*theta1)
    S2 = mr.MatrixExp6(mr.VecTose3(np.concatenate((omega_2,v_2),axis=0))*theta2)
    T_02 = S1.dot(S2).dot(M)
    return T_02

def DH(d, theta, a, alpha):
    return np.array([
        [np.cos(theta), -np.sin(theta)*np.cos(alpha), np.sin(theta)*np.sin(alpha), a*np.cos(theta)],
        [np.sin(theta), np.cos(theta) * np.cos(alpha), -np.cos(theta) * np.sin(alpha), a*np.sin(theta)],
        [0, np.sin(alpha), np.cos(alpha), d],
        [0, 0, 0, 1]
    ])

def DH_2D(theta1, theta2):
    Ts0 = np.array([
        [1,0,0,0],
        [0,0,1,-98/1000.0],
        [0,-1,0,100/1000.0],
        [0,0,0,1]
    ])

    T2b = np.array([
        [1,0,0,0],
        [0,0,1,-12.31/1000.0],
        [0,-1,0,0],
        [0,0,0,1]
    ])

    T01 = DH(-15/1000.0, theta1, 105/1000.0, np.pi)
    T12 = DH(0, theta2, (55.95+57.75)/1000.0, 0)

    return Ts0.dot(T01).dot(T12).dot(T2b)

# poe = PoE_2D(np.pi/3, np.pi/2)
# dh = DH_2D(np.pi/3, np.pi/2)

poe = PoE_2D(0/180*np.pi, 2/180*np.pi)
dh = DH_2D(90/180*np.pi, 90/180*np.pi)
print(poe)
# print(dh)

def PoE_3D(theta1, theta2, theta3):
    M = np.array([
        [1, 0, 0, (105 + 55.95 + 57.75)/1000],
        [0, 1, 0, (-98 + (-15))/1000],
        [0, 0, 1, (100 - 12.31)/1000],
        [0, 0, 0, 1]
    ])
    omega_1 = np.array([0,1,0])
    omega_2 = np.array([0,0,1])
    omega_3 = np.array([0,-1,0])

    q_1 = np.array([0,-98,100])/1000
    q_2 = np.array([0,-98,100])/1000
    q_3 = np.array([105,-113,100])/1000

    v_1 = - np.cross(omega_1,q_1)
    v_2 = - np.cross(omega_2,q_2) 
    v_3 = - np.cross(omega_3,q_3) 
    
    S1 = mr.MatrixExp6(mr.VecTose3(np.concatenate((omega_1,v_1),axis=0))*theta1)
    S2 = mr.MatrixExp6(mr.VecTose3(np.concatenate((omega_2,v_2),axis=0))*theta2)
    S3 = mr.MatrixExp6(mr.VecTose3(np.concatenate((omega_3,v_3),axis=0))*theta3)

    T_05 = S1.dot(S2).dot(S3).dot(M)
    return T_05

