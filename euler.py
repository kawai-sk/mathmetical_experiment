import matplotlib.pyplot as plt
from math import sin,cos

g = 9.8
l_1 = 0.5
l_2 = 0.5
m_1 = 0.5
m_2 = 0.3
M = m_2/(m_1+m_2)
l = l_2/l_1
w = (g/l_1)**0.5

def t1(t,X):#X = [theta1,theta2,thetadot1,thetadot2]
    return X[2]
def t2(t,X):
    return X[3]
def td1(t,X):
    return (w**2*l*(-sin(X[0])+M*cos(X[0]-X[1])*sin(X[1]))-M*l*(X[2]**2*cos(X[0]-X[1])+l*X[3]**2)*sin(X[0]-X[1]))/(l-M*l*cos(X[0]-X[1])**2)
def td2(t,X):
    return (w**2*cos(X[0]-X[1])*sin(X[0])-w**2*sin(X[1])+(X[2]**2+M*l*X[3]**2*cos(X[0]-X[1]))*sin(X[0]-X[1]))/(l-M*l*cos(X[0]-X[1])**2)

def energy(X):#エネルギー(ハミルトニアン)
    t = 0.5*m_1*l_1**2*X[2]**2+0.5*m_2*(l_1**2*X[2]**2+l_2**2*X[3]**2+2*l_1*l_2*X[2]*X[3]*cos(X[0]-X[1]))
    u = -m_1*l_1*g*cos(X[0])-m_2*g*(l_1*cos(X[0])+l_2*cos(X[1]))
    return t+u

def euler(n,a,h,Y_a):#陽的Euler法
    Y = Y_a
    F = [t1,t2,td1,td2]
    T1 = [Y[0]]
    T2 = [Y[1]]
    Td1 = [Y[2]]
    Td2 = [Y[3]]
    E = [energy(Y_a)]
    for i in range(0,n):
        Y = [Y[j]+h*F[j](a+h*i,Y) for j in range(4)]
        T1.append(cos(Y[0]))
        T2.append(cos(Y[1]))
        Td1.append(Y[2])
        Td2.append(Y[3])
        E.append(energy(Y))
    x = [a+h*i for i in range(n+1)]
    plt.plot(x,T1,"-")
    plt.xlabel("time")
    plt.ylabel("theta1")
    plt.show()
    plt.plot(x,T2,"-")
    plt.xlabel("time")
    plt.ylabel("theta2")
    plt.show()
    plt.plot(x,Td1,"-")
    plt.xlabel("time")
    plt.ylabel("thetadot1")
    plt.show()
    plt.plot(x,Td2,"-")
    plt.xlabel("time")
    plt.ylabel("thetadot2")
    plt.show()
    plt.plot(x,E,"-")
    plt.xlabel("time")
    plt.ylabel("Energy")
    plt.show()
n = 3000
a = 0
h = 0.01
Y = [0.1,0,0,0]
euler(n,a,h,Y)
n = 30000
h = 0.001
euler(n,a,h,Y)
