import sympy
from sympy import *
import matplotlib.pyplot as plt
x = Symbol("x",real=True)

l_1 = l_2 = 0.5
m_1 = 0.5
m_2 = 0.3
g = 9.8

def newton(f,a_0,x):#Newton法
    g = diff(f,x)
    x_0 = a_0
    while abs(f.subs(x,x_0)) >= 10**(-9):
        x_0 -= f.subs(x,x_0)/g.subs(x,x_0)
    return x_0

def energy(X):#エネルギー(ハミルトニアン)
    t = 0.5*m_1*l_1**2*X[2]**2+0.5*m_2*(l_1**2*X[2]**2+l_2**2*X[3]**2+2*l_1*l_2*X[2]*X[3]*cos(X[0]-X[1]))
    u = -m_1*l_1*g*cos(X[0])-m_2*g*(l_1*cos(X[0])+l_2*cos(X[1]))
    return t+u

def symplectic_euler(n,a,h,T):#T=[t1,t2,td1,td2]
    q1 = T[0]
    q2 = T[1]
    p1 = (m_1+m_2)*l_1**2*T[2]+m_1*l_1*l_2*cos(T[0]-T[1])*T[3]
    p2 = m_2*l_1*l_2*cos(T[0]-T[1])*T[2]+m_2*l_2**2*T[3]
    P1 = []
    P2 = []
    E = []
    for i in range(n):
        F_1 = l_2*p1-l_1*p2*cos(x) #Newton法の準備
        F_2 = (m_1*l_1+m_2*l_1)*p2 - m_2*l_2*p1*cos(x)
        G = m_1+m_2-m_2*cos(x)**2
        f = x - (q1-q2+h*F_1/(l_1**2*l_2*G)-h*F_2/(m_2*l_1*l_2**2*G))

        q = newton(f,q1-q2,x) #Newton法を適用
        F1 = F_1.subs(x,q)
        F2 = F_2.subs(x,q)
        GG = G.subs(x,q)

        q1 = q1 + h*F1/(l_1**2*l_2*GG) #q,pを求める
        q2 = q1 - q
        p1 = p1 + h*(-F1*F2*sin(q)/((l_1*l_2*GG)**2)-(m_1+m_2)*l_1*g*sin(q1))
        p2 = p2 + h*(F1*F2*sin(q)/((l_1*l_2*GG)**2)-m_2*l_2*g*sin(q2))

        dq = p1 - p2 #dtを求める
        det = (m_1+m_2)*m_2*(l_1*l_2)**2-(m_2*l_1*l_2*cos(dq))**2
        td1 = (m_2*l_2**2*p1-m_2*l_1*l_2*cos(dq)*p2)/det
        td2 = (-m_2*l_1*l_2*cos(dq)*p1+(m_1+m_2)*l_1**2*p2)/det
        E.append(energy([q1,q2,td1,td2]))

        if -0.01<td1<0.01:#(d/dt)theta_1=0の場合
            P1.append(q1)
            P2.append(q2)
    #plt.scatter(P1,P2,marker=".",s=5)#(d/dt)theta_1=0の点の描画
    #plt.xlabel("theta_1")
    #plt.ylabel("theta_2")
    #plt.show()

    x1 = [a+h*i for i in range(1,n+1)] #エネルギー(ハミルトニアン)の描画
    plt.plot(x1,E,"-")
    plt.xlabel("time")
    plt.ylabel("Energy")
    plt.show()

n = 10000
a = 0
h = 0.01
T = [0.1,0,0,0]
symplectic_euler(n,a,h,T)
