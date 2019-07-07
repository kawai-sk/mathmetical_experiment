#%matplotlib notebook #jupyterで描画する場合
#%matplotlib inline
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math

C_1 = 1/9
C_2 = 1
L = 1/7
m_0 = -0.5
m_1 = -0.8
B_p = 1

def coalitioned_rungekutta(n,a,h,Y_a,F):#一般的な状況に対するRunge_Kutta法
    m = len(Y_a) #Y=[y_a1,..,y_am],F=[f1,..,fm],fi=fi(x,[y1,..,ym])
    Y = [[Y_a[j] for j in range(0,m)] for i in range(0,n+1)]
    for i in range(0,n):
        theta = coalitioned_rungekutta_theta(a+h*i,Y[i],h,F)
        for j in range(0,m):
            Y[i+1][j] = Y[i][j] + h*theta[j]
    return Y
def coalitioned_rungekutta_theta(x,Y,h,f):
    m = len(Y)
    k_1 = [f[i](x,Y) for i in range(0,m)]
    k_2 = [f[j](x+h/2,[Y[i]+(h/2)*k_1[i]for i in range(0,m)]) for j in range(0,m)]
    k_3 = [f[j](x+h/2,[Y[i]+(h/2)*k_2[i]for i in range(0,m)]) for j in range(0,m)]
    k_4 = [f[j](x+h,[Y[i]+h*k_3[i]for i in range(0,m)]) for j in range(0,m)]
    return [(1/6)*(k_1[i]+2*k_2[i]+2*k_3[i]+k_4[i])for i in range(0,m)]

def plt_chua(G,n,t_0,h,Y_0):#2次元平面におけるChua回路の(v_C_1,v_C_2)の描画
    def vC1(t,X):#変数Gの値に応じて微分方程式に対応する関数を定義し直す
        g = m_0*X[0]+(m_1-m_0)*abs(X[0]+B_p)/2+(m_0-m_1)*abs(X[0]-B_p)/2
        return (G*(X[1]-X[0])-g)/C_1
    def vC2(t,X):
        return (G*(X[0]-X[1])+X[2])/C_2
    def iL(t,X):
        return -X[1]/L
    res = coalitioned_rungekutta(n,t_0,h,Y_0,[vC1,vC2,iL])
    v1 = [res[i][0] for i in range(n//2,n+1)]#定常軌道に入るまでの情報は切り捨てる
    v2 = [res[i][1] for i in range(n//2,n+1)]
    plt.scatter(v1,v2,marker=".",s=1)
    plt.xlabel("v_1")
    plt.ylabel("v_2")
    plt.show()

def plt_chua_3d(G,n,t_0,h,Y_0):#3次元平面におけるChua回路の(v_C_1,v_C_2,i_L)の描画
    def vC1(t,X):
        g = m_0*X[0]+(m_1-m_0)*abs(X[0]+B_p)/2+(m_0-m_1)*abs(X[0]-B_p)/2
        return (G*(X[1]-X[0])-g)/C_1
    def vC2(t,X):
        return (G*(X[0]-X[1])+X[2])/C_2
    def iL(t,X):
        return -X[1]/L
    res = coalitioned_rungekutta(n,t_0,h,Y_0,[vC1,vC2,iL])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x = [res[i][0] for i in range(n//2,n+1)]
    y = [res[i][1] for i in range(n//2,n+1)]
    z = [res[i][2] for i in range(n//2,n+1)]
    ax.scatter(x, y, z, marker='.',s=1)
    ax.set_xlabel('v_1')
    ax.set_ylabel('v_2')
    ax.set_zlabel('i')
    plt.show()

def plt_chua_op(G,n,t_0,h,dt,Y_0):#時間遅延座標を用いたアトラクタの再構成
    def vC1(t,X):
        g = m_0*X[0]+(m_1-m_0)*abs(X[0]+B_p)/2+(m_0-m_1)*abs(X[0]-B_p)/2
        return (G*(X[1]-X[0])-g)/C_1
    def vC2(t,X):
        return (G*(X[0]-X[1])+X[2])/C_2
    def iL(t,X):
        return -X[1]/L
    dn = math.floor(n*dt)#時間遅延dtに応じて余分に数値計算する
    res = coalitioned_rungekutta(n+2*dn,t_0,h,Y_0,[vC1,vC2,iL])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x = [res[i][2] for i in range(n//2,n+1)]
    y = [res[i][2] for i in range(n//2+dn,n+1+dn)]
    z = [res[i][2] for i in range(n//2+2*dn,n+1+2*dn)]
    ax.scatter(x, y, z, marker='.',s=1)
    ax.set_xlabel('i(t)')
    ax.set_ylabel('i(t+dt)')
    ax.set_zlabel('i(t+2dt)')
    plt.show()

G = 0.66
n = 10000
a = 0
h = 0.01
T = [0.01,0.01,0.01]
plt_chua(G,n,a,h,T)
plt_chua_3d(G,n,a,h,T)
plt_chua_op(G,n,a,h,0.004,T)#dt = 0.004の場合にアトラクタはほぼ再構成された
