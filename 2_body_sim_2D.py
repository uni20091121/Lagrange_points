import numpy as np
import matplotlib.pyplot as plt
import spiceypy as spi
import heapq


G=6.67e-11
M_1=6e24
M_2=7e22
d=3.84e8

r_1_0=np.array([0,0,0])
r_2_0=np.array([d,0,0])

L=1.5*d


plt.style.use('dark_background')

fig=plt.figure()
ax=fig.add_subplot()

ax.set_xlim(-L,L)
ax.set_ylim(-L,L)
ax.set_aspect('equal')


Body_1=plt.Circle(
    (r_1_0[0],r_1_0[1]),
    6.3e6,
    color="#0627B8"
)

Body_2=plt.Circle(
    (r_2_0[0],r_2_0[1]),
    6.3e6/6,
    color="#B0B0B0"
)


x=np.linspace(-L,L,100)
y=np.linspace(-L,L,100)

X,Y=np.meshgrid(x,y)


ax.add_artist(Body_1)
ax.add_artist(Body_2)


r_com=(M_1*r_1_0+M_2*r_2_0)/(M_1+M_2)


r_1=np.sqrt(
    (X-r_1_0[0])**2+
    (Y-r_1_0[1])**2
)

r_2=np.sqrt(
    (X-r_2_0[0])**2+
    (Y-r_2_0[1])**2
)


k=np.sqrt(
    (r_2_0[0]-r_com[0])**2+
    (r_2_0[1]-r_com[1])**2
)


omega=np.sqrt(
    G*(M_1+M_2)/(k)**3
)


V=(
    -G*M_1/r_1
    -G*M_2/r_2
    -0.5*omega**2*
    ((X-r_com[0])**2+(Y-r_com[1])**2)
)


do_V_x=(
    -G*M_1/r_1**3*(X-r_1_0[0])
    -G*M_2/r_2**3*(X-r_2_0[0])
    -omega**2*(X-r_com[0])
)


do_V_y=(
    -G*M_1/r_1**3*(Y-r_1_0[1])
    -G*M_2/r_2**3*(Y-r_com[1])
    -omega**2*(Y-r_com[1])
)


do_V_r=np.sqrt(do_V_x**2+do_V_y**2)


del_V=[]

for do in do_V_r:
    del_V.append(do)

del_V=np.array(del_V)


V_scaled=V/(G*M_1/d)

X_scaled=X/d
Y_scaled=Y/d


levels=np.linspace(
    np.min(V),
    np.max(V),
    15000
)


c=ax.contourf(
    X,
    Y,
    V,
    levels=levels,
    cmap="inferno",
    alpha=0.5
)

plt.colorbar(c)
plt.title('A potential map of the 2-body system')

plt.axis('off')


plt.show()