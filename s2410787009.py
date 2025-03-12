import math
import argparse
import pandas as pd
import scipy.constants as sc
import numpy as np
import matplotlib.pyplot as plt

# For Arguments from users
parser = argparse.ArgumentParser()
parser.add_argument('weight', type=float,nargs='+')
parser.add_argument('slip', type=float)
parser.add_argument('mu', type=float)
args = parser.parse_args()

# Declaration of inputs
mass=args.weight
alpha=args.slip
mu=args.mu

# Reading constants from csv
table = pd.read_csv('resources/table_Code.csv', index_col=0)
a1_x=table.loc['a1', 'Fx']
a2_x=table.loc['a2', 'Fx']
a3_x=table.loc['a3', 'Fx']
a4_x=table.loc['a4', 'Fx']
a5_x=table.loc['a5', 'Fx']
a6_x=table.loc['a6', 'Fx']
a7_x=table.loc['a7', 'Fx']
a8_x=table.loc['a8', 'Fx']

# Declaring constants required for Brake force
C=1.65
g=sc.g
phi=sc.pi

# Calculation and plot of Brake force vs slip
for m in mass:
    Fz = (m * g) / 4000
    Dx = a1_x * Fz**2 + a2_x*Fz
    Ex = a6_x * Fz**2 + a7_x * Fz + a8_x
    Bx = (a3_x * Fz**2 + a4_x * Fz) / (C * Dx * math.exp(a5_x * Fz))

    Fx_Int = []
    k_Int = []

    for i in range(0, 101):
        phi_value = (1-Ex)*i+(Ex/Bx)*math.atan(Bx*i)*180/phi
        Sigma_x=-i/(1+i)
        Sigma_y=-math.tan(alpha*phi/180)/(1+i)
        Overall_Sigma=math.sqrt(Sigma_x**2+Sigma_y**2)
        fx = Dx * math.sin(C * math.atan(Bx * phi_value))
        fx = -(Sigma_x/Overall_Sigma)*mu*fx

        Fx_Int.append(fx)
        k_Int.append(i)

    plt.plot(k_Int, Fx_Int, label=f'Fy = {Fz}')
plt.xlabel('Longitudinal Slip k[%]')
plt.title('Side force and brake force vs longitudinal')
plt.legend()
plt.grid(True)

# Reading constants required for side force from csv
a1_y=table.loc['a1', 'Fy']
a2_y=table.loc['a2', 'Fy']
a3_y=table.loc['a3', 'Fy']
a4_y=table.loc['a4', 'Fy']
a5_y=table.loc['a5', 'Fy']
a6_y=table.loc['a6', 'Fy']
a7_y=table.loc['a7', 'Fy']
a8_y=table.loc['a8', 'Fy']

# Declaring constants required for Side force
CY = 1.30
Y=0

# Calculation and plot of Side force vs Slip
for m in mass:
    Fz = (m * g) / 4000
    Dy = a1_y * Fz**2 + a2_y * Fz
    Ey = a6_y * Fz**2 + a7_y * Fz + a8_y
    By= (a3_y*np.sin(a4_y*math.atan(a5_y*Fz)))/(CY*Dy)
    Fy_Int = []
    k_Int = []

    for i in range(0,101):
        Sigma_x=-i/(1+i)
        Sigma_y=-math.tan(alpha*phi/180)/(1+i)
        Overall_Sigma=math.sqrt(Sigma_x**2+Sigma_y**2)
        phi_y=(1-Ey)*alpha+(Ey/By)*math.atan(By*(alpha*phi/180))
        Fy = Dy * np.sin(CY * math.atan(By * phi_y))
        Fy = -(Sigma_y/Overall_Sigma)*mu*Fy
        Fy_Int.append(Fy)
        k_Int.append(i)

    plt.plot(k_Int, Fy_Int, label=f'alpha = {alpha}', linestyle=':')
plt.show()
