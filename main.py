import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from shapely.geometry import Polygon
import geopandas as gpd
import matplotlib as mpl

mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['axes.spines.top'] = False

# Load data
data_t15 = np.loadtxt("T15MD_imp_898_1000_135_05_al110-4.dat")
data_var1 = np.loadtxt("poloidal_Flux.txt")

r_t15 = data_t15[0, 1:]
r_var1 = data_var1[0, 1:]
z_t15 = data_t15[1:, 0]
z_var1 = data_var1[1:, 0]

data_t15 = data_t15[1:, 1:]
data_var1 = data_var1[1:, 1:]
delta_r_files = np.sum(np.abs(r_t15 - r_var1))
delta_z_files = np.sum(np.abs(z_t15 - z_var1))
F = data_t15
F2 = data_var1
R = r_t15
Z = z_t15
dr = R[1] - R[0]
dz = Z[1] - Z[0]

# Read camera data
camera_df = pd.read_excel('камера.xlsx', sheet_name='camera', header=None)
r_camera = camera_df.iloc[3:13, 1].to_numpy()
z_camera = camera_df.iloc[3:13, 2].to_numpy()

pf_coils_df = pd.read_excel('камера.xlsx', sheet_name='PF-coils', header=None)
pf_coils_cells = pd.read_excel('камера.xlsx', sheet_name='PF-coils', header=None)
inductor_df = pd.read_excel('камера.xlsx', sheet_name='Inductor', header=None)
r_inductor = np.concatenate((inductor_df.iloc[:, 1], inductor_df.iloc[:, 5], inductor_df.iloc[:, 9]))
z_inductor = np.concatenate((inductor_df.iloc[:, 2], inductor_df.iloc[:, 6], inductor_df.iloc[:, 10]))

probes_df = pd.read_excel('zondu.xlsx', header=None)
coils_df = pd.read_excel('камера.xlsx', sheet_name='coils', header=None)

probes = np.zeros((37, 2))
for i in range(37):
    for j in range(2):
        probes[i, j] = float(probes_df.iloc[i, j])

coils = np.zeros((11, 2))
for i in range(17, 28):
    for j in range(3, 5):
        coils[i - 17, j - 3] = float(coils_df.iloc[i, j])

pf_coils_df = pf_coils_df.drop(0,axis=1)
pf_coils_df = pf_coils_df.drop([0,1],axis=0)
# print(pf_coils_df.head())
r_pf_coils = pf_coils_df.iloc[:, 0].to_numpy()
z_pf_coils = pf_coils_df.iloc[:, 1].to_numpy()
dr_pf_coils = pf_coils_df.iloc[:, 3].to_numpy()
dz_pf_coils = pf_coils_df.iloc[:, 4].to_numpy()

r_div_plast = camera_df.iloc[3:7, 17].to_numpy()
z_div_plast = camera_df.iloc[3:7, 18].to_numpy()
d_div_plast = camera_df.iloc[3:7, 19].to_numpy()

r_camera_inner = camera_df.iloc[15:25, 1].to_numpy()
z_camera_inner = camera_df.iloc[15:25, 2].to_numpy()
r_camera_outer = camera_df.iloc[15:25, 3].to_numpy()
z_camera_outer = camera_df.iloc[15:25, 4].to_numpy()

RR, ZZ = np.meshgrid(R, Z)

# Read sensors data
sensors_df = pd.read_excel('sensors.xlsx', header=None)
sensors = np.zeros((30, 2))
for i in range(30):
    for j in range(2):
        sensors[i, j] = float(sensors_df.iloc[i, j])

# Read coils data
coils1_df = pd.read_excel('coils.xlsx', header=None)
coils1 = np.zeros((10, 2))
for i in range(10):
    for j in range(2):
        coils1[i, j] = float(coils1_df.iloc[i, j])

# Read diaphragm data
diafragma_df = pd.read_excel('difragma.xlsx', header=None)
diaf = np.zeros((20, 2))
for i in range(20):
    for j in range(2):
        diaf[i, j] = float(diafragma_df.iloc[i, j])

# Plotting
fig = plt.figure(figsize=(12,12))
plt.clf()
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')

# Plot camera polygon
camera_poly = Polygon(zip(np.concatenate((r_camera_inner, r_camera_outer[::-1])), 
                         np.concatenate((z_camera_inner, z_camera_outer[::-1]))))
camera_patch = gpd.GeoSeries([camera_poly])
camera_patch.plot(ax=plt.gca(), color='m')

# Plot PF coils
for i in range(len(r_pf_coils)):
    color = 'b' if i != 3 and i!=8 else 'g'
    rect = Rectangle((r_pf_coils[i] - dr_pf_coils[i]/2, z_pf_coils[i] - dz_pf_coils[i]/2), 
                     dr_pf_coils[i], dz_pf_coils[i], facecolor=color, ec='black')
    plt.gca().add_patch(rect)
    plt.text(pf_coils_cells.iloc[i+2, 1]+pf_coils_cells.iloc[i+2, 4]/2+0.01, pf_coils_cells.iloc[i+2, 2], pf_coils_cells.iloc[i+2, 0], fontsize=10, ha='left', weight='medium', va='center')
    
# Plot probes
for i in range(len(probes)):
    r_loc, z_loc = probes[i]
    rad = 0.02
    rect = Rectangle((r_loc-rad, z_loc-rad), 2*rad, 2*rad, facecolor='black', ec='black')
    plt.gca().add_patch(rect)

# Plot sensors
for i in range(len(sensors)):
    r_loc, z_loc = sensors[i]
    rad = 0.02
    rect = Rectangle((r_loc-rad, z_loc-rad), 2*rad, 2*rad, facecolor='yellow', ec='black')
    plt.gca().add_patch(rect)

# Plot diaphragm
r_d = diaf[:, 0]
z_d = diaf[:, 1]
plt.plot(r_d, z_d, 'black', linewidth=1.5)

# Plot additional elements
plt.gca().add_patch(Rectangle((0.776589-0.01, 0.880668), 0.02, 1.08082 - 0.880668, facecolor='c', ec='black'))
plt.gca().add_patch(Rectangle((0.776589-0.01, -1.0568), 0.02, abs(-0.848644 - -1.03278), facecolor='c', ec='black'))
pgon = Polygon([(1.8259-0.01+0.01, -1.20091+0.01), (1.93747-0.01, -1.03278), 
                (1.93747+0.01-0.01, -1.03278-0.01), (1.8259+0.01, -1.20091)])
pgon_patch = gpd.GeoSeries([pgon])
pgon_patch.plot(ax=plt.gca(), color='c', ec='black', lw=0.8)

pgon = Polygon([(1.93747-0.01, 1.04079), (1.8259-0.01+0.01, 1.19291-0.01), 
                (1.8259+0.01, 1.19291), (1.93747+0.01-0.01, 1.04079+0.01)])
pgon_patch = gpd.GeoSeries([pgon])
pgon_patch.plot(ax=plt.gca(), color='c', ec='black', lw=0.8)

# Plot div_plast
plt.gca().add_patch(Rectangle((r_div_plast[0], z_div_plast[0]-d_div_plast[0]/2), 
                              r_div_plast[1]-r_div_plast[0], d_div_plast[0], facecolor=[0.7, 0.7, 0.7], ec='black'))
plt.gca().add_patch(Rectangle((r_div_plast[2], z_div_plast[2]-d_div_plast[2]/2), 
                              r_div_plast[3]-r_div_plast[2], d_div_plast[2], facecolor=[0.7, 0.7, 0.7], ec='black'))

sep = 0.0285771
x = -0.53256

levels = np.linspace(x, 13.5 * x, 75)
levels = np.sort(levels)  # Ensure levels are sorted in increasing order

plt.contour(RR, ZZ, F, levels, linewidths=0.2)
plt.contour(RR, ZZ, F, [x], colors='g', linewidths=0.2)
# plt.contour(RR, ZZ, F2, [sep], colors='r')

# GRAPH FIELD

# Загрузка данных из .dat файла
data = np.loadtxt('T15MD_imp_898_1000_135_05_al110-4.dat')

# Извлечение координат x и y, а также значений z
x = data[0, 1:]
y = data[1:, 0]
z = data[1:, 1:]

# Создание сетки для координат x и y
X, Y = np.meshgrid(x, y)

# print(np.min(z), np.max(z))
# print(X,Y,z)


# Установим порог
special_level = -0.24811

# Определяем уровни для контуров, ниже порога
levels = np.linspace(-1, special_level, 19)[:-1]


contour = plt.contour(X, Y, z, levels=levels, colors='y', linestyles='solid', linewidths=0.5)
special_contour = plt.contour(X, Y, z, levels=[special_level], colors='g', linestyles='solid', linewidths=0.5)
# plt.colorbar(contour)

# plt.clabel(contour, inline=True, fontsize=8)
plt.title('T-15MD')
plt.xlabel('R')
plt.ylabel('Z')
plt.xticks(np.arange(0, 4, 0.5))
plt.yticks(np.arange(-3.5, 3.0, 0.5))
plt.ylim(-3.28 - 0.01, 2.411 + 0.01)
plt.xlim(0.26 - 0.01, 3.32 + 0.01)
# bag: special_contour_ax = special_contour.collections[0].set_label('-0.24811 inverse') don't show green line in legend because plot fake line below
plt.plot([5, 6], [5, 8], label='-0.24811 inverse', color='g', linewidth=0.5)
plt.legend(loc=(1.25, 1))
# special_contour_ax.legend()
# plt.colorbar(contour, label='Уровень')
# Создание сетки для координат x и y
# X, Y = np.meshgrid(x, y)

plt.show()
