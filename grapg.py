import numpy as np
import matplotlib.pyplot as plt

# Загрузка данных из .dat файла
data = np.loadtxt('data_file.dat')

# Извлечение координат x и y, а также значений z
x = data[0, 1:]
y = data[1:, 0]
z = data[1:, 1:]

# Создание сетки для координат x и y
X, Y = np.meshgrid(x, y)

# Построение контурных линий
plt.figure(figsize=(10, 8))
contour = plt.contour(X, Y, z, levels=150, cmap='viridis')
plt.clabel(contour, inline=True, fontsize=8)
plt.title('Контурные линии')
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar(contour, label='Уровень')

# Установка соотношения сторон
plt.gca().set_aspect('auto')  # Параметр 'auto' позволяет вручную задать соотношение
plt.gca().set_aspect(1)  # Увеличивает растяжение по Y относительно X
plt.xlim(0, X.max())
# plt.ylim(0, Y.max())
# plt.tight_layout(pad=2.0)

plt.grid()
plt.show()
