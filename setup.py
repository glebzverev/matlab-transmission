import sys
from cx_Freeze import setup, Executable

# Определение базовой опции для Windows
base = None
if sys.platform == "win32":
    base = "Win32GUI"  # для GUI-приложений, используйте "Win32GUI"

# Определение целевого исполняемого файла
executables = [
    Executable("main.py", base=base, target_name="MyApp")
]

# Настройки сборки
setup(
    name="MyApp",
    version="0.1",
    description="My cross-platform app",
    executables=executables
)
