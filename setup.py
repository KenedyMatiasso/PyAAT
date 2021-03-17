from setuptools import setup, find_packages

setup(
    name="PyAAT",
    description="Python Aerospace Analysis Toolbox",
    version="0.0.dev1",
    author="Kenedy Matiassp Portella",
    author_email="kenedyportella@hotmail.com",
    download_url="https://github.com/KenedyMatiasso/PyAAT",
    license="MIT",
    packages=find_packages("pyaat"),
    keywords=[
      "aerospace", "aeronauticd", "dynamics", "flight mechanics",
      "simulation", "space systems", "aero", "control"
    ],
    python_requires=">=3.5",
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "control"
    ],
)
