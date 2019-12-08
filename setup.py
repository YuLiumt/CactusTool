from distutils.core import setup

setup(
    name="CactusTool", 
    version="0.0.1",
    description = "Cactus Tool",
    author = "Yu Liu",
    author_email = "yuliumutian@gmail.com",
    url = "https://cactustool.readthedocs.io/en/latest/",
    packages=['CactusTool'],
    install_requires=[
        "numpy", 
        "scipy", 
        "h5py",
        "matplotlib",
    ])