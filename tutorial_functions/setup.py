from setuptools import setup

setup(
    name="tutorial_functions",
    version="0.1.0",
    description="Tutorials for unfolding",
    url="https://github.com/bjorneju/pyphi_units",
    author="Bjørn E juel",
    author_email="bjorneju@gmail.com",
    license="BSD 2-clause",
    packages=["tutorial_functions"],
    install_requires=[
        "tqdm",
        "numpy",
        "pandas",
        "pyphi",
        "matplotlib",
    ],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.9",
    ],
)
