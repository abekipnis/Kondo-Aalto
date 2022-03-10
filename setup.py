from setuptools import setup
from os.path import join, dirname

requirementstxt = join(dirname(__file__), "requirements.txt")
requirements = [ line.strip() for line in open(requirementstxt, "r") if line.strip() ]

setup(
    name='AaltoAtoms',
    version='0.0.0',
    packages=['AaltoAtoms'],
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "Operating System :: Microsoft :: Windows",
        "Topic :: Scientific/Engineering :: Data Analysis",

    ])
