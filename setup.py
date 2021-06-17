import sys
from setuptools import setup

if "clean" in sys.argv[1:]:
    pass
else:
    try:
        import fenics # noqa
    except ImportError:
        raise Exception("FEniCS needs to be installed and activated. "
                        "Please visit fenicsproject.org")
setup(
    name='IRKsome',
    version='0.0.1',
    author='Rob Kirby, Jorge Marchena Menendez',
    author_email='Robert_Kirby@baylor.edu',
    description='A library for fully implicit Runge-Kutta methods in FEniCS (port of the Firedrake version of IRKsome)',
    long_description='',
    packages=['irksome'],
    zip_safe=False,
)
