from setuptools import setup

with open('README.rst', 'r', 'utf-8') as f:
    readme = f.read()

setup(
    name="fa2l",
    version="0.1",
    author="Eugene Bosiakov",
    author_email="eugenebosyakov@gmail.com",
    url='https://github.com/bosiakov/fa2l',
    description="Force Atlas 2 graph layout",
    long_description=readme,
    packages=['fa2l'],
    install_requires=[
        "networkx<2.0.0",
        "numpy"
    ],
)
