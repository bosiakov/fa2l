from setuptools import setup

setup(
    name="fa2l",
    version="0.1",
    author="Eugene Bosiakov",
    author_email="eugenebosyakov@gmail.com",
    description="Force Atlas 2 graph layout",
    packages=['fa2l'],
    install_requires=[
        "networkx<2.0.0",
        "numpy"
    ],
)