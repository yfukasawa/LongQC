import os
from setuptools import setup


def read(fname):
    try:
        with open(os.path.join(os.path.dirname(__file__), fname)) as f:
            return f.read()
    except:
        return ""

setup(
    name="mixem",
    version="0.1.3",

    author="Stefan Seemayer",
    author_email="stefan@seemayer.de",

    description="Expectation-Maximization (EM) algorithm for fitting mixtures of probability distributions",
    long_description=read("README.rst"),

    license="MIT",
    keywords="numeric em expectation maximization probability statistics distribution",
    url="https://github.com/sseemayer/mixem",

    packages=['mixem', 'mixem.distribution'],

    install_requires=[
        "numpy>=1.7.0",
        "scipy>=0.14.0",
    ],

    classifiers=[
        "Development Status :: 3 - Alpha",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: MIT License",

        "Intended Audience :: Science/Research",
        "Intended Audience :: Financial and Insurance Industry",

        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
)
