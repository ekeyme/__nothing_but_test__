from setuptools import setup
from pm import __version__ 

setup(
    name = "biopm",
    version = __version__,
    url = "https://github.com/ekeyme/bio-pm",
    license = "MIT",
    author = "Ekeyme Mo",
    author_email = "ekeyme@gmail.com",
    description = "point mutation pattern analyzing tool for nucleotide sequence",
    packages = ["pm", ],
    platforms = "any",
    install_requires = ['biopython', ],
    classifiers= [
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3 :: Only",
    ]
)