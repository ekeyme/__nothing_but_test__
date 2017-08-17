from setuptools import setup

setup(
    name = "bio-pm",
    version = '1.1.0b5',
    url = "https://github.com/ekeyme/bio-pm",
    license = "MIT",
    author = "Ekeyme Mo",
    author_email = "ekeyme@gmail.com",
    description = "point mutation pattern analyzing tool for nucleotide sequence",
    keywords = 'python3 point mutation analyzing tool',
    packages = ["pm", ],
    install_requires = ['biopython', ],
    classifiers= [
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3 :: Only",
    ]
)