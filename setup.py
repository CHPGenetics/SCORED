from setuptools import setup, find_packages

setup(
    name="SCORED",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A package for single-cell RNA sequencing data imputation using SimRank and Random Walk with Restart.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/scored",
    packages=find_packages(),
    install_requires=[
        "scikit-learn>=1.0",
        "scipy>=1.7",
        "matplotlib>=3.4",
        "seaborn>=0.11",
        "tqdm>=4.0",
        "networkx>=2.6",
        "numpy>=1.20",
        "torch>=1.9",
        "scanpy>=1.8",
        "pandas>=1.3",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)