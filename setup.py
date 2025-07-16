from setuptools import setup, find_packages

setup(
    name="fastacheck",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "colorama>=0.4.4",
        "numpy>=1.21.0",
        "scipy>=1.7.0",
        "tqdm>=4.62.0",
    ],
    entry_points={
        "console_scripts": [
            "fastacheck=fastacheck.main:main",
        ],
    },
    python_requires=">=3.8",
)