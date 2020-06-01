import os, setuptools

dir_path = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(dir_path, "requirements.txt")) as f:
    required_packages = f.read().splitlines()
with open(os.path.join(dir_path, "README.md"), "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="tsib",
    version="0.1.3",
    author="Leander Kotzur",
    author_email="l.kotzur@fz-juelich.de",
    description="Time Series Initialization for Buildings",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/FZJ-IEK3-VSA/tsib",
    include_package_data=True,
    packages=setuptools.find_packages(),
    install_requires=required_packages,
    setup_requires=["setuptools-git"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: End Users/Desktop",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    keywords=["buildings", "thermal load", "electricity load", "optimization"],
)
