import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="zaya",
    version="0.1",
    author="Thomas Titscher",
    author_email="thomas.titscher@gmail.com",
    description="Some helper scripts and classes",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=["tabulate", "pyyaml"],
    extras_require={  # Optional
        "dev": ["black"],
        "test": ["pytest", "pytest-cov", "flake8"],
    },
    entry_points={"console_scripts": ["texify=zaya.texify:main"]}
)
