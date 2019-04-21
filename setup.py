import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="Lindel",
    version="0.0.1",
    author="Wei Chen",
    author_email="wchen108@uw.edu",
    description="A package for Lindel",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/shendurelab/CRISPR_NHEJ_prediction/tree/master/scripts",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)