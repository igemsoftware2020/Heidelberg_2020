import setuptools 

with open("README.md", "r") as fh:
  long_description = fh.read() 

setuptools.setup(
  name="3DOC",
  version="0.0.1",
  author="Lucas Arnoldt",
  author_email="lucas.arnoldt@stud.uni-heidelberg.de",
  description="Pipeline for PDB generation of fusion protein domains.",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/igemsoftware2020/Heidelberg_2020.git",
  packages=setuptools.find_packages(),
  classifiers=("Programming Language :: Python :: 3", "License :: OSI Approved :: MIT License", "Operating System :: OS Independent",
  ),
)
