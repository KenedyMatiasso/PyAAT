from setuptools import setup, find_packages
import pathlib

HERE = pathlib.Path(__file__).parent

VERSION = "0.0.dev10"
PACKAGE_NAME = "PyAAT"
AUTHOR = "Kenedy Matiassp Portella"
AUTHOR_EMAIL = "kenedyportella@hotmail.com"
URL = "https://github.com/KenedyMatiasso/PyAAT"

LICENSE = "MIT"
DESCRIPTION = "Python Aerospace Analysis Toolbox"
LONG_DESCRIPTION = (HERE / "docs/LongDescription.MD").read_text()
LONG_DESC_TYPE = "text/markdown"

INSTALL_REQUIRES = [
	"numpy",
        "scipy",
        "matplotlib",
        "control"
]

setup(name=PACKAGE_NAME,
      version=VERSION,
      description=DESCRIPTION,
      long_description=LONG_DESCRIPTION,
      long_description_content_type=LONG_DESC_TYPE,
      author=AUTHOR,
      license=LICENSE,
      author_email=AUTHOR_EMAIL,
      url=URL,
      install_requires=INSTALL_REQUIRES,
      packages=find_packages()
      )
