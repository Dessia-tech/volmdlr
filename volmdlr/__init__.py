# -*- coding: utf-8 -*-


import pkg_resources
__version__ = pkg_resources.require("volmdlr")[0].version

from .core import *
