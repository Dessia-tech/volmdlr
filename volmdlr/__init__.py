# -*- coding: utf-8 -*-


import pkg_resources
__version__ = pkg_resources.require("genmechanics")[0].version

from .core import *
