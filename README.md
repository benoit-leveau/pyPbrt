Overview
========

pyPbrt is a Python port/rewrite of pbrt, the physically based renderer by Matt Pharr and Greg Humphreys.

Original (C++) source code can be downloaded at http://pbrt.org/downloads.php

Coding Style
============

Whenever possible, I tried to adhere to "PEP 8 -- Style Guide for Python Code" for coding style and to "PEP 257 -- Docstring Conventions" for docstrings.

Links:
* PEP 8: http://www.python.org/dev/peps/pep-0008/
* PEP 257: http://www.python.org/dev/peps/pep-0257/

This significantly differs from the conventions adopted by the pbrt C++ source code.

                       |   pbrt                  |    pyPbrt
---------------------------------------------------------------------------------
Class Names            |   CapitalizedWords      |    CapitalizedWords
                       |   eg: Vector, Sphere    |    eg: Vector, Sphere
---------------------------------------------------------------------------------
Class Methods          |   CapitalizedWords      |    mixedCase
                       |   eg: vector.Length()   |    eg: vector.length()
---------------------------------------------------------------------------------
Global Variable Names  |                         |    lowercase_with_underscores
                       |   eg:                   |    eg:
---------------------------------------------------------------------------------
Function Names         |   CapitalizedWords      |    lowercase_with_underscores
                       |   eg: AbsDot()          |    eg: abs_dot()
---------------------------------------------------------------------------------
