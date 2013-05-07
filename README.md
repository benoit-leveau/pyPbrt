Overview
========

pyPbrt is a Python port/rewrite of <b>pbrt</b>, the physically based renderer by Matt Pharr and Greg Humphreys, which is the companion software for their book <i>Physically Based Rendering, Second Edition: From Theory To Implementation</i>.
(<a href="http://www.amazon.com/Physically-Based-Rendering-Second-Implementation/dp/0123750792">amazon page for the book</a>)

The pbrt website (for the book and the software) is located at http://pbrt.org .
The original (C++) source code can be downloaded at http://pbrt.org/downloads.php


Rationale
=========

The idea behind this project is mainly for me to go along the book and write code as I progress through it.
Having a raytracer or path tracer entirely written in Python doesn't make much sense as far as performance is concerned,
but I thought it would be a good way to make sure I understand all the algorithms as I have to write them. Writing it in C++
would make copy-paste from the book (or from the source code) too tempting, and inevitably I would come up with the same exact code
as pbrt.
There is one benefit to having it written in Python though. It becomes very easy to test new algorithms and concepts interactively in the
interpreter, reloading code when making changes and testing them directly.

Parts of pbrt that are closely tied to the C++ implementation, like the MemoryArena class, don't make sense in Python and will therefore be skipped.

A more clever approach, if this project wasn't meant for learning purposes, would be to have a Python bindings to a C++ implementation of pbrt, something that can even be 
automated with a framework like py++ (http://sourceforge.net/projects/pygccxml/)

Once the project is finished, I'll also try to see if I can get decent running times by switching the low-level parts (geometry and transformation classes) to C++, or by using PyCuda.


Installation
============

Just cd to the project folder and run <i>pyPbrt</i>.
Use the flag <i>--help</i> for a description of all available parameters.

Unit Tests
==========

Unit tests are implemented using the <i>unittest</i> module and are located in the tests folder.
For convenience, I use <i>nose</i> to run the tests. More information about nose can be found there: https://pypi.python.org/pypi/nose/1.3.0
Once nose is installed, running the tests is just a matter of cd'ing to the pyPbrt folder and running:
<pre>
[~/Documents/Prog/pyPbrt] - nosetests
----------------------------------------------------------------------
Ran X test in 0.018s

OK
</pre>

Coding Style
============

Whenever possible, I tried to adhere to <i>PEP 8 -- Style Guide for Python Code</i> for coding style and to <i>PEP 257 -- Docstring Conventions</i> for docstrings.
Code is checked with the command line tools <i>pep8</i> and <i>pylint</i>.

Links:
* PEP 8: http://www.python.org/dev/peps/pep-0008/
* PEP 257: http://www.python.org/dev/peps/pep-0257/
* pep8 command line tool: https://pypi.python.org/pypi/pep8/
* pylint command line tool: http://www.logilab.org/857

This significantly differs from the conventions adopted by the pbrt C++ source code.
The differences are highlighted in the table below:

<table>
  <tr>
    <th></th><th>pbrt</th><th>pyPbrt</th>
  </tr>
  <tr>
    <td>Global Variable Names</td><td></td><td>lowercase_with_underscores</td>
  </tr>
  <tr>
    <td>Function Names</td><td>CapWords<br>eg: <i>AbsDot(v1, v2)</i></td><td>lowercase_with_underscores<br>eg: <i>abs_dot(v1, v2)</i></td>
  </tr>
  <tr>
    <td>Class Names</td><td>CapWords<br>eg: <i>Vector, Sphere</i></td><td>CapWords<br>eg: <i>Vector, Sphere</i></td>
  </tr>
  <tr>
    <td>Class Methods</td><td>CapWords<br>eg: <i>vector.LengthSquared()</i></td><td>lowercase_with_underscores<br>eg: <i>vector.length_squared()</i></td>
  </tr>
</table>
