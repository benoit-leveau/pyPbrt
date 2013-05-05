Overview
========

pyPbrt is a Python port/rewrite of <b>pbrt</b>, the physically based renderer by Matt Pharr and Greg Humphreys, which is the companion software for their book <i>Physically Based Rendering, Second Edition: From Theory To Implementation</i>.
(<a href="http://www.amazon.com/Physically-Based-Rendering-Second-Implementation/dp/0123750792">amazon page for the book</a>)

The pbrt website (for the book and the software) is located at http://pbrt.org .
The original (C++) source code can be downloaded at http://pbrt.org/downloads.php

Coding Style
============

Whenever possible, I tried to adhere to <i>PEP 8 -- Style Guide for Python Code</i> for coding style and to <i>PEP 257 -- Docstring Conventions</i> for docstrings.

Links:
* PEP 8: http://www.python.org/dev/peps/pep-0008/
* PEP 257: http://www.python.org/dev/peps/pep-0257/

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
    <td>Class Methods</td><td>CapWords<br>eg: <i>vector.Length()</i></td><td>mixedCase<br>eg: <i>vector.length()</i></td>
  </tr>
</table>
