# Welcome to pythroughput!

![Python-version](https://img.shields.io/badge/Python-3.5-green.svg)
![BSD-3-LICENSE](https://img.shields.io/badge/licence-BSD--3--Clause-blue.svg)

Python module (component) to perform high-throughput first-principles 
calculation in [XenonPy](http://github.com/yoshida-lab/XenonPy) package.

## Requirement

### python

[XenonPy](https://github.com/yoshida-lab/XenonPy) requires
[python](https://www.python.org/) 3.5 or later.
Therefore, you should prepare python with that version.

### pymatgen

In our module, we use [pymatgen](https://github.com/materialsproject/pymatgen)
package in order to handle crystal structures on your computer.

If you have pip, you can install pymatgen as follows:

~~~~
> pip install --user pymatgen
~~~~

### First-principles calculation package

Our module (will) support calculation packages of GPAW and VASP.
GPAW released under the GNU General Public License as published by the Free Software Foundation,
so you can use it freely for any purposes with no licence fee.
However, at present, VASP has higher performance and analysis software is more substantial.

Now, we are implementing the option to perform high-throughput calculation using VASP,
and we intend to be able to select this when you want higher performances.

In order to use GPAW, you must install ASE, Libxc,
BLAS/LAPACK and some high performance calculation package in python.
(If you installed python via Anaconda or Homebrew
packages, you should already have it!)

For more details, see Installation section in
[GPAW wiki](https://wiki.fysik.dtu.dk/gpaw/install.html).

In convenience, we show the link to the packages
you need to perform first-principles calculation in our module.

**GPAW**

1. [GPAW 1.4.0 or later](https://wiki.fysik.dtu.dk/gpaw/install.html)
2. [ASE 3.16.0 or later](https://wiki.fysik.dtu.dk/ase/)
3. [Libxc 2.0.1 or later](http://www.tddft.org/programs/libxc/)

## Usage

In preparation. (See samples, please.)

## Install

Clone this repository.

## Licence

This module will be released under the BSD 3-Clause License.

~~~~
Copyright (c) 2018-2019, Taku MURAKAMI
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
~~~~

## Author

Taku MURAKAMI[@murakami17](https://github.com/murakami17/),
master course student at Shizuoka university.

Please contant me via:
[e-mail](mailto:murakami.taku.17@shizuoka.ac.jp) or
[github issues](https://github.com/murakami17/pythoughput/issues).
