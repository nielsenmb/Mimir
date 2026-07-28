Installation
============

Published releases
------------------

Install Mimir from PyPI:

.. code-block:: bash

   python -m pip install mimir-astro

The distribution is named ``mimir-astro`` because an unrelated project already
uses ``mimir`` on PyPI. The Python import remains:

.. code-block:: python

   import mimir

MAST support
------------

Lightkurve-based MAST access is an optional dependency:

.. code-block:: bash

   python -m pip install "mimir-astro[mast]"

Source and development installations
------------------------------------

Until the first PyPI release, install the latest public source directly:

.. code-block:: bash

   python -m pip install "mimir-astro @ git+https://github.com/nielsenmb/Mimir.git"

For development, clone the repository and install all checking dependencies:

.. code-block:: bash

   git clone https://github.com/nielsenmb/Mimir.git
   cd Mimir
   python -m pip install -e ".[test,mast,docs]"
