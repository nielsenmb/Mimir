Installation
============

Development installation
------------------------

Install the package and its test dependencies from a checkout:

.. code-block:: bash

   python -m pip install -e ".[test]"

MAST support
------------

Lightkurve-based MAST access is an optional dependency:

.. code-block:: bash

   python -m pip install -e ".[mast]"

The distribution will be published as ``mimir-astro`` while the Python import
name is ``mimir``.
