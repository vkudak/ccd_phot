.. _start_scripts:

How to start scripts
====================

There are few ways to start scripts from this project and begin processing FITS files.

Instruction for **Windows users**:
----------------------------------

#. Start script from ``CCD_PHOT`` directory (install dir).
    If there is only one Python environment in you PC ecosystem.
    That means you have installed Python only for processing photometry files and will not install any other Python packages.


    Navigate to the installation folder ``CCD_PHOT`` and start any script with expected arguments. Example:

    .. code-block:: console

     python3 sat_phot.py path\to\fits\files

    You can also add ``CCD_PHOT`` directory to your system PATH (how to: env_link_) and use scripts from any directory. For example if you are located in directory with your FITS file, than just do:

    .. code-block:: console

       python3 sat_phot.py .

    Where ``.`` means - use current working directory

#. Use virtual environment.
    If you create and setup virtual environment, as it described in :ref:`install section <installation>`, you can use it to lunch the scripts.
    In this case you will never had experience of Python packages version conflicts.

    For example you call your virtual environment  ``venv`` as it proposed in :ref:`install section <installation>`.

    You can make executable `bat` file (for example name it ``ccd_phot_python.bat``) that will use Python from your ``venv``.
    This file should be created in ``CCD_PHOT`` folder where also ``venv`` directory is located.
    Command inside ``ccd_phot_python.bat`` file should look like:

    .. code-block:: console

      venv\Scripts\python %*

    And now you should be able to start any script with command:

    .. code-block:: console

      ccd_phot_python.bat script_name.py arguments

    Example:

    .. code-block:: console

      ccd_phot_python.bat sat_phot.py path/to/fits/files and any other arguments if needed


Instruction for **Linux users**:
----------------------------------

Help yourself in a way you like the most :)

.. _env_link: https://www.eukhost.com/kb/how-to-add-to-the-path-on-windows-10-and-windows-11/
