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

    You can activate your virtual env from ``ccd_phot`` dir

    .. code-block:: console

      .\venv\Scripts\activate.bat

    and than start any script

    .. code-block:: console

      python script_name.py arguments
      python sat_phot.py arguments


    This is enough or,
    if you need to start your virtual env from any other dir, you can make bat file (example ``ccd_python.bat`` with code below).
    Change ``PROJECT_ROOT`` to find your ccd_phot directory. Add ``ccd_phot`` directory to your system path (see: env_link_).
    Then, in any directory, you can type:

    .. code-block:: console

        ccd_python script_name.py arguments

    For example:

    .. code-block:: console

        ccd_python sat_phot.py . -c config_name.ini


    Code of ``ccd_python.bat``:

    .. code-block:: text

      @echo off

      :: Define your project root path
      set "PROJECT_ROOT=path\to\ccd_phot"
      set "PYTHON_EXE=%PROJECT_ROOT%\venv\Scripts\python.exe"
  
      :: Check if script name is present (original %1)
      if "%1"=="" (
         echo Error: No Python script name was specified.
         goto :eof
      )
  
      :: Construct the full script path using the initial %1
      set "SCRIPT_PATH=%PROJECT_ROOT%\%1"
  
      :: --- Validation Checks ---
      if not exist "%PYTHON_EXE%" (
         echo Error: Python interpreter not found at: "%PYTHON_EXE%"
         goto :eof
      )
      if not exist "%SCRIPT_PATH%" (
         echo Error: Script not found at path: "%SCRIPT_PATH%"
         goto :eof
      )
  
  
      setlocal ENABLEDELAYEDEXPANSION
       set "_args=%*"
       set "_args=!_args:*%1 =!"
  
  
      "%PYTHON_EXE%" "%SCRIPT_PATH%" %_args%
  
      endlocal


Instruction for **Linux users**:
----------------------------------

Help yourself in a way you like the most :).

You know how to deal with it.

.. _env_link: https://www.eukhost.com/kb/how-to-add-to-the-path-on-windows-10-and-windows-11/
