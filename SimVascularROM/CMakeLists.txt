#-----------------------------------------------------------------------------
set(MODULE_NAME SimVascularROM)

#-----------------------------------------------------------------------------
# `svromsetup` does the work; this module is an MRML adapter over it. Listing the package's files
# here installs them beside the module, which is already on Slicer's Python path, as Mesh Prep
# does with `svmeshcomplete`. `sv_rom_simulation` is not listed: it is a package of its own,
# found where it is installed or where the panel's Tools section says a checkout is.
set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  svromsetup/__init__.py
  svromsetup/boundary_conditions.py
  svromsetup/case.py
  svromsetup/results.py
  svromsetup/solver.py
  svromsetup/testing.py
  svromsetup/visualization.py
  )

set(MODULE_PYTHON_RESOURCES
  Resources/UI/${MODULE_NAME}.ui
  LICENSE.txt
  )

#-----------------------------------------------------------------------------
slicerMacroBuildScriptedModule(
  NAME ${MODULE_NAME}
  SCRIPTS ${MODULE_PYTHON_SCRIPTS}
  RESOURCES ${MODULE_PYTHON_RESOURCES}
  WITH_GENERIC_TESTS
  )

#-----------------------------------------------------------------------------
if(BUILD_TESTING)

  # Register the unittest subclass in the main script as a ctest.
  # Note that the test will also be available at runtime.
  slicer_add_python_unittest(SCRIPT ${MODULE_NAME}.py)

  # Additional build-time testing
  add_subdirectory(Testing)
endif()
