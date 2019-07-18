import os
import re
import sys
import sysconfig
import platform
import subprocess

from distutils.version import LooseVersion
from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
  def __init__(self, name, sourcedir=''):
    Extension.__init__(self, name, sources=[])
    self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
  def run(self):
    try:
      out = subprocess.check_output(['cmake', '--version'])
    except OSError:
      raise RuntimeError(
        "CMake must be installed to build the following extensions: " +
        ", ".join(e.name for e in self.extensions))

    if platform.system() == "Windows":
      cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)',
                                             out.decode()).group(1))
      if cmake_version < '3.1.0':
        raise RuntimeError("CMake >= 3.1.0 is required on Windows")

    pathtmp = self.build_temp
    for ext in self.extensions:
      pathdir,name0 = os.path.split(pathtmp)
      self.build_temp = os.path.join(pathdir,ext.name+"_"+name0)
      self.build_extension(ext)

  def build_extension(self, ext):
    extdir = os.path.abspath(
      os.path.dirname(self.get_ext_fullpath(ext.name)))
    cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + '.',
                  '-DPYTHON_EXECUTABLE=' + sys.executable]

    cfg = 'Debug' if self.debug else 'Release'
    build_args = ['--config', cfg]

    if platform.system() == "Windows":
      cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(
        cfg.upper(),
        extdir)]
      if sys.maxsize > 2 ** 32:
        cmake_args += ['-A', 'x64']
      build_args += ['--', '/m']
    else:
      cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
      build_args += ['--', '-j2']

    env = os.environ.copy()
    env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(
      env.get('CXXFLAGS', ''),
      self.distribution.get_version())
    if not os.path.exists(self.build_temp):
      os.makedirs(self.build_temp)
    print("$$$", self.build_temp)
    print("###cmake ext", ext.name, ext.sourcedir, extdir)
    print("@@@", cmake_args)
    subprocess.check_call(['cmake', ext.sourcedir] + cmake_args,
                          cwd=self.build_temp, env=env)
    subprocess.check_call(['cmake', '--build', '.'] + build_args,
                          cwd=self.build_temp)
    #subprocess.run(["rm",'-r',self.build_temp])
    print()  # Add an empty line for cleaner output


setup(name='PyDelFEM2',
      version='0.0.4',
      description='Handy toolset for implementing geometry processing and finite element simulation',
      author='Nobuyuki Umetani',
      author_email='n.umetani@gmail.com',
      url='https://github.com/nobuyuki83/delfem2',
      packages=['PyDelFEM2'],
      classifiers=['Topic :: Scientific/Engineering :: Physics',
                   'Topic :: Games/Entertainment :: Simulation',
                   'License :: OSI Approved :: MIT License'],
      keywords=['fem','simulation'],
      install_requires=[ 'numpy', 'PyOpenGL', 'glfw', 'PySide2' ],
      license="MIT",
      include_package_data=True,
      py_modules = ['PyDelFEM2','PyDelFEM2.gl','PyDelFEM2.eigen','PyDelFEM2.qt' ],      
      package_data = { '': ['*.so'], },      
      ext_modules=[CMakeExtension('c_core','src_pybind/core'),
                    CMakeExtension('c_gl','src_pybind/gl'),
                    CMakeExtension('c_eigen','src_pybind/eigen'),
                    ],  # location of *.so file, cmake file
      cmdclass=dict(build_ext=CMakeBuild),
      )