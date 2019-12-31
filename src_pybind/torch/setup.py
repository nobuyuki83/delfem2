from setuptools import setup, Extension
from torch.utils import cpp_extension

setup(name='torch_delfem2',
      ext_modules=[cpp_extension.CppExtension('torch_delfem2', ['main.cpp'])],
      cmdclass={'build_ext': cpp_extension.BuildExtension})