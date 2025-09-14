from setuptools import setup, Extension

module = Extension(
'symnmfmodule',
sources=['symnmfmodule.c', 'symnmf.c'],
extra_compile_args=['-O2', '-Wall', '-Wextra', '-Werror', '-pedantic-errors'],
extra_link_args=['-lm'], 
)

setup(
name='symnmfmodule',
version='1.0',
description='Symmetric NMF C extension module',
ext_modules=[module],
)