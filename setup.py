from setuptools import setup, Extension

module = Extension('mykmeanssp',
                    sources=[
                        'kmeansmodule.c', 
                        'kmeans.c'
                        ])

setup(
    name='mykmeanssp',
    version='1.0',
    description='KMeans C extension module',
    ext_modules=[module]
)
