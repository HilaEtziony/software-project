from setuptools import setup, Extension

module = Extension('mykmeanspp',
                    sources=[
                        'kmeansmodule.c', 
                        'kmeans.c'
                        ])

setup(
    name='mykmeanspp',
    version='1.0',
    description='KMeans C extension module',
    ext_modules=[module]
)
