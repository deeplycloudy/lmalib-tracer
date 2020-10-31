from setuptools import setup, find_packages

setup(
    name='lmalibtracer',
    version='0.1',
    description='LMA processing support software for the TRACER field program',
    packages=find_packages(),
    author='Eric Bruning',
    author_email='eric.bruning@gmail.com',
    url='https://github.com/deeplycloudy/lmalib-tracer/',
    license='BSD-3-Clause',
    long_description=open('README.md').read(),
    include_package_data=True,
)
