from setuptools import setup

setup(name='pyTDFSDK',
      version='0.1.6',
      description='Python wrapper for Bruker TDF-SDK',
      url='https://github.com/gtluu/pyTDFSDK',
      author='Gordon T. Luu',
      author_email='gtluu912@gmail.com',
      license='Apache License',
      packages=['pyTDFSDK', 'TDF-SDK'],
      include_package_data=True,
      package_data={'': ['*.dll', '*.so']},
      install_requires=['numpy'])
