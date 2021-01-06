from setuptools import setup


# Read in requirements.txt
requirements = open('requirements.txt').readlines()
requirements = [r.strip() for r in requirements]

setup(name='jwcas',
      version='0.1',
      description='This package does simple Jordon Wigner transformation of the provided fermionic Hamiltonian',
      url='http://github.com/vibinabraham/jwcas',
      author='Vibin Abraham',
      author_email='vibin1@vt.edu',
      packages=['jwcas'],
      install_requires=requirements,
      license='MIT',
      zip_safe=False)
