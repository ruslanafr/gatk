from distutils.core import setup

setup(
    name='gcnv',
    version='0.2dev',
    author='Mehrtash Babadi',
    author_email='mehrtash@broadinstitute.org',
    packages=['gcnv'],
    license='LICENSE.txt',
    description='GATK gCNV-theano computational core',
    long_description=open('README.txt').read(),
    install_requires=[
        "theano >= 0.9.0",
        "pymc3 >= 3.1",
        "numpy >= 1.13.1",
        "scipy >= 0.19.1",
        "tqdm >= 4.15.0"
    ],
)