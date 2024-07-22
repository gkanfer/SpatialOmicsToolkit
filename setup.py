
from setuptools import setup, find_packages

setup(
    name='SPATIALOMICSTOOLKIT',
    version='0.0.0.0.0.2',
    packages=find_packages(include=['spatialomicstoolkit','spatialomicstoolkit.analysis','spatialomicstoolkit.qc','spatialomicstoolkit.report']),
    entry_points={
        'console_scripts': [
            'run_vizium=scripts.run_vizum:main',
        ],
    },)
