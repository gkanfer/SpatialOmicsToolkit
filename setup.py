
from setuptools import setup, find_packages

setup(
    name='SPATIALOMICSTOOLKIT',
    version='0.0.0.0.0.1',
    packages=find_packages(include=['utils','utils.vizium','utils.vizium.viziumHD','utils.vizium.viziumHD','utils.vizium.writeAnndataFile']),
    entry_points={
        'console_scripts': [
            'run_vizium=scripts.run_vizum:main',
        ],
    },
)