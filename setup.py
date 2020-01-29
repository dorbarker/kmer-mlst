from setuptools import setup, find_packages

setup(
        author='Dillon Barker',
        author_email='dillon.barker@canada.ca',
        name='acm',
        version='0.1',
        #packages=find_packages(),
        #include_package_data=True,
        py_modules=['acm'],
        install_requires=[
            'Click',
            'biopython',
            'python-magic'
        ],
        entry_points={
            'console_scripts': [
                'acm=acm:cli'
            ]
        }
)
