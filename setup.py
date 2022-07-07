try:
    import pybind11
except ImportError:
    # Install pybind11 if building package without isolation
    import subprocess
    import sys
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pybind11"])
    import pybind11

from setuptools import setup, Extension, find_packages

with open("README.md", "r") as readme_file:
    readme = readme_file.read()

requirements = ['pybind11', 'pandas', 'numpy', 'scipy', 'tqdm']

__version__ = "0.1.0"


ext_modules = [
    Extension("dcona.core.extern.utils",
        ["dcona/native/src/utils/utils.cpp",
         "dcona/native/src/lib/utils.cpp"],
         include_dirs=[pybind11.get_include()],
         language='c++',
         extra_compile_args=['-std=c++11', '-O3', '-Wall']
    ),

    Extension("dcona.core.extern.correlations",
        ["dcona/native/src/utils/utils.cpp",
         "dcona/native/src/lib/correlations.cpp",
         "dcona/native/src/correlations/correlations.cpp"],
         include_dirs=[pybind11.get_include()],
         language='c++',
         extra_compile_args=['-std=c++11', '-O3', '-Wall']
    ),
        
    Extension("dcona.core.extern.tests",
        ["dcona/native/src/utils/utils.cpp",
         "dcona/native/src/correlations/correlations.cpp",
         "dcona/native/src/tests/tests.cpp",
         "dcona/native/src/lib/tests.cpp"],
         include_dirs=[pybind11.get_include()],
         language='c++',
         extra_compile_args=['-std=c++11', '-O3', '-Wall']
    ),
        
    Extension("dcona.core.extern.scores",
        ["dcona/native/src/utils/utils.cpp",
         "dcona/native/src/scores/scores.cpp",
         "dcona/native/src/lib/scores.cpp"],
         include_dirs=[pybind11.get_include()],
         language='c++',
         extra_compile_args=['-std=c++11', '-O3', '-Wall']
    ),
        
    Extension("dcona.core.extern.pipelines",
        ["dcona/native/src/utils/utils.cpp",
         "dcona/native/src/correlations/correlations.cpp",
         "dcona/native/src/tests/tests.cpp",
         "dcona/native/src/scores/scores.cpp",
         "dcona/native/src/pipelines/pipelines.cpp",
         "dcona/native/src/lib/pipelines.cpp"],
         include_dirs=[pybind11.get_include()],
         language='c++',
         extra_compile_args=['-std=c++11', '-O3', '-Wall']
    ),
]


setup(
    name="dcona",
    version=__version__,
    author="Anton Zhiyanov, Narek Engibaryan",
    author_email="zhiyanovap@gmail.com",
    url="https://github.com/zhiyanov/DCoNA",
    description="Differential Correlation Network Analysis",
    long_description=readme,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=requirements,
    ext_modules=ext_modules,
    python_requires=">=3.7",
    entry_points={
        'console_scripts': [
            'dcona = dcona.__main__:main',
        ],
    },
)
