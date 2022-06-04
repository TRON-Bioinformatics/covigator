from setuptools import find_packages, setup
import covigator


# parses requirements from file
with open("requirements.txt") as f:
    required = f.read().splitlines()

with open("README.md", "r") as f:
    long_description = f.read()

# Build the Python package
setup(
    name="covigator",
    version=covigator.VERSION,
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "covigator-ena-accessor=covigator.command_line:ena_accessor",
            "covigator-ena-downloader=covigator.command_line:ena_downloader",
            "covigator-gisaid-accessor=covigator.command_line:gisaid_accessor",
            "covigator-processor=covigator.command_line:processor",
            "covigator-dashboard=covigator.dashboard.dashboard:main",
            "covigator-pipeline=covigator.command_line:pipeline",
            "covigator-precompute=covigator.command_line:precompute_queries",
            "covigator-cooccurrence-matrix=covigator.command_line:cooccurrence"
        ],
    },
    author_email="patrick.sorn@tron-mainz.de",
    author="TRON - Translational Oncology at the University Medical Center of the Johannes Gutenberg University Mainz "
           "- Computational Medicine group",
    description="",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    requires=[],
    # NOTE: always specify versions to ensure build reproducibility
    # NOTE2: sklearn==0.19.0 is a hidden dependency as it is required by Classifier.pickle
    install_requires=required,
    setup_requires=[],
    classifiers=[
        "Development Status :: 3 - Alpha",  # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        "Intended Audience :: Healthcare Industry",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3.8",
    ],
)
