# CoVigator (Corona Navigator)

CoVigator is a fully automatized Corona Analysis Pipeline integrating the following major steps:

* Generation/Download of Reference DBs
* Starting of different workflows depending of input data type (e.g. RNA-Seq, Assembly)

# Configuration

`config.py` contains all the essential configuration paths

# Usage

```
python processing.py -i <path_to_input_fastqs> -w <working_dir>
```

whereas <input_paths> should be a list of FASTQs or a path where the FASTQs are stored and
<working_dir> describes the path to store the results in.

## Developer guide

Build the binary: `python setup.py bdist_wheel`
Install covigator: `pip install dist/covigator-x.y.z-py3-none-any.whl`
After installation there will be two endpoints available in the path: `covigator-download` and `covigatpr-pipeline`