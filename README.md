# Anolis carolinensis in Ogasawara Archipelago

## Preparation

1.  Install `pip` and `virtualenv` for macOS system python 2
    ```sh
    % /usr/bin/python -m ensurepip -v --user
    % ~/Library/Python/2.7/bin/pip install -U setuptools pip virtualenv
    ```

1.  Create a new virtual environment, and activate it:
    ```sh
    % ~/Library/Python/2.7/bin/virtualenv ~/.virtualenv/dadi
    % source ~/.virtualenv/dadi/bin/activate
    (dadi) % pip install -U setuptools pip
    (dadi) % pip install -U flake8 seaborn
    ```

1.  Fetch and install [dadi](https://bitbucket.org/gutenkunstlab/dadi):
    ```sh
    (dadi) % git clone https://heavywatal@bitbucket.org/gutenkunstlab/dadi.git
    (dadi) % cd dadi/
    (dadi) % python setup.py install
    ```

## Procedure

1. Convert [SNP data](https://bitbucket.org/gutenkunstlab/dadi/wiki/DataFormats) to `fs` format:
   `convert2fs.py ***.txt`
1. Calculate likelihood on grid:
   `run_dadi.py -e -b 60 ***.fs`
1. Extract the MLE parameters:
   `headjson.py ***.json`
1. Try dadi optiomization from there:
   `run_dadi.py -o -l pexh-***.json  ***.fs`
1. Plot likelihood landscape: `loglik.R`
1. Make simulated samples with the MLE parameters:
   `run_ms.py -r100000 popt-***.json | gzip >msout-***.txt.gz`
1. Calculate summary statistics:
   `gunzip -c msout-***.txt.gz | sample_stats++ | gzip >stats-***.tsv.gz`
1. Summarize and plot null distributions: `nulldist.R`
