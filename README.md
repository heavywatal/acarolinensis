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
